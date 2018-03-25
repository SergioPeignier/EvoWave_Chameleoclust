#-*- coding: utf-8 -*-
import os
import numpy as np
import re
import json
import sys
import math
import difflib
import argparse
from datetime import datetime as datetime
import csv
import StringIO

import WifiObjects as wo
from WifiObjects import WifiObjectRawTypeError

class WifiLogFile:
    """Base class for Log files manipulations"""
    
    # Expose different WifiLogFile types
    class Type:
        Log, Raw, Xml, Csv = range(4)
    
    # Specific format
    class Format:
        DateTime = "%Y/%m/%d %H:%M:%S"
        Date = "%Y/%m/%d"
    
    # Version
    VERSION = "2.0"
    
    # Precision
    class Precision:
        Floats      = 4
        Distances   = 2
        Stats       = 3
    
    # Verbosity Level
    class VerbosityLevel:
        No, Error, Warning, Debug = range(4)
    
    # Wifi Objects categories
    CAT_CONTAINER = 0
    CAT_ACQUISITION = 1
    CAT_CAPTURE = 2
    CAT_SOURCE = 3
    CAT_AUTHOR = 4
    CAT_LABEL = 5
    
    CAT_FROM_TYPE = {
        wo.WifiAcquisitionsContainer:CAT_CONTAINER,
        wo.WifiRawAcquisitionsContainer:CAT_CONTAINER,
        wo.WifiAcquisition:CAT_ACQUISITION,
        wo.WifiRawAcquisition:CAT_ACQUISITION,
        wo.WifiCapture:CAT_CAPTURE,
        wo.WifiRawCapture:CAT_CAPTURE,
        wo.WifiSource:CAT_SOURCE,
        wo.WifiRawSource:CAT_SOURCE,
        wo.WifiAcquisitionAuthor:CAT_AUTHOR,
        wo.WifiAcquisitionLabel:CAT_LABEL
    }
    
    def __init__(self, logFilePath, verbosityLevel = VerbosityLevel.No):
        self.path = logFilePath
        self.verbosityLevel = verbosityLevel
        
    def _debug(self, message):
        if self.verbosityLevel >= self.VerbosityLevel.Debug:
            print "DBG: {}".format(message)
        
    def _error(self, message):
        if self.verbosityLevel >= self.VerbosityLevel.Error:
            print "ERR: {}".format(message)
            
    def _warning(self, message):
        if self.verbosityLevel >= self.VerbosityLevel.Warning:
            print "WARN: {}".format(message)

class WifiLogFileLoader(WifiLogFile):
    """Loader of WifiLogFiles (open files and store their contents in WifiObjects)"""
        
    # Regex Patterns
    PATTERN_ACQ = re.compile(r"(?P<mode>(?:RAW|LOG))ACQ{(?P<acquisition_header>.*)}", re.IGNORECASE)
    PATTERN_CAP = re.compile(r"(?P<mode>(?:RAW|LOG))CAP{(?P<capture_header>.*)}", re.IGNORECASE)
    PATTERN_SOURCE = re.compile(r"(?P<source_macAddress>(?:[a-f0-9]{2}:){5}[a-f0-9]{2})\|(?P<source_distancesAvg>(?:[0-9]+\.[0-9]+)|None)\|(?P<source_packets>(?:[0-9]+)|None)", re.IGNORECASE)
    PATTERN_RAWSOURCE = re.compile(r"(?P<source_macAddress>(?:[a-f0-9]{2}:){5}[a-f0-9]{2})\|\[(?P<source_distances>(?:[0-9]+\.[0-9]+(?:, )?)*)\]", re.IGNORECASE)
    
    # Default Values
    DEFAULT_VALUES_ACQ = {"authorName":None, "authorMac":None, "building":None, "floor":int(), "room":None, "colleagues":list(), "comment":"", "duration":None, "spacing":None, "minDistance":None, "medDistance":None, "maxDistance":None, "distancesAvg":None, "distancesStd":None, "version":None}
    DEFAULT_VALUES_CAP = {"start":None, "minDistance":None, "medDistance":None, "maxDistance":None, "distancesAvg":None, "distancesStd":None}
    
    # Load Type
    TYPE_AUTO = 0
    TYPE_LOG = 1
    TYPE_RAW = 2
    
    def __init__(self, logFilePath, verbosityLevel = WifiLogFile.VerbosityLevel.Warning):
        WifiLogFile.__init__(self, logFilePath, verbosityLevel)
        self._acquisitions = None
        
    def _loadDictionary(self, defaultValues, definitionStr):
        loadedValues = dict()
        
        try:
            loadedValues = json.loads("{{ {} }}".format(definitionStr), encoding=sys.getdefaultencoding())
        except ValueError:
            self._error("Can't load dictionary")
            
        return dict(defaultValues.items() + loadedValues.items())
    
    def _loadList(self, listJSONStr):
        loadedValues = list()
        
        try:
            loadedValues = json.loads("[ {} ]".format(listJSONStr), encoding=sys.getdefaultencoding())
        except ValueError:
            self._error("Can't load list")
            
        return loadedValues
        
    def _isRawObject(self, mode):
        return mode.lower() == "raw"
        
    def _loadAcquisition(self, line):
        """Loading method expecting to parse an acquisition line"""
        result = re.search(self.PATTERN_ACQ, line)
        
        # WifiAcquisition found : parsing definition and entering in acquisition bloc => next line should be a capture
        if result:
            self._debug("Entering {}ACQ...".format(result.group("mode")))
            
            values = self._loadDictionary(self.DEFAULT_VALUES_ACQ, result.group("acquisition_header"))
            
            if values["version"]:
                if float(values["version"]) > float(self.VERSION):
                    self._debug("<ACQUISITION IGNORED DUE TO TOO RECENT VERSION>")
                    return self.CAT_ACQUISITION
                if float(values["version"]) < float(self.VERSION):
                    self._warning("Caution, acquisition version is older ({} < {}).".format(float(values["version"]), self.VERSION));
            
            wifiAuthor = wo.WifiAcquisitionAuthor(values["authorName"], values["authorMac"])
            wifiLabel = wo.WifiAcquisitionLabel(values["building"], values["floor"], values["room"], values["colleagues"], values["comment"])
            wifiStats = wo.WifiStats(values["minDistance"], values["medDistance"], values["maxDistance"], values["distancesAvg"], values["distancesStd"])
            
            if not self._isRawObject(result.group("mode")):
                if not self._acquisitions:
                    self._acquisitions = wo.WifiAcquisitionsContainer()
                wifiAcq = wo.WifiAcquisition(values["duration"], values["spacing"], wifiAuthor, wifiLabel, wifiStats)
            else:
                if not self._acquisitions:
                    self._acquisitions = wo.WifiRawAcquisitionsContainer()
                wifiAcq = wo.WifiRawAcquisition(values["duration"], values["spacing"], wifiAuthor, wifiLabel, wifiStats)
            
            self._acquisitions.addAcquisition(wifiAcq)
            
            return self.CAT_CAPTURE
        
        # WifiAcquisition not found : ignoring line
        self._debug("<IGNORED LINE>")
        return self.CAT_ACQUISITION
        
    def _loadCapture(self, line):
        """Loading method expecting to parse a capture line"""
        result = re.search(self.PATTERN_CAP, line)
        
        # WifiCapture found : parsing definition and entering in capture bloc => next line should be a source
        if result:
            self._debug("Entering {}CAP...".format(result.group("mode")))
            
            values = self._loadDictionary(self.DEFAULT_VALUES_CAP, result.group("capture_header"))
            
            try:
                start = datetime.strptime(str(values["start"]), self.Format.DateTime)
            except ValueError:
                self._debug("Can't convert capture start datetime ({})".format(str(values["start"])))
                start = self.DEFAULT_VALUES_CAP["start"]
            
            wifiStats = wo.WifiStats(values["minDistance"], values["medDistance"], values["maxDistance"], values["distancesAvg"], values["distancesStd"])
            
            if not self._isRawObject(result.group("mode")):
                wifiCapture = wo.WifiCapture(start, wifiStats)
            else:
                wifiCapture = wo.WifiRawCapture(start, wifiStats)
            
            self._acquisitions.getLastAcquisition().addCapture(wifiCapture)
            
            return self.CAT_SOURCE
        
        # WifiCapture not found : is this an acquisition ?
        return self._loadAcquisition(line)
    
    def _loadSource(self, line):
        """Loading method expecting to parse a source line"""
        
        # WifiSource found : parsing line
        result = re.search(self.PATTERN_SOURCE, line)
        if result:
            self._debug("Reading SRC. ({})".format(result.group("source_macAddress")))

            try:
                packets = int(result.group("source_packets"))
            except ValueError:
                self._debug("Can't get source packets' count ({})".format(result.group("source_packets")))
                packets = None
            
            wifiStats = wo.WifiStats(distancesAvg=result.group("source_distancesAvg"))
            wifiSource = wo.WifiSource(result.group("source_macAddress"), packets, wifiStats)
            self._acquisitions.getLastAcquisition().getLastCapture().addSource(wifiSource)
            
            return self.CAT_SOURCE
        
        # WifiRawSource found : parsing line
        result = re.search(self.PATTERN_RAWSOURCE, line)
        if result:
            self._debug("Reading RAWSRC. ({})".format(result.group("source_macAddress")))
            
            distances = self._loadList(result.group("source_distances"))
            wifiSource = wo.WifiRawSource(result.group("source_macAddress"), distances)
            self._acquisitions.getLastAcquisition().getLastCapture().addSource(wifiSource)
            
            return self.CAT_SOURCE
        
        # WifiSource not found : is this an capture ?
        return self._loadCapture(line)
        
    # Map load function's references in order to apply them as a switch-statement 
    _LOAD_BLOC = {
        WifiLogFile.CAT_ACQUISITION    : _loadAcquisition,
        WifiLogFile.CAT_CAPTURE        : _loadCapture,
        WifiLogFile.CAT_SOURCE         : _loadSource
    }
        
    def load(self, loadType=TYPE_AUTO):
        """Load files at path and return a corresponding WifiAcquisitionContainer"""
        if os.path.exists(self.path):
            if not os.path.isfile(self.path):
                raise IOError("""Given path isn't a file ("{}")""".format(self.path))
        else:
            raise IOError("""Given path doesn't exist ("{}")""".format(self.path))
        
        self._debug("Loading file '{}'".format(self.path))
        
        currentBloc = self.CAT_ACQUISITION
        
        if loadType == self.TYPE_LOG:
            self._acquisitions = wo.WifiAcquisitionsContainer()
        elif loadType == self.TYPE_RAW:
            self._acquisitions = wo.WifiRawAcquisitionsContainer()
        else:
            self._acquisitions = None
        
        with open(self.path, "r") as logFile:
            for line in logFile:
                currentBloc = self._LOAD_BLOC[currentBloc](self, line)
        
        return self._acquisitions

class WifiLogFileDumper(WifiLogFile):
    """Dumper of WifiLogFiles (saves WifiObjects into a file)"""
    
    # Dump file conflict resolution
    class ConflictResolution:
        Error, Append, Override = range(3)
        
        @classmethod
        def getFileMode(cls, conflictResolution):
            return "a" if conflictResolution == cls.Append else "w"
    
    # Dump tags' names
    _DUMPTAG_NAMES = {
        WifiLogFile.Type.Log: {
            WifiLogFile.CAT_ACQUISITION : "LOGACQ",
            WifiLogFile.CAT_CAPTURE : "LOGCAP"
        },
        WifiLogFile.Type.Raw: {
            WifiLogFile.CAT_ACQUISITION : "RAWACQ",
            WifiLogFile.CAT_CAPTURE : "RAWCAP"
        }
    }
    
    class XmlOptions:
        """Options available in XML mode only.
        'indentSymbol' is symbol used for indentation (set to "" in order to avoid indentation) 
        Set 'displayNone' to true if you want to display attributes with None value."""
        def __init__(self, displayNone=False, indentSymbol='\t'):
            self.indentSymbol = indentSymbol
            self.displayNone = displayNone
            
    class CsvOptions:
        def __init__(self, delimiter=',', escapechar="\\", doublequote=False, quoting=csv.QUOTE_MINIMAL):
            self.delimiter = delimiter
            self.escapechar = escapechar
            self.doublequote = doublequote
            self.quoting = quoting
    
    def __init__(self, logFilePath=None, verbosityLevel=WifiLogFile.VerbosityLevel.Warning, xmlOptions=XmlOptions(), csvOptions=CsvOptions(), macKeyLen=2):
        """WifiObjects will be saved at 'logFilePath'"""
        WifiLogFile.__init__(self, logFilePath, verbosityLevel)
        self._dumpMode = None
        self._xmlOptions = xmlOptions
        self._csvOptions = csvOptions
        self._writer = None
        #default options
        self.macKeyLen = macKeyLen
    
    def _xmlTag(self, indent=0, tagName="TAG", attributes=dict(), content=str()):
        attributesContent = str()
        for key in attributes.keys():
            if self._xmlOptions.displayNone or attributes[key]:
                attributesContent += ' {}="{}"'.format(key, attributes[key])

        formatDict = {
            "indent":(self._xmlOptions.indentSymbol * indent),
            "tagName":tagName,
            "attributes":attributesContent,
            "content":content
        }
    
        if content:
            return """{indent}<{tagName}{attributes}>\n{content}{indent}</{tagName}>\n""".format(**formatDict)
        else:
            return """{indent}<{tagName}{attributes}/>\n""".format(**formatDict)
    
    def _logTag(self, tagName="TAG", attributes=dict(), content=str()):
        formatDict = {
            "tagName":tagName,
            "attributes":json.dumps(attributes, skipkeys=True, encoding=sys.getdefaultencoding(), allow_nan=False),
            "content":content
        }
        
        return """{tagName}{attributes}\n{content}""".format(**formatDict)
    
    def _formatAttribute(self, attribute):
        if type(attribute) is datetime:
            return attribute.strftime(WifiLogFile.Format.DateTime)
        
        if type(attribute) in [float, np.float64]:
            if math.isinf(attribute) or math.isnan(attribute):
                return None
            return round(attribute, self.Precision.Floats)
        
        return attribute
    
    def _formatDictAttributes(self, attributes = dict()):
        for key in attributes.keys():
            attributes[key] = self._formatAttribute(attributes[key])
        
    def _formatListAttributes(self, attributes = list()):
        for i in range(0, len(attributes)):
            attributes[i] = self._formatAttribute(attributes[i])
    
    def _dumpTag(self, tagName, indent = 0, attributes = dict(), content = str()):
        self._formatDictAttributes(attributes)
                
        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw]:
            return self._logTag(tagName, attributes, content)
        
        if self._dumpMode == WifiLogFile.Type.Xml:
            return self._xmlTag(indent, tagName, attributes, content)
        
    def _write(self, content):
        if content :
            if self._dumpMode == WifiLogFile.Type.Csv:
                self._writer.writerow(content)
            else:
                self._writer.write(content)
        
    def _dumpAcquisitionsContainer(self, acquisitionsContainer, indent = 0):
        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw]:
            content = str()
            
            for acquisition in acquisitionsContainer.acquisitions:
                content += self._dumpAcquisition(acquisition, )
                
            return content
        
        if self._dumpMode == WifiLogFile.Type.Csv:
            acqNumber = 0
            for acquisition in acquisitionsContainer.acquisitions:
                self._dumpAcquisition(acquisition, acqNumber)
                acqNumber += 1
            return None
                
        if self._dumpMode == WifiLogFile.Type.Xml:
            firstAcq = acquisitionsContainer.getFirstAcquisition(True)
            firstCapt = firstAcq.getFirstCapture(True) if firstAcq else None 
            
            lastAcq = acquisitionsContainer.getLastAcquisition(True)
            lastCapt = lastAcq.getLastCapture(True) if lastAcq else None
            
            attributes = {
                "acquisitions":len(acquisitionsContainer.acquisitions),
                "captures":acquisitionsContainer.getCapturesCount(),
                "sources":acquisitionsContainer.getSourcesCount(),
                "packets":acquisitionsContainer.getPacketsCount(),
                "firstCapture":firstCapt.start if firstCapt else None,
                "lastCapture":lastCapt.start if lastCapt else None,
            }
            attributes = dict(acquisitionsContainer.getStatistics().items() + attributes.items())
            
            content = str()
            for acquisition in acquisitionsContainer.acquisitions:
                content += self._dumpAcquisition(acquisition, indent + 1)
    
            return self._dumpTag(type(acquisitionsContainer).__name__, indent, attributes, content)
    
    def _dumpAcquisition(self, acquisition, indent=0, container=None):
        if self._dumpMode == WifiLogFile.Type.Csv:
            for capture in acquisition._captures:
                self._dumpCapture(capture, indent, acquisition)
            return None
        
        firstCapture = acquisition.getFirstCapture(True)
        lastCapture = acquisition.getLastCapture(True)
        
        attributes = {key : acquisition.__dict__[key] for key in ["duration", "spacing"]}
        attributes = dict({
            "version":self.VERSION,
            "sources":acquisition.getSourcesCount(),
            "packets":acquisition.getPacketsCount(),
            "captures":acquisition.getCapturesCount(),
            "firstCapture":firstCapture.start if firstCapture else None,
            "lastCapture":lastCapture.start if lastCapture else None
        }.items() + attributes.items() + acquisition.stats.getStatistics(self.Precision.Stats, self.Precision.Distances).items())
        
        content = str()
        for capture in acquisition._captures:
            content += self._dumpCapture(capture, indent + 1)
                
        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw]:
            attributes = dict(self._dumpAuthor(acquisition.author).items() + attributes.items())
            attributes = dict(self._dumpLabel(acquisition.label).items() + attributes.items())
            
            return self._dumpTag(self._DUMPTAG_NAMES[self._dumpMode][self.CAT_FROM_TYPE[type(acquisition)]], indent, attributes, content)
        
        if self._dumpMode == WifiLogFile.Type.Xml:
            content = self._dumpAuthor(acquisition.author, indent + 1) + self._dumpLabel(acquisition.label, indent + 1) + content

            return self._dumpTag(type(acquisition).__name__, indent, attributes, content)
        
    def _dumpCapture(self, capture, indent=0, container=None):      
        attributes = dict({"start":capture.start, "sources":capture.getSourcesCount(), "packets":capture.getPacketsCount()}.items() + capture.stats.getStatistics(self.Precision.Stats, self.Precision.Distances).items())

        if self._dumpMode == WifiLogFile.Type.Csv:
            rows = {key:None for key in self._writer.fieldnames}
            rows = dict(rows.items() + attributes.items()) #rows.update(attributes) ?
            rows = dict(rows.items() + self._dumpAuthor(container.author).items())
            rows = dict(rows.items() + self._dumpLabel(container.label).items())
            rows = dict(rows.items() + {
                "acquisition":indent,
                "duration":container.duration,
                "captureUID":int((capture.start - datetime(2015,1,1)).total_seconds())
            }.items())
            
            macConflicts = dict()
            for source in capture._sources:
                mac, avg = self._dumpSource(source, indent + 1)
                if not np.isnan(avg):
                    key = self._macToKey(mac)
                
                    if rows[key]:
                        macConflicts[key] += 1
                        rows[key] += round(avg, self.Precision.Stats)
                    else:
                        macConflicts[key] = 1
                        rows[key] = round(avg, self.Precision.Stats)
                    
            for key in macConflicts.keys():
                rows[key] = round(rows[key] / float(macConflicts[key]), self.Precision.Stats)

            self._write(rows)
            return None
        
        content = str()
        for source in capture._sources:
            content += self._dumpSource(source, indent + 1)

        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw]:
            return self._dumpTag(self._DUMPTAG_NAMES[self._dumpMode][self.CAT_FROM_TYPE[type(capture)]], indent, attributes, content)
        
        if self._dumpMode == WifiLogFile.Type.Xml:
            return self._dumpTag(type(capture).__name__, indent, attributes, content)
        
    def _dumpSource(self, source, indent=0, container=None):  
        if self._dumpMode == WifiLogFile.Type.Log:
            stats = source.stats.getStatistics(self.Precision.Stats, self.Precision.Distances)
            values = [source.macAddress, stats["distancesAvg"], source.packets] # Values will be printed in same order as those keys
            self._formatListAttributes(values)
            return "|".join("{}".format(k) for k in values) + "\n"
        
        if self._dumpMode == WifiLogFile.Type.Raw:
            values = [source.macAddress, map(lambda x: round(x, self.Precision.Distances), source._distances) if hasattr(source, "_distances") else list()] # Values will be printed in same order as those keys
            self._formatListAttributes(values)
            return "|".join("{}".format(k) for k in values) + "\n"
        
        if self._dumpMode == WifiLogFile.Type.Xml:
            attributes = {key : source.__dict__[key] for key in ["packets", "macAddress"]}
            attributes = dict(source.stats.getStatistics(self.Precision.Stats, self.Precision.Distances).items() + attributes.items())
            return self._dumpTag(type(source).__name__, indent, attributes)
        
        if self._dumpMode == WifiLogFile.Type.Csv:
            return (source.macAddress, source.stats.getStatistics(self.Precision.Stats, self.Precision.Distances)["medDistance"])
    
    def _dumpAuthor(self, author, indent=0, container=None):
        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw, WifiLogFile.Type.Csv]:
            return {
                "authorName":author.name,
                "authorMac":author.macAddress
            }
            
        if self._dumpMode == WifiLogFile.Type.Xml:
            attributes = {key : author.__dict__[key] for key in ["name", "macAddress"]}
            return self._dumpTag(type(author).__name__, indent, attributes)
    
    def _dumpLabel(self, label, indent=0, container=None):
        attributes = {key : label.__dict__[key] for key in ["building", "room", "floor", "comment"]}
        
        if self._dumpMode in [WifiLogFile.Type.Log, WifiLogFile.Type.Raw]:
            attributes["colleagues"] = label.colleagues
            return attributes
        
        if self._dumpMode == WifiLogFile.Type.Csv:
            attributes["colleagues"] = ";".join(label.colleagues) if label.colleagues else ""
            return attributes
        
        if self._dumpMode == WifiLogFile.Type.Xml:
            content = str()
            if label.colleagues:
                for colleague in label.colleagues:
                    content += self._dumpTag("WifiAcquisitionColleague", indent + 1, {"name":colleague})
                    
            return self._dumpTag(type(label).__name__, indent, attributes, content)
        
    # Map dump function's references in order to apply them as a switch-statement
    _DUMPERS = {
        WifiLogFile.CAT_CONTAINER  : _dumpAcquisitionsContainer,
        WifiLogFile.CAT_ACQUISITION : _dumpAcquisition,
        WifiLogFile.CAT_CAPTURE : _dumpCapture,
        WifiLogFile.CAT_SOURCE : _dumpSource,
        WifiLogFile.CAT_AUTHOR : _dumpAuthor,
        WifiLogFile.CAT_LABEL : _dumpLabel
    }
    
    def _macToKey(self, mac):
        mac = ''.join(mac.split(':'))
        return "MAC-" + str(int(mac[-self.macKeyLen:], 16))
        # key = "MAC-" + str(int(mac[-2:], 16))
        # if mac[-4]%2==0:
        #     key += 'p' 
        # return key
        # return "MAC-" + str(int(mac[-2:], 16))
        
    def _initWriter(self, logStream):
        if self._dumpMode == WifiLogFile.Type.Csv:
            fieldnames = ["captureUID", "acquisition", "start", "duration", "authorName", "authorMac", "building", "floor", "room", "colleagues", "comment", "packets", "sources", "distancesAvg", "distancesStd", "maxDistance", "medDistance", "minDistance"]
            nb_macKeys = 16**self.macKeyLen
            fieldnames += map(lambda x:"MAC-"+str(x), range(nb_macKeys))
            self._writer = csv.DictWriter(logStream, fieldnames=fieldnames, delimiter=self._csvOptions.delimiter, escapechar=self._csvOptions.escapechar, doublequote=self._csvOptions.doublequote, quoting=self._csvOptions.quoting)
            self._writer.writeheader()
        else:
            self._writer = logStream
        
    def dump(self, wifiObject, fileResolution=ConflictResolution.Error, mode=WifiLogFile.Type.Log):
        self._dumpMode = mode
        self._debug("Dumping file '{}'".format(self.path))

        if mode == WifiLogFile.Type.Raw and not wifiObject.isWifiRawObject():
            raise WifiObjectRawTypeError("WifiRawObject", type(wifiObject))
        
        if self.path:
            if os.path.isfile(self.path) and fileResolution == self.ConflictResolution.Error:
                raise IOError("""Given path already exists ("{}"). Use 'fileResolution' parameter to specify if you to overwrite it or append it.""".format(self.path))
        
            directory = os.path.dirname(self.path);
            if not os.path.isdir(directory):
                os.makedirs(directory);
        
        logStream = StringIO.StringIO()
        self._initWriter(logStream)

        # Call right specific dumper based on given wifiObject type
        try:
            method = self._DUMPERS[self.CAT_FROM_TYPE[type(wifiObject)]]
        except KeyError:
            raise TypeError("WifiLogFileDumper.dump() requires valid type for parameter wifiObject. (Given type : '{}')".format(type(wifiObject)))
        
        self._write(method(self, wifiObject))
        
        if self.path:
            with open(self.path, self.ConflictResolution.getFileMode(fileResolution)) as logFile:
                logFile.write(logStream.getvalue())
        else:
            return logStream.getvalue()

class Test:
    def __init__(self, printDebug = False, catchExceptions = True, logFileVerbosityLevel=WifiLogFile.VerbosityLevel.Warning):
        self.printDebug = printDebug
        self.catchExceptions = catchExceptions
        self.logFileVerbosityLevel = logFileVerbosityLevel
        
    def _getBaseFileName(self, fileName):
        return '.'.join(fileName.split('.')[:-1])
    
    def _compareFiles(self, fileNameA, fileNameB):
        self._debug("""Comparing "{}" and "{}"...""".format(fileNameA, fileNameB))
        
        if os.path.getsize(fileNameA) == os.path.getsize(fileNameB) and open(fileNameA,'rb').read() == open(fileNameB,'rb').read():
            self._ok()
            return True
            
        else:
            diff = difflib.ndiff(open(fileNameA).readlines(), open(fileNameB).readlines())
            changes = [l for l in diff if l.startswith('+ ') or l.startswith('- ') or l.startswith('? ')]
            
            differences = "\n"
            for change in changes:
                differences += "\t\t" + change
            
            self._fail(differences + "\t")
            return False
        
    def _debug(self, message):
        if self.printDebug:
            print "\t{}".format(message),
        
    def _fail(self, exception):
        if self.printDebug:
            print """FAIL ("{}")""".format(exception)
        
    def _ok(self):
        if self.printDebug:
            print "OK"
        
    def _try(self, function, kwargs = dict()):
        try:
            result = function(**kwargs)
        except Exception, e:
            if self.catchExceptions:
                self._fail(e)
            else:
                raise
        else:
            self._ok()
            return result if result else True
        
        return False
    
    def _printTestResult(self, testName, success):
        if self.printDebug:
            print """Test "{}" : {}""".format(testName, "SUCCESS" if success else "FAIL")
        else:
            print "SUCCESS" if success else "FAIL"
        
    def _printTest(self, message):
        if self.printDebug:
            print """Testing "{}"...""".format(message)
        else:
            print """Testing "{}"...""".format(message),
        
    def _loadFile(self, fileName):
        self._debug("""Loading "{}"...""".format(fileName))
        return self._try(WifiLogFileLoader(fileName, verbosityLevel=self.logFileVerbosityLevel).load)
    
    def _printFile(self, container, fileName = ""):
        self._debug("""Printing "{}"...""".format(fileName) if fileName else "Printing...")
        content = self._try(WifiLogFileDumper(xmlOptions=WifiLogFileDumper.XmlOptions(True), verbosityLevel=self.logFileVerbosityLevel).dump, {"wifiObject":container, "mode":WifiLogFile.Type.Xml})
        
        if content :
            print content
            return True
        
        return False
    
    def _dumpFile(self, container, fileName, dumpMode = WifiLogFile.Type.Log):
        self._debug("""Dumping "{}"...""".format(fileName))
        return self._try(WifiLogFileDumper(fileName, verbosityLevel=self.logFileVerbosityLevel).dump, {"wifiObject":container, "fileResolution":WifiLogFileDumper.ConflictResolution.Override, "mode":dumpMode})
    
    def _extractObject(self, wifiObject, kwargs = dict()):
        self._debug("Extracting...")
        return self._try(wifiObject.extract, kwargs)
    
    def loadNPrint(self, fileName):
        self._printTest("loadNPrint")
        result = False
        container = self._loadFile(fileName)
        
        if container:
            result = self._printFile(container, fileName)
        
        self._printTestResult("loadNPrint", result)
        return result
        
    def loadNDump(self, fileName):
        self._printTest("loadNDump")
        result = False
        container = self._loadFile(fileName)
        
        if container:
            baseName = self._getBaseFileName(fileName)
            outLog = baseName + ".out.log"
            
            result = self._dumpFile(container, outLog)
            
            if container.isWifiRawObject():
                outRawLog = baseName + ".out.rawlog"
                result &= self._dumpFile(container, outRawLog, WifiLogFile.Type.Raw)
            
        self._printTestResult("loadNDump", result)
        return result
        
    def reload(self, fileName):
        self._printTest("reload")
        result = False
        originContainer = self._loadFile(fileName)
        
        if originContainer:
            fileNames = "{}.out{{}}.{{}}log".format(self._getBaseFileName(fileName))
            outLog1 = str.format(fileNames, "1", "")
            outRawLog1 = str.format(fileNames, "1", "raw")
            outLog2 = str.format(fileNames, "2", "")
            outRawLog2 = str.format(fileNames, "2", "raw")
        
            logResult = False
            if self._dumpFile(originContainer, outLog1):
                outContainer = self._loadFile(outLog1)
                logResult = self._dumpFile(outContainer, outLog2) and self._compareFiles(outLog1, outLog2)
            
            rawResult = True
            if originContainer.isWifiRawObject():
                rawResult = False
                if self._dumpFile(originContainer, outRawLog1, dumpMode=WifiLogFile.Type.Raw):
                    outContainer = self._loadFile(outRawLog1)
                    
                    if outContainer:
                        rawResult = self._dumpFile(outContainer, outRawLog2, dumpMode=WifiLogFile.Type.Raw) and self._compareFiles(outRawLog1, outRawLog2)
        
            result = logResult and rawResult
            
        self._printTestResult("reload", result)
        return result
        
    def extractAcquisition(self, fileName):
        self._printTest("extractAcquisition")
        result = False
        
        container = self._loadFile(fileName)
        
        if container:
            acquisition = container.getFirstAcquisition()
            
            if acquisition:
                extraction = self._extractObject(acquisition, {"minSources":5})
                
                if extraction:
                    self._debug("Comparing extraction and original...")
                    if extraction.getCapturesCount() < acquisition.getCapturesCount():
                        self._ok()
                        result = True
                    else:
                        self._fail("{} < {}".format(extraction.getCapturesCount(), acquisition.getCapturesCount()))
            
        self._printTestResult("extractAcquisition", result)
        return result
        
    def extractCapture(self, fileName):
        self._printTest("extractCapture")
        result = False
        
        container = self._loadFile(fileName)
        
        if container:
            capture = container.getFirstAcquisition().getFirstCapture()
            
            if capture:
                extraction = self._extractObject(capture, {"minPackets":5})
                
                if extraction:
                    self._debug("Comparing extraction and original...")
                    if extraction.getSourcesCount() < capture.getSourcesCount():
                        self._ok()
                        result = True
                    else:
                        self._fail("{} < {}".format(extraction.getSourcesCount(), capture.getSourcesCount()))
            
        self._printTestResult("extractCapture", result)
        return result
    
    def main(self, fileName):
        self.loadNPrint(fileName)
        self.loadNDump(fileName)
        self.extractAcquisition(fileName)
        self.extractCapture(fileName)
        self.reload(fileName)
        
class WifiLog(object):
    """WifiLog command line interface"""
    
    VERSION = "2.0"
    
    def __init__(self):
        self._cmdLineArgParser = self._initArgParser()
        self._options = self._cmdLineArgParser.parse_args()
    
    def _wifiVerbosityLevel(self):
        return WifiLogFile.VerbosityLevel.Debug if self._options.verbose else WifiLogFile.VerbosityLevel.Warning
    
    def _getDumpMode(self):
        if self._options.type in ["log"]:
            return WifiLogFile.Type.Log
        
        if self._options.type in ["raw"]:
            return WifiLogFile.Type.Raw
        
        if self._options.type in ["xml"]:
            return WifiLogFile.Type.Xml
        
        if self._options.type in ["csv"]:
            return WifiLogFile.Type.Csv
        
        return None
    
    def _getLoadType(self):
        if self._options.type in ["log"]:
            return WifiLogFileLoader.TYPE_LOG
        
        if self._options.type in ["raw"]:
            return WifiLogFileLoader.TYPE_RAW
        
        return WifiLogFileLoader.TYPE_AUTO
        
    DEFAULT_EXTENSIONS = {
        WifiLogFile.Type.Log:"log",
        WifiLogFile.Type.Raw:"rawlog",
        WifiLogFile.Type.Xml:"xml",
        WifiLogFile.Type.Csv:"csv"
    }
    
    def _convertFile(self, inFilePath):
        outFilePath = self._getOutFilePath(inFilePath, self._options.outPath)

        self._msg("""Converting file "{}" to "{}"...""".format(inFilePath, outFilePath))
        
        if not self._options.override and os.path.exists(outFilePath):
            return self._error("""Output file "{}" already exists.""".format(outFilePath));
        
        try:
            wifiContainer = WifiLogFileLoader(inFilePath, verbosityLevel=self._wifiVerbosityLevel()).load()
        except Exception, e:
            self._error("""Can't load file "{}" ({}).""".format(inFilePath, str(e)))
            raise
        
        if wifiContainer:           
            self._msg("""File "{}" successfully loaded !""".format(inFilePath))
            
            if not wifiContainer.isWifiRawObject() and self._getDumpMode() == WifiLogFile.Type.Raw:
                return self._error("""Can't convert file "{}" : it's impossible to convert WifiObjects to WifiRawObjects.""".format(inFilePath))

            try:
                print type(wifiContainer)
                WifiLogFileDumper(outFilePath, 
                                    verbosityLevel=self._wifiVerbosityLevel(),
                                    macKeyLen=self._options.macKeyLen,
                                    ).dump(wifiContainer, 
                                            fileResolution=WifiLogFileDumper.ConflictResolution.Override if self._options.override else WifiLogFileDumper.ConflictResolution.Error, 
                                            mode=self._getDumpMode(),
                                            )
            except Exception, e:
                self._error("""Can't convert file "{}" ({}).""".format(inFilePath, str(e)))
                raise
            else:
                self._msg("""File "{}" successfully converted to "{}" !""".format(inFilePath, outFilePath))
            
        else:
            self._error("""File "{}" doesn't contain any valid wifi acquisition.""".format(inFilePath))
            
    def convert(self):        
        if not os.path.exists(self._options.path):
            return self._error("""Given path doesn't exist ({}).""".format(self._options.path))
        
        if os.path.isdir(self._options.path):
            for filename in os.listdir(self._options.path):
                path = os.path.join(self._options.path, filename)
                if os.path.isfile(path):
                    self._convertFile(path)
        else:
            self._convertFile(self._options.path)
        
    def test(self):
        Test(printDebug=True, catchExceptions=False, logFileVerbosityLevel=self._wifiVerbosityLevel()).main(self._options.path)
    
    def _loadFile(self, inFilename):
        try:
            wifiContainer = WifiLogFileLoader(inFilename, verbosityLevel=self._wifiVerbosityLevel()).load(self._getLoadType())
        except wo.WifiObjectRawTypeError:
            self._msg("""File "{}": {} found instead of {}.""".format(inFilename, "WifiRawAcquisitions" if self._getLoadType() == WifiLogFileLoader.TYPE_LOG else "WifiAcquisitions", "WifiAcquisitions" if self._getLoadType() == WifiLogFileLoader.TYPE_LOG else "WifiRawAcquisitions"))
            return None
        except Exception, e:
            self._error("""File "{}": An error occurs while loading ({}).""".format(inFilename, str(e)))
            return None
        
        if wifiContainer and (self._options.keepEmpty or not wifiContainer.isEmpty()):
            self._msg("""File "{}": {} {} has been successfully loaded !""".format(inFilename, wifiContainer.getAcquisitionsCount(), "WifiRawAcquisitions" if wifiContainer.isWifiRawObject() else "WifiAcquisitions"))
            return wifiContainer
        else:
            if self._options.type == "log":
                acq = "WifiAcquisitions"
            elif self._options.type == "raw":
                acq = "WifiRawAcquisitions"
            else:
                acq = "acquisitions"
                
            self._msg("""File "{}": no {} has been found.""".format(inFilename, acq))
            return None
        
    def _getExtractArgs(self):
        kwargs = dict()
        
        if self._options.fromDate:
            try:
                kwargs["fromDatetime"] = datetime.strptime(self._options.fromDate, WifiLogFile.Format.Date)
            except ValueError, e:
                dateError = str(e)

            if "fromDatetime" not in kwargs:
                try:
                    kwargs["fromDatetime"] = datetime.strptime(self._options.fromDate, WifiLogFile.Format.DateTime)
                except ValueError, e:
                    dateTimeError = str(e)
            
            if "fromDatetime" not in kwargs:
                self._error("""Parameter "fromDate" will be ignored (DATE:"{}" - DATETIME:"{}")""".format(dateError, dateTimeError))
        
        if self._options.toDate:
            try:
                kwargs["toDatetime"] = datetime.strptime(self._options.toDate, WifiLogFile.Format.Date)
            except ValueError, e:
                dateError = str(e)

            if "toDatetime" not in kwargs:
                try:
                    kwargs["toDatetime"] = datetime.strptime(self._options.toDate, WifiLogFile.Format.DateTime)
                except ValueError, e:
                    dateTimeError = str(e)
            
            if "toDatetime" not in kwargs:
                self._error("""Parameter "toDate" will be ignored (DATE:"{}" - DATETIME:"{}")""".format(dateError, dateTimeError))
                
        if self._options.allowedMac or self._options.forbiddenMac:
            kwargs["macList"] = list()
            
            if self._options.allowedMac:
                macListFileName = self._options.allowedMac
                kwargs["macForbidden"] = False
            else:
                macListFileName = self._options.forbiddenMac
                kwargs["macForbidden"] = True
            
            try:
                with open(macListFileName, "r") as macListFile:
                    kwargs["macList"].append(macListFile.readlines())
            except IOError, e:
                self._error("""Parameter "macList" will be ignored ({})""".format(e))

        kwargs["author"] = wo.WifiAcquisitionAuthor(self._options.authorName, self._options.authorMac)
        kwargs["label"] = wo.WifiAcquisitionLabel(self._options.building, self._options.floor, self._options.room, self._options.colleagues, self._options.comment)

        kwargs["minSpacing"] = self._options.minSpc
        kwargs["maxSpacing"] = self._options.maxSpc
        kwargs["minDuration"] = self._options.minDur
        kwargs["maxDuration"] = self._options.maxDur

        kwargs["minCaptures"] = self._options.minCap
        kwargs["maxCaptures"] = self._options.maxCap
        kwargs["minSources"] = self._options.minSrc
        kwargs["maxSources"] = self._options.maxSrc
        kwargs["minPackets"] = self._options.minPck
        kwargs["maxPackets"] = self._options.maxPck
        
        kwargs["minStats"] = wo.WifiStats(self._options.minMin, self._options.minMed, self._options.minMax, self._options.minAvg, self._options.minStd)
        kwargs["maxStats"] = wo.WifiStats(self._options.maxMin, self._options.maxMed, self._options.maxMax, self._options.maxAvg, self._options.maxStd)

        return kwargs

    def _getOutFilePath(self, inPath, outPath, outBaseFilenameAppend=""):
        ### So we want to construct output file path as following : outPathDir/outBaseFilename.outExtFilename
        
        ### 1) Get the output path directory (outPathDir)
        # Get outPath as output base path is specified
        outBasePath = outPath if outPath else inPath

        # Extract output path directory from output base path (ex: "./a/b/c.ext" => "./a/b")
        outPathDir = outBasePath if os.path.isdir(outBasePath) else os.path.dirname(outBasePath)
        
        # Set output path directory to default (".") if empty
        outPathDir = "." if not outPathDir else outPathDir
        
        ### 2) Get the output filename extension (outExtFilename) 
        # Split extension and root path from output base path (ex: "./a/b/c.ext" => ("./a/b/c", ".ext"))
        outRootPath, outExtFilename = os.path.splitext(outBasePath)
        
        # Keep current extension only if outPath is specified by user
        outExtFilename = outExtFilename[1:] if outPath else ""
    
        # Set extension to default if empty
        outExtFilename = outExtFilename if outExtFilename else self.DEFAULT_EXTENSIONS[self._getDumpMode()]
        
        ### 3) Get the output base filename (outBaseFilename)
        # Get the base file name from root path (ex: "./a/b/c.ext" => "c.ext"). Note that at this step we can assume we have no extension in outRootPath, so in outBaseFilename.
        outBaseFilename = os.path.basename(outRootPath)
        
        # Set output base filename to default ("out") if empty
        outBaseFilename = "out" if not outBaseFilename or os.path.isdir(outBasePath) else outBaseFilename
        
        # Add ".out" to output base filename if outPath is not specified by user
        outBaseFilename += outBaseFilenameAppend if not outPath and outBaseFilename != "out" else ""
        
        self._debug("outPathDir={}".format(outPathDir))
        self._debug("outBaseFilename={}".format(outBaseFilename))
        self._debug("outExtFilename={}".format(outExtFilename))
        
        ### Here we are, now we got the right output file path !
        return os.path.join(outPathDir, outBaseFilename + "." + outExtFilename)

    def extract(self):
        wifiContainer = None
        
        if not os.path.exists(self._options.path):
            return self._error("""Given path doesn't exist ({}).""".format(self._options.path))
            
        if os.path.isdir(self._options.path):
            self._msg("""Extracting acquisitions from files in "{}"...""".format(self._options.path))

            for filename in os.listdir(self._options.path):
                path = os.path.join(self._options.path, filename)
                
                if os.path.isfile(path):
                    self._debug("""Loading file {}...""".format(path))
                    
                    content = self._loadFile(path)
                    
                    if content:
                        if not wifiContainer:
                            wifiContainer = content
                            if self._options.type == "auto":
                                self._options.type = "raw" if wifiContainer.isWifiRawObject() else "log"
                        else:
                            try:
                                wifiContainer.append(content)
                            except Exception, e:
                                self._error("""Can't append acquisitions from file "{}" ({})""".format(path, e))
        else:
            self._msg("""Extracting acquisitions from file "{}"...""".format(self._options.path))

            wifiContainer = self._loadFile(self._options.path)
            if self._options.type == "auto":
                self._options.type = "raw" if wifiContainer.isWifiRawObject() else "log"

            
        if wifiContainer:
            self._msg("""{} WifiAcquisitions has been successfully loaded !""".format(len(wifiContainer.acquisitions)))

            outPath = self._getOutFilePath(self._options.path, self._options.outFile, ".out");
            
            if not self._options.keepDuplicate:
                wifiContainer.removeDuplicates()
            
            extraction = wifiContainer.extract(**self._getExtractArgs())
                
            if extraction:
                self._msg("""{} WifiAcquisitions has been successfully extracted !""".format(len(extraction.acquisitions)))
                
                if extraction.hasAcquisitions():
                    fileResolution = WifiLogFileDumper.ConflictResolution.Override if self._options.override else WifiLogFileDumper.ConflictResolution.Error
                    fileResolution = WifiLogFileDumper.ConflictResolution.Append if self._options.append else fileResolution
                    
                    try:
                        WifiLogFileDumper(outPath, verbosityLevel=self._wifiVerbosityLevel()).dump(extraction, fileResolution=fileResolution, mode=self._getDumpMode())
                    except Exception, e:
                        return self._error("""Can't save extracted WifiAcquisitions in file "{}" ({}).""".format(outPath, str(e)))
                    else:
                        self._msg("""Extracted WifiAcquisitions successfully saved to "{}" !""".format(outPath))
            else:
                self._error("""Extraction failed""")

        else:
            self._error("""No acquisition loaded to extract from.""")
    
    _COMMANDS = {
        "convert":  {"func":convert, "help":"convert from a WifiLog format into another one"},
        "test":     {"func":test,    "help":"perform a test"},
        "extract":  {"func":extract, "help":"extract acquisition(s) from WifiLog(s) and save result to another WifiLog"}
    }
    
    def main(self):
        self._COMMANDS[self._options.command]["func"](self)
    
    def _initArgParser(self):
        _cmdLineArgParser = argparse.ArgumentParser(prog="WifiLog")
        _cmdLineArgParser.add_argument('-V', '--version', action="version", version="""%(prog)s {}""".format(self.VERSION))
        
        # Positional exclusive arguments
        sp = _cmdLineArgParser.add_subparsers(dest='command')
 
        for cmd in self._COMMANDS.keys():
            self._COMMANDS[cmd]["parser"] = sp.add_parser(cmd, help=self._COMMANDS[cmd]["help"])
            self._COMMANDS[cmd]["parser"].add_argument('path', help="path to WifiLog file(s) (only accept rawlog or log format)")
            self._COMMANDS[cmd]["parser"].add_argument("-v", "--verbose", action="store_true", help='debug mode (verbose)')
            
        # Convert command arguments
        self._COMMANDS["convert"]["parser"].add_argument("type", choices=("raw", "log", "xml", "csv"), help="type to convert file(s) into")
        self._COMMANDS["convert"]["parser"].add_argument("--outPath", "-o", type=str, help="directory where to store converted file(s)")
        self._COMMANDS["convert"]["parser"].add_argument("--override", action="store_true", help='override file(s) if already exist (default: false)')
        self._COMMANDS["convert"]["parser"].add_argument("--macKeyLen", type=int, default=2, help='number of characters from the end of the mac address used to identify signal (2 char = 256 dimensions) (default: 2)')

        # Extract command arguments
        self._COMMANDS["extract"]["parser"].add_argument("--outFile", "-o", type=str, help="filename where to store extractions (default: same as path suffixed by 'out')")
        self._COMMANDS["extract"]["parser"].add_argument("--type", "-t", choices=("raw", "log", "auto"), default="auto", help="type of input WifiLog file(s) to read (default: 'auto')")
        self._COMMANDS["extract"]["parser"].add_argument("--keepDuplicate", action="store_true", help='do not remove duplicated acquisitions (default: false)')
        self._COMMANDS["extract"]["parser"].add_argument("--keepEmpty", action="store_true", help='do not remove empty acquisitions (default: false)')
         
        fileResolution = self._COMMANDS["extract"]["parser"].add_mutually_exclusive_group()
        fileResolution.add_argument("--override", action="store_true", help='override output file if already exist (default: false)')
        fileResolution.add_argument("--append", action="store_true", help='append output file if already exist (default: false)')
        
        filters = self._COMMANDS["extract"]["parser"].add_argument_group(description="Extraction filters : only acquisitions that satisfy those filters will be extracted.")
        filters.add_argument("--fromDate", type=str, help="Keep acquisitions with a more recent date than FROMDATE. (format : '{}' or '{}')".format(WifiLogFile.Format.Date, WifiLogFile.Format.DateTime).replace("%", "%%"))
        filters.add_argument("--toDate", type=str, help="Keep acquisitions with a more ancient date than TODATE. (format : '{}' or '{}')".format(WifiLogFile.Format.Date, WifiLogFile.Format.DateTime).replace("%", "%%"))
         
        macList = filters.add_mutually_exclusive_group()
        macList.add_argument("--allowedMac", type=str, help="Keep acquisitions which only contains WifiSources from an allowed MAC address listed in ALLOWEDMAC file. ALLOWEDMAC file must only contains a MAC address per line.")
        macList.add_argument("--forbiddenMac", type=str, help="Keep acquisitions which doesn't contain any WifiSources from a forbidden MAC address listed in FORBIDDENMAC file. FORBIDDENMAC file must only contains a MAC address per line.")
         
        # Label
        filters.add_argument("--building", type=str, help="Keep acquisitions containing label BUILDING.")
        filters.add_argument("--colleagues", "-c", action="append", type=str, help="Keep acquisitions containing label COLLEAGUES")
        filters.add_argument("--comment", type=str, help="Keep acquisitions containing label COMMENT")
        filters.add_argument("--room", type=str, help="Keep acquisitions containing label ROOM")
        filters.add_argument("--floor", type=str, help="Keep acquisitions containing label FLOOR")
        
        # Author
        filters.add_argument("--authorName", type=str, help="Keep acquisitions containing label AUTHORNAME")
        filters.add_argument("--authorMac", type=str, help="Keep acquisitions containing label AUTHORMAC")

        filters.add_argument("--minSpc", type=float, help="Keep acquisitions with a space between two captures greater than MINSPC.")
        filters.add_argument("--maxSpc", type=float, help="Keep acquisitions with a space between two captures lower than MAXSPC.")
        filters.add_argument("--minDur", type=int, help="Keep acquisitions with a capture duration greater than MINDUR.")
        filters.add_argument("--maxDur", type=int, help="Keep acquisitions with a capture duration lower than MAXDUR.")
        filters.add_argument("--minCap", type=int, help="Keep acquisitions with at least MINCAP captures.")
        filters.add_argument("--maxCap", type=int, help="Keep acquisitions with at most MAXCAP captures.")
        filters.add_argument("--minSrc", type=int, help="Keep acquisitions with at least MINSRC sources.")
        filters.add_argument("--maxSrc", type=int, help="Keep acquisitions with at most MAXSRC sources.")
        filters.add_argument("--minPck", type=int, help="Keep acquisitions with at least MINPCK packets.")
        filters.add_argument("--maxPck", type=int, help="Keep acquisitions with at most MAXPCK packets.")
        filters.add_argument("--minMed", type=float, help="Keep acquisitions with a median distance greater than MINMED.")
        filters.add_argument("--maxMed", type=float, help="Keep acquisitions with a median distance lower than MAXMED.")
        filters.add_argument("--minMin", type=float, help="Keep acquisitions with a minimum distance greater than MINMIN.")
        filters.add_argument("--maxMin", type=float, help="Keep acquisitions with a minimum distance greater than MAXMIN.")
        filters.add_argument("--minMax", type=float, help="Keep acquisitions with a maximum distance greater than MINMAX.")
        filters.add_argument("--maxMax", type=float, help="Keep acquisitions with a maximum distance greater than MAXMAX.")
        filters.add_argument("--minAvg", type=float, help="Keep acquisitions with an average distance greater than MINAVG.")
        filters.add_argument("--maxAvg", type=float, help="Keep acquisitions with an average distance greater than MAXAVG.")
        filters.add_argument("--minStd", type=float, help="Keep acquisitions with a standard distance greater than MINSTD.")
        filters.add_argument("--maxStd", type=float, help="Keep acquisitions with a standard distance greater than MAXSTD.")

        return _cmdLineArgParser
    
    def _getBaseFileName(self, fileName):
        baseFileName = '.'.join(fileName.split('.')[:-1])
        return baseFileName if baseFileName else fileName
    
    def _msg(self, message):
        print(message)
    
    def _error(self, message):
        print("ERROR: " + message)
    
    def _debug(self, message):
        if self._options.verbose:
            print("DEBUG: " + message)
    

if __name__ =='__main__':
    
    wifiLog = WifiLog()
    wifiLog.main()
    
    #cont = WifiLogFileLoader("test/2015-01-06.rawlog").load()
    #print WifiLogFileDumper(csvOptions=WifiLogFileDumper.CsvOptions(",")).dump(cont, mode=WifiLogFile.Type.Csv)
    
    if False:
        with open("out.csv", "w") as csvfile:
            fieldnames = ['first_name', 'last_name','third']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writeheader()
            writer.writerow({'first_name' : "toto", "last_name":'tata'})
            writer.fieldnames
             
    if False:
        with open('out.csv') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                print row
                
