#-*- coding: utf-8 -*-

import subprocess
import platform
import math
import time 
import os
import re
import commands
from datetime import datetime as datetime
import argparse

import WifiObjects as wo
import WifiLog as wl

__GLOBAL_VERSION = "2.0"

class WiresharkBridge():
    """Bridge between Wireshark and WifiObjects"""
    
    # Usefull for parsing TShark files
    _SOURCE      = "mac_source"
    _DESTINATION = "mac_destination"
    _POWER       = "power"
    _FREQUENCY   = "frequency"
    
    # Lines TShark can return :
    # - Type 1 (accept): 9c:4e:36:cb:6c:24 -> 01:00:5e:7f:ff:fa -48 5240
    # - Type 2 (refuse):             -> 60:f8:1d:a8:6e:3a (RA) -48 5240
    # - Type 3 (accept): 9c:1c:12:2b:c4:79 (TA) -> 60:f8:1d:a8:6e:3a (RA) -48 5240
    _TSHARK_REGEX = re.compile(r"(?P<" + _SOURCE + ">(?:[a-f0-9]{2}:){5}[a-f0-9]{2})(?: \(TA\))? -> (?P<" + _DESTINATION + ">(?:[a-f0-9]{2}:){5}[a-f0-9]{2})(?: \(RA\))? (?P<" + _POWER + ">-?[0-9]+) (?P<" + _FREQUENCY + ">[0-9]+)", re.IGNORECASE)

    _TSHARK_CMD = """{tshark_path} -i {capturing_interface} -I -n -o gui.column.format:'"Source", "%s", "Destination", "%d", "Pwr", "%Cus:radiotap.dbm_antsignal", "Freq", "%Cus:radiotap.channel.freq"' -a duration:{capture_duration} -T text > {tshark_log_file}"""

    def __init__(self, tsharkPath, tmpLogFile, excludeMac, interface, capDuration, forbiddenMac=True, debug=False):
        """Bridge to use TShark in python.
        'tsharkPath' is the path to TShark binary.
        'tmpLogFile' is the name of the temp file which will store TShark output. 
        'excludeMac' is a list of forbidden mac addresses.
        'interface' is the network interface to use to capture.
        'capDuration' is the duration in seconds of a capture.
        Set 'forbiddenMac' to false to use 'excludeMac' as a list of exclusively allowed mac.
        Set 'debug' to true to use verbose mode."""
        self._commandLine = self._TSHARK_CMD.format(tshark_path=tsharkPath, capturing_interface=interface, capture_duration=capDuration, tshark_log_file=tmpLogFile)
        self._tmpLogFile = str(tmpLogFile)
        self._capDuration = int(capDuration)
        self._excludedMacList = list()
        self._forbiddenMac = forbiddenMac
        self._verbose = debug
                
        if excludeMac:
            with open(excludeMac, "r") as excludedMacFile:
                self._excludedMacList = excludedMacFile.readlines()
    
    def _capture(self):
        start = datetime.now()
        self._msg("""New capture launched at "{}" for {}s :""".format(start, self._capDuration))
        
        # Call TShark
        self._debug("""Calling following command line : "{}" """.format(self._commandLine))
        tsharkResult = subprocess.call(self._commandLine, shell=True)
        
        if tsharkResult == 0:
            wifiCapture = wo.WifiRawCapture(start)
            
            # Parse TShark file
            (_, sources) = self._parseTSharkFile()
            
            for source in sources.keys():
                wifiCapture.addSource(wo.WifiRawSource(source, sources[source]))
            
            self._msg("Capture succeed !")
            return wifiCapture
        
        self._error("Capture failed !")
        return None
        
    def _promptUserQuestion(self, question):
        answer = raw_input("{} (default: Yes) : ".format(question))
        
        if answer.strip().lower() in ["no", "n"]:
            return False
        
        return True

    def acquire(self, capNumber, capSpacing, wifiAuthor, wifiLabel):
        """Launch a Wifi acquisition.
        'capNumber' is the number of captures to perform. If 'capNumber' < 1, acquisition is manual (after each capture, user will be prompted to continue).
        'capSpacing' is time in seconds between two captures.
        'wifiAuthor' is the author of this acquisition (meta-data).
        'wifiLabel' is the label of this acquisition (meta-data)."""
        assert(isinstance(capNumber, int))
        capturesCount = 0
        continueAcquisition = True  
        
        # In manual acquisition, capSpacing makes no sense
        if capNumber < 1:
            capSpacing = -1.
            
        # A negative capSpacing also means a manual acquisition
        if capSpacing < 0.:
            capNumber = 0
            
        wifiAcquisition = wo.WifiRawAcquisition(self._capDuration, capSpacing, wifiAuthor, wifiLabel)
        
        self._msg("""New {manual} acquisition launched by "{authorName}"...""".format(manual="manual" if capNumber < 1 else "automatic", authorName=wifiAuthor.name))
        
        while (capNumber > 0 and capturesCount < capNumber) or (capNumber < 1 and continueAcquisition):
            wifiCapture = self._capture()
            
            if not wifiCapture:
                break
            
            capturesCount += 1
            wifiAcquisition.addCapture(wifiCapture)

            if capSpacing >= 0. and capNumber > 0 and capturesCount < capNumber:
                self._msg("Capture Spacing : waiting {}s before next capture...".format(capSpacing))
                time.sleep(capSpacing)
                
            if capNumber < 1:
                continueAcquisition = self._promptUserQuestion("Would you like to launch another capture ?")
        
        self._msg("""Acquisition launched by "{authorName}" ended.""".format(authorName=wifiAuthor.name))
        
        return wifiAcquisition
    
    def _hashMacAddress(self, macAddress):
        return macAddress #not yet implemented
    
    def _parseTSharkFile(self):
        # Open TShark file and copy its content in a list (line by line)
        with open(self._tmpLogFile, "r") as tSharkFile:
            lines = tSharkFile.readlines()

        dictSourceDistance = {}
        
        for line in lines:
            result = self._TSHARK_REGEX.search(line)
            
            if result:
                self._debug(("Line matched : source={" + self._SOURCE + "}, dest={" + self._DESTINATION + "}, pwr={" + self._POWER + "}, freq={" + self._FREQUENCY + "}").format(**result.groupdict()))
                                
                if (self._forbiddenMac and result.group(self._SOURCE) not in self._excludedMacList) or (not self._forbiddenMac and result.group(self._SOURCE) in self._excludedMacList):
                    source = self._hashMacAddress(result.group(self._SOURCE))
                    dictSourceDistance.setdefault(source, []).append(float(WifiDistanceEstimator(result.group(self._POWER), result.group(self._FREQUENCY))))
            else:
                self._debug("""Following line didn't match: "{}" """.format(line.rstrip()))
            
        return (len(lines), dictSourceDistance)

    def _msg(self, message):
        print(message)
    
    def _error(self, message):
        print("ERROR: " + message)
    
    def _debug(self, message):
        if self._verbose:
            print("DEBUG: " + message)

class WifiDistanceEstimator:
    """Estimate a transmitter distance using Friis's formula. Recurrent values are cached."""

    # Cache frequencies coefficients
    frequencyCoefs = {}

    def __init__(self, transmitterPower, frequency):
        """transmitterPower (dB), frequency (MHz)"""
        
        self.transmitterPower = abs(float(transmitterPower))
        self.frequency = float(frequency)
        
    def compute(self):
        """Compute transmitter distance"""
    
        if self.frequency not in WifiDistanceEstimator.frequencyCoefs.keys():
            WifiDistanceEstimator.frequencyCoefs[self.frequency] = 27.55 - (20 * math.log10(self.frequency))
            
        return 10**((self.transmitterPower + WifiDistanceEstimator.frequencyCoefs[self.frequency]) / 20)
            
    def __float__(self):
        return self.compute()

@wo.Singleton
class WifiCapture:
    """WifiCapture command line interface"""
    
    VERSION = "2.0"
    
    def __init__(self):
        self._cmdLineArgParser = self._initArgParser()
        self._options = self._cmdLineArgParser.parse_args()
        self._authorMac = None
                
    def _resume(self):
        try:
            container = wl.WifiLogFileLoader(self._options.filename).load()
        except IOError, e:
            self._error("""Can't load given filename ("{}") in order to resume capture. ("{}")""".format(self._options.filename, e))
            return
        
        if container and container.hasAcquisitions():
            lastAcquisition = container.getLastAcquisition(True)
            self._launch(lastAcquisition.duration, lastAcquisition.spacing, lastAcquisition.getCapturesCount(), lastAcquisition.label, self._options.filename, wl.WifiLogFileDumper.ConflictResolution.Append, os.path.dirname(self._options.filename))
        else:
            self._error("""Can't resume because given file doesn't contain any capture. ("{}")""".format(self._options.filename))
    
    def _new(self):
        filename = self._getOutputFileName()
        fileResolution = wl.WifiLogFileDumper.ConflictResolution.Override if self._options.override else wl.WifiLogFileDumper.ConflictResolution.Append
        label = wo.WifiAcquisitionLabel(self._options.building, self._options.floor, self._options.room, self._options.colleague, self._options.comment)
        
        self._launch(self._options.capDuration, self._options.capSpacing, self._options.capNumber, label, filename, fileResolution, self._options.pathToLog)
        
    def _launch(self, duration, spacing, capNumber, label, filename, fileResolution, pathToLog):
        (macListFileName, forbiddenMac) = self._getMacListRestriction()
        author = wo.WifiAcquisitionAuthor(self._options.author, self._authorMac)
        
        bridge = WiresharkBridge(
            tsharkPath=self._options.tsharkBin,
            tmpLogFile=self._options.tsharkOut,
            excludeMac=macListFileName,
            forbiddenMac=forbiddenMac,
            interface=self._options.interface,
            debug=self._options.verbose,
            capDuration=duration
        )
        
        wifiAcquisition = bridge.acquire(capNumber, spacing, author, label)
        
        # Writing raw file
        outPath = os.path.join(pathToLog, filename);
        wl.WifiLogFileDumper(outPath).dump(wifiAcquisition, fileResolution=fileResolution, mode=wl.WifiLogFile.Type.Raw)
        self._msg("Acquisition saved ({}) !".format(outPath))
    
    def main(self):
        if not self._getSystemDefaultsOptions():
            return False
        
        if self._options.command == "resume":
            self._resume()
        elif self._options.command == "new":
            self._new()
            
        return True
    
    def _addCommonArguments(self, parser):
        parser.add_argument('-v', '--verbose', action="store_true", help='debug mode')
        parser.add_argument('--tsharkOut', type=str, default="tshark.log", help="Temporary TShark output file (default: 'tshark.log')")
        parser.add_argument('--tsharkBin', type=str, default=None, help='Path to TShark binary (leave blank to use system default)')
        parser.add_argument('--interface', '-i', type=str, default=None, help='Wireless interface (leave blank to use system default)')
        parser.add_argument('author', type=str, help='Name of the author of this capture')
        
        macList = parser.add_mutually_exclusive_group()
        macList.add_argument("--allowedMac", type=str, help="Only captures WifiSources with an allowed MAC address listed in ALLOWEDMAC file. ALLOWEDMAC file must only contains a MAC address per line.")
        macList.add_argument("--forbiddenMac", type=str, help="Ignore WifiSource with a forbidden MAC address listed in FORBIDDENMAC file. FORBIDDENMAC file must only contains a MAC address per line.")
        
        
    def _initArgParser(self):
        _cmdLineArgParser = argparse.ArgumentParser(prog="WifiCapture")
        _cmdLineArgParser.add_argument('-V', '--version', action="version", version="""%(prog)s {}""".format(self.VERSION))
        
        subParsers = _cmdLineArgParser.add_subparsers(dest='command')
        
        ### Resume command ###
        subParserResume = subParsers.add_parser("resume", help="Resume a previous acquisition")
        self._addCommonArguments(subParserResume)
        
        # Resume command positional arguments
        subParserResume.add_argument('filename', type=str, help='File to resume')
        
        
        ### New command ###
        subParserNew = subParsers.add_parser("new", help="Perform a new acquisition")
        self._addCommonArguments(subParserNew)
        
        # New command label arguments
        label = subParserNew.add_argument_group(description="Label : describe current context and identify acquisition.")
        label.add_argument('building', type=str, help='Current building location.')
        label.add_argument('floor', type=int, help='Current floor location.')
        label.add_argument('room', type=str, help='Current room location.')
        label.add_argument('-c', '--colleague', action="append", type=str, default=None, help='List of current co-workers. Use multiple times to specify many colleagues.')
        label.add_argument('-cmt', '--comment', type=str, default=None, help='Comment.')
        
        # New command arguments
        subParserNew.add_argument('--capDuration', '-d',    type=int,   default=3, help='Duration in seconds of each capture (default: 3)')
        subParserNew.add_argument('--capNumber', '-n',      type=int,   default=0, help='Number of capture(s). Set CAPNUMBER to 0 to switch in manual mode and manually trigger each acquisition. If CAPNUMBER > 0, acquisition is in automatic mode and CAPSPACING is used. (default: 0)')
        subParserNew.add_argument('--capSpacing', '-s',     type=float, default=0., help='Time to wait in seconds between two captures in automatic mode (see also CAPNUMBER, default: 0.0)')
        subParserNew.add_argument('--pathToLog', '-p',      type=str,   default="./", help='Path to the folder (where the log will be saved), default : current folder')
        subParserNew.add_argument('--filename', '-o',       type=str,   default=None, help='How to name the file')
        subParserNew.add_argument('--override', '-w',       action="store_true", help='If output file already exists, override it instead of appending it with new acquisition. (default: false)')
        
        return _cmdLineArgParser
    
    _OS_SPECIFICITY = {
        "Windows":{"platform":("windows")},
        "Linux":{
            "platform":("linux", "bsd", "kali", "solaris"),
            "tsharkBin":"tshark",
            "interface":"wlan0",
            "ifconfigMAC":"HWaddr"
        },
        "Mac":{
            "platform":("mac", "darwin", "osx"),
            "tsharkBin":"/Applications/Wireshark.app/Contents/Resources/bin/tshark",
            "interface":"en0",
            "ifconfigMAC":"ether"
        }
    }
    
    def _getCurrentOS(self):
        currentOS = platform.system().lower()
        
        for os in self._OS_SPECIFICITY.keys():
            if any(currentOS in s for s in self._OS_SPECIFICITY[os]["platform"]):
                return os
            
        return None
    
    def _getSystemDefaultsOptions(self):
        if self._options.tsharkBin == None or self._options.interface == None:
            self._debug("Missing argument(s)")
            self._debug("""Trying to guess parameters from OS version, use "python WifiCapture.py -h" for help""")
            
            currentOS = self._getCurrentOS()
            
            self._debug("""Detected OS : "{}" """.format(currentOS))
            
            if not currentOS:
                self._debug("os not in database, assuming it's linux based")
                currentOS = "Linux"
            
            if currentOS == "Windows":
                self._error("sorry, monitoring is not supported on windows, install linux")
                return False
                
            if self._options.tsharkBin == None:
                self._options.tsharkBin = self._OS_SPECIFICITY[currentOS]["tsharkBin"]
                self._debug("""no TShark path specified, using most common path for os "{}" ({})""".format(currentOS, self._options.tsharkBin))

            if self._options.interface == None:
                self._options.interface = self._OS_SPECIFICITY[currentOS]["interface"]
                self._debug("""no interface specified, using most common interface for os "{}" ({})""".format(currentOS, self._options.interface))
                
            self._authorMac = self._getInterfaceMac(currentOS, self._options.interface)
        
        return True
    
    def _getInterfaceMac(self, currentOS, interface):
        ifconfigOutput = commands.getoutput("ifconfig " + interface).split()
        
        if self._OS_SPECIFICITY[currentOS]["ifconfigMAC"] in ifconfigOutput:
            return ifconfigOutput[ ifconfigOutput.index(self._OS_SPECIFICITY[currentOS]["ifconfigMAC"]) + 1 ]

        return None
    
    def _getOutputFileName(self):
        formatDict = {
            "date":time.strftime("%Y-%m-%d", time.localtime()),
            "name":"_" + self._options.filename if self._options.filename else ""
        }
        return "{date}{name}.rawlog".format(**formatDict)
    
    def _getMacListRestriction(self):
        macListFileName = self._options.forbiddenMac
        forbiddenMac = True
        
        if self._options.allowedMac:
            macListFileName = self._options.allowedMac
            forbiddenMac = False
            
        return (macListFileName, forbiddenMac)

    def _msg(self, message):
        print(message)
    
    def _error(self, message):
        print("ERROR: " + message)
    
    def _debug(self, message):
        if self._options.verbose:
            print("DEBUG: " + message)

    
if __name__ =='__main__':
    wifiCapture = WifiCapture()
    wifiCapture.main()
    