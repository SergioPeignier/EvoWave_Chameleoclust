#-*- coding: utf-8 -*-
import numpy as np
import copy
import operator
from collections import Iterable
from datetime import datetime as datetime
from types import NoneType

_VERSION = "2.0"
_DEBUG = False
THRESHOLD_NUMBER_OF_VALUES_PER_MAC_ADDRESS = 120
def Singleton(singletonClass):
    """Singleton pattern as decorator"""
    
    instances = {}
    
    def getInstance():
        if singletonClass not in instances:
            instances[singletonClass] = singletonClass()
        return instances[singletonClass]
    return getInstance

def Debug(message):
    """Print debug message if _DEBUG is true"""
    
    if _DEBUG:
        if isinstance(message, Iterable) and not isinstance(message, basestring):
            for item in message:
                Debug(item)
        else:
            print "DBG: {}".format(message)
        
class WifiObject(object):
    """Base class of every WifiObject"""
    
    def isWifiRawObject(self):
        return False
    
    @staticmethod
    def _RemoveAllExcept(itemsList, exceptions):
        Debug([item for item in itemsList if item not in exceptions])
        itemsList[:] = [item for item in itemsList if item in exceptions]

    @staticmethod
    def _Remove(itemsList, items):
        itemsList[:] = [item for item in itemsList if item not in items]
    
    @staticmethod
    def _RemoveDuplicates(itemsList):
        itemsList = list(set(itemsList))
    
    
class WifiRawObject(WifiObject):
    """Base class of every WifiRawObject. A raw version of a WifiObject is a version containing the whole data of a WifiCapture, allowing to recompute each of its state."""
    def isWifiRawObject(self):
        return True
        
        
class WifiObjectRawTypeError(TypeError):
    """Type Error which occurs on trying to mix WifiObjects and WifiRawObjects in the same container"""
    
    def __init__(self, requiredType, givenType):
        self.requiredType = requiredType
        self.givenType = givenType
        
    def __str__(self):
        return """Type "{}" expected: "{}" found.""".format(self.requiredType, self.givenType)

        
class WifiObjectsSelector(object):
    """Usefull class used in extraction. Basically, a selection is a subset of items satisfying some criteria."""
    
    def __init__(self, objectsContainer, evalMethod = None, attributeName=None):
        """'objectsContainer' contains items which will be evaluated by 'evalMethod' or using its 'attributeName'"""
        self.objectsContainer = objectsContainer
        self.evalMethod = evalMethod
        self.attributeName = attributeName
        
    def _getValue(self, obj):
        if self.evalMethod:
            return self.evalMethod(obj)
        
        if self.attributeName:
            return getattr(obj, self.attributeName)
        
        return obj
    
    def select(self, value, equalMethod=operator.__eq__):
        """Returns items in objectsContainer equals to 'value' using 'equalMethod'"""
        if not isinstance(self.objectsContainer, Iterable):
            return None
        
        if value:
            return filter(lambda obj : equalMethod(self._getValue(obj), value), self.objectsContainer)
            
        return self.objectsContainer
    
    def selectBetween(self, lowerBound=None, upperBound=None):
        """Returns items in objectsContainer with a greater value than 'lowerBound' (if specified) and a lower value than 'upperBound' (if specified), using common operators '>' and '<'"""
        if not isinstance(self.objectsContainer, Iterable):
            return None
        
        result = self.objectsContainer
        
        if lowerBound:
            result = filter(lambda obj : self._getValue(obj) > lowerBound, result)
        
        if upperBound:
            result = filter(lambda obj : self._getValue(obj) < upperBound, result)
            
        return result 
    
    def _isObjInList(self, obj, inList):
        if isinstance(obj, Iterable) and not isinstance(obj, basestring):
            for item in obj:
                if self._isObjInList(item, inList):
                    return True
            return False
        
        return obj in inList
    
    def selectIn(self, selectList, forbiddenList=False):
        """If 'forbiddenList' is false, return items from objectsContainer for which every value is in 'selectList'.
        If 'forbiddenList' is true, return items from objectsContainer for which any value is not in 'selectList'."""
        if not isinstance(self.objectsContainer, Iterable):
            return None

        if forbiddenList:
            filterObj = lambda obj : not self._isObjInList(self._getValue(obj), selectList)
        else:
            filterObj = lambda obj : self._isObjInList(self._getValue(obj), selectList)
        
        return filter(filterObj, self.objectsContainer)
    
    
class WifiObjectsAggregator(WifiObjectsSelector):
    """Usefull class used to aggregate items in containers. 
    Aggregate means a way to group some objects contained in a tree of containers (i.e. a container of container of ... of objects to aggregate).
    Derives for WifiObjectsSelector to use mechanism of getting a value using its attribute name or an evaluating method."""
    
    def getSum(self):
        """Returns the sum of aggregated objects (using common operator '+')"""
        if not self.objectsContainer:
            return None
        
        result = self._getValue(self.objectsContainer[0])
        
        for obj in self.objectsContainer[1:]:
            objValue = self._getValue(obj)
            if objValue:
                result += objValue
            
        return result
    
    def aggregateList(self, removeDuplicates=False):
        """Returns every items in a '1-depth' list. 'objectsContainer' must be a list."""
        assert isinstance(self.objectsContainer, list)
        result = list()
        
        for obj in self.objectsContainer:
            objValue = self._getValue(obj)
            if objValue:
                if isinstance(objValue, list):
                    result.extend(objValue)
                else:
                    result.append(objValue)
            
        if removeDuplicates:
            result = list(set(result))
                    
        return result
    

class WifiStats(object):
    """Represents statistics of a set."""
    
    # Precision
    PRECISION_DIST  = 2
    PRECISION_STATS = 3
    
    def __init__(self, minDistance=None, medDistance=None, maxDistance=None, distancesAvg=None, distancesStd=None, getDataMethod=None):
        """Stats can be computed or just stored. If 'getDataMethod' is specified, stats could be computed on data returned by this function."""
        self.minDistance = minDistance
        self.medDistance = medDistance
        self.maxDistance = maxDistance
        self.distancesAvg = distancesAvg
        self.distancesStd = distancesStd
        self.getDataMethod = getDataMethod
        
    def getStatistics(self, precisionStats = PRECISION_STATS, precisionDist = PRECISION_DIST):
        """Compute stats (if 'getDataMethod' specified) and return a dictionary containing statistics"""
        if self.getDataMethod:
            data = self.getDataMethod()
            if len(data) > THRESHOLD_NUMBER_OF_VALUES_PER_MAC_ADDRESS:
                self.update(min(data), np.median(data), max(data), np.average(data), np.std(data))
            else :
                self.update(np.nan, np.nan, np.nan, np.nan, np.nan)
                    
        result = {key : (round(self.__dict__[key], precisionDist) if isinstance(self.__dict__[key], float) else self.__dict__[key]) for key in ["maxDistance", "minDistance", "medDistance"]}
        result = dict({key : (round(self.__dict__[key], precisionStats) if isinstance(self.__dict__[key], float) else self.__dict__[key]) for key in ["distancesAvg", "distancesStd"]}.items() + result.items())
        return result
    
    def update(self, minDistance = float("+inf"), medDistance = None, maxDistance = float("-inf"), distancesAvg = None, distancesStd = None):
        """Manually updates statistics""" 
        self.minDistance = minDistance
        self.medDistance = medDistance
        self.maxDistance = maxDistance
        self.distancesAvg = distancesAvg
        self.distancesStd = distancesStd
        
    def __repr__(self):
        values = ""
        stats = self.getStatistics()
        
        for stat in stats.keys():
            if stats[stat]:
                values += """{}="{}" """.format(stat, stats[stat])
                
        return "<WifiStats {}>".format(values.strip())
    
    def __le__(self, other):
        return (
            (not self.minDistance     or not other.minDistance     or self.minDistance    <= other.minDistance)   and 
            (not self.maxDistance     or not other.maxDistance     or self.maxDistance    <= other.maxDistance)   and 
            (not self.medDistance     or not other.medDistance     or self.medDistance    <= other.medDistance)   and 
            (not self.distancesAvg    or not other.distancesAvg    or self.distancesAvg   <= other.distancesAvg)  and 
            (not self.distancesStd    or not other.distancesStd    or self.distancesStd   <= other.distancesStd)
        )
    
    def __lt__(self, other):
        return (
            (not self.minDistance     or not other.minDistance     or self.minDistance    < other.minDistance)   and 
            (not self.maxDistance     or not other.maxDistance     or self.maxDistance    < other.maxDistance)   and 
            (not self.medDistance     or not other.medDistance     or self.medDistance    < other.medDistance)   and 
            (not self.distancesAvg    or not other.distancesAvg    or self.distancesAvg   < other.distancesAvg)  and 
            (not self.distancesStd    or not other.distancesStd    or self.distancesStd   < other.distancesStd)
        )
        
    def __ge__(self, other):
        return (
            (not self.minDistance     or not other.minDistance     or self.minDistance    >= other.minDistance)   and 
            (not self.maxDistance     or not other.maxDistance     or self.maxDistance    >= other.maxDistance)   and 
            (not self.medDistance     or not other.medDistance     or self.medDistance    >= other.medDistance)   and 
            (not self.distancesAvg    or not other.distancesAvg    or self.distancesAvg   >= other.distancesAvg)  and 
            (not self.distancesStd    or not other.distancesStd    or self.distancesStd   >= other.distancesStd)
        )
        
    def __gt__(self, other):
        return (
            (not self.minDistance     or not other.minDistance     or self.minDistance    > other.minDistance)   and 
            (not self.maxDistance     or not other.maxDistance     or self.maxDistance    > other.maxDistance)   and 
            (not self.medDistance     or not other.medDistance     or self.medDistance    > other.medDistance)   and 
            (not self.distancesAvg    or not other.distancesAvg    or self.distancesAvg   > other.distancesAvg)  and 
            (not self.distancesStd    or not other.distancesStd    or self.distancesStd   > other.distancesStd)
        )

class WifiAcquisitionsContainer(WifiObject):
    """Container of WifiAcquisitions"""
    
    def __init__(self):
        self.acquisitions = list()
        
        # Aggregators
        self._capturesAggregator = WifiObjectsAggregator(self.acquisitions, WifiAcquisition.getCapturesCount)
        self._sourcesAggregator = WifiObjectsAggregator(self.acquisitions, WifiAcquisition.getSourcesCount)
        self._packetsAggregator = WifiObjectsAggregator(self.acquisitions, WifiAcquisition.getPacketsCount)
        self._macAggregator = WifiObjectsAggregator(self.acquisitions, WifiAcquisition.getMacList)
        
        # Selectors
        self._statsSelector = WifiObjectsSelector(self.acquisitions, attributeName="stats")
        self._spacingSelector = WifiObjectsSelector(self.acquisitions, attributeName="spacing")
        self._durationSelector = WifiObjectsSelector(self.acquisitions, attributeName="duration")
        self._authorSelector = WifiObjectsSelector(self.acquisitions, attributeName="author")
        self._labelSelector = WifiObjectsSelector(self.acquisitions, attributeName="label")
    
    def addAcquisition(self, wifiAcquisition):
        if type(wifiAcquisition) is not WifiAcquisition:
            raise WifiObjectRawTypeError("WifiAcquisition", type(wifiAcquisition))

        self.acquisitions.append(wifiAcquisition)
        
    def append(self, wifiAcquisitionContainer):
        """Merge two WifiAcquisitions : add each WifiAcquisition of 'wifiAcquisitionContainer' into self"""
        for acquisition in wifiAcquisitionContainer.acquisitions:
            self.addAcquisition(acquisition)
        
    def hasAcquisitions(self):
        """Returns true if this container contains at least one WifiAcquisition (empty or not)"""
        return len(self.acquisitions) > 0
    
    def isEmpty(self):
        """Returns true if this container doesn't contain at least one non-empty WifiAcquisition"""
        for acquisition in self.acquisitions:
            if not acquisition.isEmpty():
                return False
        return True
    
    def getLastAcquisition(self, datetimeOrdered = False):
        """Returns the last acquisition.
        If 'datetimeOrdered' is true, the last means the last WifiAcquisition ordered by its start time.
        If 'datetimeOrdered' is false, the last means the last WifiAcquisition ordered by its index in this container (the last added)."""
        if len(self.acquisitions) < 1:
            return None
        
        if datetimeOrdered:
            return max(self.acquisitions, key=lambda acquisition : acquisition.getFirstCapture(True))
        
        return self.acquisitions[-1]
        
    def getFirstAcquisition(self, datetimeOrdered = False):
        """Returns the first acquisition.
        If 'datetimeOrdered' is true, the first means the first WifiAcquisition ordered by its start time.
        If 'datetimeOrdered' is false, the last first the first WifiAcquisition ordered by its index in this container (the first added)."""
        if len(self.acquisitions) < 1:
            return None
        if datetimeOrdered:
            return min(self.acquisitions, key=lambda acquisition : acquisition.getFirstCapture(True))
        else:
            return self.acquisitions[-1]
        
    def selectByDatetime(self, fromDatetime, toDatetime = datetime.now()):
        """Returns acquisitions which contains captures done between 'fromDatetime' and 'toDatetime'"""
        acquisitions = list()
        
        for acquisition in self.acquisitions:
            if not acquisition.hasCapture():
                continue
            
            selected = True
            
            if fromDatetime:
                firstStart = acquisition.getFirstCapture(True).start
                selected = firstStart and firstStart > fromDatetime
                
            if toDatetime:
                lastStart = acquisition.getLastCapture(True).start
                selected = selected and (lastStart and lastStart < toDatetime)
            
            if selected:
                acquisitions.append(acquisition)
                
        return acquisitions
    
    def extract(self,
                minDuration=None, maxDuration=None,
                minSpacing=None, maxSpacing=None,
                minCaptures=None, maxCaptures=None,
                author=None, label=None,
                fromDatetime=None, toDatetime=None,
                minSources=None, maxSources=None,
                minPackets=None, maxPackets=None,
                minStats=None, maxStats=None,
                macList=list(), forbiddenList=True):
        """Return a copy of the current WifiAcquisitionContainer which only contains WifiAcquisitions that satisfy all criteria"""
        container = copy.deepcopy(self)# Needed to keep internal types (Raw vs not Raw)

        Debug("""Select acquisitions with date between "{}" and "{}":""".format(fromDatetime, toDatetime))
        self._RemoveAllExcept(container.acquisitions, container.selectByDatetime(fromDatetime, toDatetime))
        
        Debug("""Select acquisitions with spacing between "{}" and "{}":""".format(minSpacing, maxSpacing))
        self._RemoveAllExcept(container.acquisitions, container._spacingSelector.selectBetween(minSpacing, maxSpacing))
        
        Debug("""Select acquisitions with duration between "{}" and "{}":""".format(minDuration, maxDuration))
        self._RemoveAllExcept(container.acquisitions, container._durationSelector.selectBetween(minDuration, maxDuration))
        
        Debug("""Select acquisitions with author "{}":""".format(author))
        self._RemoveAllExcept(container.acquisitions, container._authorSelector.select(author))
        
        Debug("""Select acquisitions with label "{}":""".format(label))
        self._RemoveAllExcept(container.acquisitions, container._labelSelector.select(label, WifiAcquisitionLabel.extractionComparision))
        
        Debug("""Select acquisitions with captures between "{}" and "{}":""".format(minCaptures, maxCaptures))
        self._RemoveAllExcept(container.acquisitions, container._capturesAggregator.selectBetween(minCaptures, maxCaptures))
        
        Debug("""Select acquisitions with MAC addresses {} in "{}":""".format("not" if forbiddenList else "", macList))
        self._RemoveAllExcept(container.acquisitions, container._macAggregator.selectIn(macList, forbiddenList))
        
        Debug("""Select acquisitions with sources between "{}" and "{}":""".format(minSources, maxSources))
        self._RemoveAllExcept(container.acquisitions, container._sourcesAggregator.selectBetween(minSources, maxSources))
        
        Debug("""Select acquisitions with packets between "{}" and "{}":""".format(minPackets, maxPackets))
        self._RemoveAllExcept(container.acquisitions, container._packetsAggregator.selectBetween(minPackets, maxPackets))
        
        Debug("""Select acquisitions with stats between "{}" and "{}":""".format(minStats, maxStats))
        self._RemoveAllExcept(container.acquisitions, container._statsSelector.selectBetween(minStats, maxStats))
        
        return container

    def removeDuplicates(self):
        """Remove duplicated WifiAcquisitions"""
        self._RemoveDuplicates(self.acquisitions)

    def getAcquisitionsCount(self):
        return len(self.acquisitions)
    
    def getCapturesCount(self):
        return self._capturesAggregator.getSum()
    
    def getPacketsCount(self):
        return self._packetsAggregator.getSum()
    
    def getSourcesCount(self):
        return self._sourcesAggregator.getSum()
    
    def getStatistics(self):
        stats = {"maxDistance":None, "minDistance":None, "medDistance":None, "distancesAvg":None, "distancesStd":None}
        
        for acquisition in self.acquisitions:
            acqStats = acquisition.stats.getStatistics()
            if acqStats["maxDistance"]:
                if not stats["maxDistance"]:
                    stats["maxDistance"] = float("-inf")
                stats["maxDistance"] = max(stats["maxDistance"], acqStats["maxDistance"])
                
            if acqStats["minDistance"]:
                if not stats["minDistance"]:
                    stats["minDistance"] = float("+inf")
                stats["minDistance"] = min(stats["minDistance"], acqStats["minDistance"])
                
        return stats
    
class WifiRawAcquisitionsContainer(WifiAcquisitionsContainer, WifiRawObject):
    """Container of WifiRawAcquisition"""
    
    def addAcquisition(self, wifiRawAcquisition):
        if type(wifiRawAcquisition) is not WifiRawAcquisition:
            raise WifiObjectRawTypeError("WifiRawAcquisition", type(wifiRawAcquisition))

        self.acquisitions.append(wifiRawAcquisition)

    def getStatistics(self):
        stats = WifiAcquisitionsContainer.getStatistics(self)
        
        allDistances = self.getAllDistances()
        if len(allDistances) > 0:
            stats["medDistance"] = np.median(allDistances)
            stats["distancesAvg"] = np.average(allDistances)
            stats["distancesStd"] = np.std(allDistances)
            
        return stats

    def getAllDistances(self):
        """Returns a list of all distances contained in this acquisition (including each source of each capture)"""
        return WifiObjectsAggregator(self.acquisitions, WifiRawAcquisition.getAllDistances).aggregateList()
    
class WifiAcquisition(WifiObject):
    """Represents a sequence of WifiCapture"""
    
    def __init__(self, duration, spacing, wifiAuthor, wifiLabel, wifiStats = None):
        if not isinstance(wifiAuthor, WifiAcquisitionAuthor):
            raise TypeError("WifiRawAcquisition requires type 'WifiAcquisitionAuthor' for parameter wifiAuthor. (Given type : '{}')".format(type(wifiAuthor)))
        
        if not isinstance(wifiLabel, WifiAcquisitionLabel):
            raise TypeError("WifiRawAcquisition requires type 'WifiAcquisitionLabel' for parameter wifiLabel. (Given type : '{}')".format(type(wifiLabel)))
        
        self.duration = duration
        self.spacing = spacing
        self.author = wifiAuthor
        self.label = wifiLabel
        self.stats = wifiStats if wifiStats else WifiStats()
        self._captures = list()
        
        # Aggregators
        self._sourcesAggregator = WifiObjectsAggregator(self._captures, WifiCapture.getSourcesCount)
        self._packetsAggregator = WifiObjectsAggregator(self._captures, WifiCapture.getPacketsCount)
        self._macAggregator = WifiObjectsAggregator(self._captures, WifiCapture.getMacList)
        
        # Selectors
        self._statsSelector = WifiObjectsSelector(self._captures, attributeName="stats")
        self._dateSelector = WifiObjectsSelector(self._captures, attributeName="start")

    def addCapture(self, wifiCapture):
        if type(wifiCapture) is not WifiCapture:
            raise WifiObjectRawTypeError("WifiCapture", type(wifiCapture))

        self._captures.append(wifiCapture)
    
    def hasCapture(self):
        return len(self._captures) > 0
    
    def isEmpty(self):
        for capture in self._captures:
            if not capture.isEmpty():
                return False
        return True
        
    def getLastCapture(self, datetimeOrdered = False):
        if len(self._captures) < 1:
            return None
        
        if datetimeOrdered:
            return max(self._captures, key=lambda capture : capture.start)
        
        return self._captures[-1]
        
    def getFirstCapture(self, datetimeOrdered = False):
        if len(self._captures) < 1:
            return None
        
        if datetimeOrdered:
            return min(self._captures, key=lambda capture : capture.start)

        return self._captures[0]
        
    def extract(self,
                fromDatetime=None, toDatetime=None,
                minSources=None, maxSources=None,
                minPackets=None, maxPackets=None,
                minStats=None, maxStats=None,
                macList=list(), forbiddenList=True):
        """Return a copy of the current WifiAcquisition which only contains WifiCaptures that satisfy all criteria"""
        acquisition = copy.deepcopy(self)# Needed to keep internal types (Raw vs not Raw)
        
        Debug("""Select captures with date between "{}" and "{}":""".format(fromDatetime, toDatetime))
        self._RemoveAllExcept(acquisition._captures, acquisition._dateSelector.selectBetween(fromDatetime, toDatetime))
        
        Debug("""Select captures with MAC addresses {} in "{}":""".format("not" if forbiddenList else "", macList))
        self._RemoveAllExcept(acquisition._captures, acquisition._macAggregator.selectIn(macList, forbiddenList))
        
        Debug("""Select captures with sources between "{}" and "{}":""".format(minSources, maxSources))
        self._RemoveAllExcept(acquisition._captures, acquisition._sourcesAggregator.selectBetween(minSources, maxSources))
        
        Debug("""Select captures with packets between "{}" and "{}":""".format(minPackets, maxPackets))
        self._RemoveAllExcept(acquisition._captures, acquisition._packetsAggregator.selectBetween(minPackets, maxPackets))
        
        Debug("""Select captures with stats between "{}" and "{}":""".format(minStats, maxStats))
        self._RemoveAllExcept(acquisition._captures, acquisition._statsSelector.selectBetween(minStats, maxStats))
        
        return acquisition
                    
    def getCapturesCount(self):
        return len(self._captures)
    
    def getPacketsCount(self):
        return self._packetsAggregator.getSum()
    
    def getSourcesCount(self):
        return self._sourcesAggregator.getSum()
    
    def getMacList(self, removeDuplicates=False):
        return WifiObjectsAggregator(self._captures, WifiCapture.getMacList).aggregateList(removeDuplicates)

    def __hash__(self):
        result = self.getCapturesCount() * 1
        
        pckCnt = self.getPacketsCount()
        result += (pckCnt * 10) if pckCnt else 0
        
        srcCnt = self.getSourcesCount()
        result += (srcCnt * 100) if srcCnt else 0

        result += (self.duration * 1000) if self.duration else 0
        
        result += (self.spacing * 10000) if self.spacing else 0

        return result

    def __eq__(self, other): 
        return (
            self._captures == sorted(other._captures) and
            self.duration == other.duration and
            self.spacing == other.spacing and
            self.author == other.author and
            self.label == other.label
        )
        
    def __repr__(self):
        return """<{} duration="{}" spacing="{}" captures="{}">""".format(self.__class__.__name__, self.duration, self.spacing, self.getCapturesCount())
    
class WifiRawAcquisition(WifiAcquisition, WifiRawObject):
    """Represents a sequence of WifiCapture"""
    
    def __init__(self, duration, spacing, wifiAuthor, wifiLabel, wifiStats=None):
        WifiAcquisition.__init__(self, duration, spacing, wifiAuthor, wifiLabel, wifiStats)
        self.stats.getDataMethod = self.getAllDistances
        
    def addCapture(self, wifiRawCapture):
        if type(wifiRawCapture) is not WifiRawCapture:
            raise WifiObjectRawTypeError("WifiRawCapture", type(wifiRawCapture))
        
        self._captures.append(wifiRawCapture)
        
    def getAllDistances(self):
        """Returns a list of all distances contained in this capture (from each source)"""
        return WifiObjectsAggregator(self._captures, WifiRawCapture.getAllDistances).aggregateList()
        
        
class WifiCapture(WifiObject):
    """Represents a set of WifiSource"""
    
    def __init__(self, startDateTime = None, wifiStats = None):
        if type(startDateTime) not in [datetime, NoneType]:
            raise TypeError("WifiCapture requires type 'datetime.datetime' for parameter startDateTime. (Given type : '{}')".format(type(startDateTime)))

        self.start = startDateTime
        self.stats = wifiStats if wifiStats else WifiStats()
        self._sources = list()
        
        # Aggregators
        self._packetsAggregator = WifiObjectsAggregator(self._sources, attributeName="packets")
        self._macAggregator = WifiObjectsAggregator(self._sources, attributeName="macAddress")
        
        # Selectors
        self._statsSelector = WifiObjectsSelector(self._sources, attributeName="stats")
        
    def addSource(self, wifiSource):
        if type(wifiSource) is not WifiSource:
            raise WifiObjectRawTypeError("WifiSource", type(wifiSource))

        self._sources.append(wifiSource)
        
    def isEmpty(self):
        for source in self._sources:
            if not source.isEmpty():
                return False
        return True
    
    def extract(self,
                minPackets=None, maxPackets=None,
                minStats=None, maxStats=None,
                macList=list(), forbiddenList=True):
        """Return a copy of the current WifiCapture which only contains WifiSources that satisfy all criteria"""
        capture = copy.deepcopy(self)# Needed to keep internal types (Raw vs not Raw)
        
        Debug("""Select sources with MAC addresses {} in "{}":""".format("not" if forbiddenList else "", macList))
        self._RemoveAllExcept(capture._sources, capture._macAggregator.selectIn(macList, forbiddenList))
        
        Debug("""Select sources with packets between "{}" and "{}":""".format(minPackets, maxPackets))
        self._RemoveAllExcept(capture._sources, capture._packetsAggregator.selectBetween(minPackets, maxPackets))
        
        Debug("""Select sources with stats between "{}" and "{}":""".format(minStats, maxStats))
        self._RemoveAllExcept(capture._sources, capture._statsSelector.selectBetween(minStats, maxStats))
        
        return capture
        
    def getSourcesCount(self):
        return len(self._sources)
        
    def getPacketsCount(self):
        return self._packetsAggregator.getSum()
    
    def getMacList(self, removeDuplicates=False):
        return self._macAggregator.aggregateList(removeDuplicates)
    
    def __eq__(self, other): 
        return (
            self._sources == sorted(other._sources) and
            self.start == other.start
        )
        
    def __repr__(self):
        return """<{} start="{}" sources="{}">""".format(self.__class__.__name__, self.start, self.getSourcesCount())
    
class WifiRawCapture(WifiCapture, WifiRawObject):
    """Represents a set of WifiSource"""
    
    def __init__(self, startDateTime, wifiStats = None):
        WifiCapture.__init__(self, startDateTime, wifiStats)
        self.stats.getDataMethod = self.getAllDistances
    
    def addSource(self, wifiRawSource):
        if type(wifiRawSource) is not WifiRawSource:
            raise WifiObjectRawTypeError("WifiRawSource", type(wifiRawSource))

        self._sources.append(wifiRawSource)
        
    def getAllDistances(self):
        return WifiObjectsAggregator(self._sources, attributeName="_distances").aggregateList()
    
        
class WifiSource(WifiObject):
    """Source wifi defined by its madAdress, an amount of packets received, and statistics about distances evaluated for each packet between emitter and receiver."""
    
    def __init__(self, macAddress, nbPackets = None, wifiStats = None):
        self.macAddress = macAddress
        self.packets = nbPackets
        self.stats = wifiStats if wifiStats else WifiStats()
        
    def isEmpty(self):
        return self.macAddress == None
    
    def __eq__(self, other): 
        return (
            self.macAddress == other.macAddress and
            self.packets == other.packets
        )
        
    def __repr__(self):
        return """<{} mac="{}" packets="{}">""".format(self.__class__.__name__, self.macAddress, self.packets)
    
class WifiRawSource(WifiSource, WifiRawObject):
    """WifiSource with all its evaluated distances between emitter and receiver for each packet received"""
    
    def __init__(self, macAddress, distances = list()):
        WifiSource.__init__(self, macAddress=macAddress, nbPackets=0)
        self._distances = distances
        self.packets = len(self._distances)
        self.stats.getDataMethod = self.getDistances
        
    def getDistances(self):
        return self._distances
    
    def addDistance(self, wifiDistance):
        self._distances.append(float(wifiDistance))
        self.packets += 1    
    
class WifiAcquisitionAuthor(object):
    """Represents the author of a WifiAcquisition"""
    
    def __init__(self, name = None, macAddress = None):
        self.name = str(name) if name else None
        self.macAddress = macAddress
        
    def extractionComparision(self, other):
        return (
            (not other.name or self.name == other.name) and
            (not other.macAddress or self.macAddress == other.macAddress)
        )
        
    def __eq__(self, other):
        return (
            (not self.name or not other.name or self.name == other.name) and
            (not self.macAddress or not other.macAddress or self.macAddress == other.macAddress)
        )

    def __repr__(self):
        return """<{} name="{}" mac="{}">""".format(self.__class__.__name__, self.name, self.macAddress)
    
class WifiAcquisitionLabel(object):
    """Represents a label of a WifiAcquisition"""
    
    def __init__(self, building=None, floor=None, room=None, colleagues=list(), comment=None):
        self.building = str(building) if building else None
        self.floor = int(floor) if floor else None
        self.room = str(room) if room else None
        self.colleagues = colleagues
        self.comment = str(comment) if comment else None
    
    def extractionComparision(self, other):
        return (
            (not other.building or self.building == other.building) and
            (not other.floor or self.floor == other.floor) and
            (not other.room or self.room == other.room) and
            (not other.colleagues or (self.colleagues and set(other.colleagues).issubset(set(self.colleagues)))) and
            (not other.comment or self.comment == other.comment)
        )
        
    def __eq__(self, other):
        return (
            (self.building == other.building) and
            (self.floor == other.floor) and
            (self.room == other.room) and
            ((not self.colleagues and not other.colleagues) or (self.colleagues and other.colleagues and sorted(self.colleagues) == sorted(other.colleagues))) and
            (self.comment == other.comment)
        )
        
    def __repr__(self):
        return """<{} building="{}" floor="{}" room="{}" colleagues="{}" comment="{}">""".format(self.__class__.__name__, self.building, self.floor, self.room, self.colleagues, self.comment)

if __name__ =='__main__':
    
    import WifiLog as wo
    content = wo.WifiLogFileLoader("test/2015-01-06.rawlog", verbosityLevel = wo.WifiLogFile.VerbosityLevel.No).load()
    print wo.WifiLogFileDumper().dump(content, mode=wo.WifiLogFile.Type.Xml)
    #content.extract(minDuration, maxDuration, minSpacing, maxSpacing, minCaptures, maxCaptures, author, label, fromDatetime, toDatetime, minSources, maxSources, minPackets, maxPackets, minStats, maxStats, macList, forbiddenList)
    print "ACQUISITION"
    print content.getFirstAcquisition().getMacList()
    print content.getFirstAcquisition().getAllDistances()
    print "CAPTURE"
    print content.getFirstAcquisition().getFirstCapture().getMacList()
    print content.getFirstAcquisition().getFirstCapture().getAllDistances()
    print "OTHER"
    print content.getFirstAcquisition().extract(macList=["50:17:ff:dd:95:af"], forbiddenList=True).getMacList()
    print content.getFirstAcquisition().getFirstCapture().extract(minStats=WifiStats(maxDistance=.8)).getMacList()
    #content = wo.WifiLogFileLoader("20150112.rawlog", verbosityLevel = wo.WifiLogFile.VERBOSITY_NONE).load()
    content.removeDuplicates()
    print wo.WifiLogFileDumper().dump(content.extract(author=WifiAcquisitionAuthor("Lefebvre Leo")), mode=wo.WifiLogFile.Type.Xml)
    