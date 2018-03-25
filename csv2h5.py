#-*- coding: utf-8 -*-
# 12/02/2016 : This script import evowave aggregated data (csv), build class labels, apply a metric, and export to hdf5 pandas table.

import pandas as pd
import numpy as np
import random

stream = 1
 #import datas
input_path = "../20160208_EVOWAVE_BASE1/collect"+str(stream)+"/finalGlobal"+str(stream)+".csv"
data = pd.read_csv(input_path)
data = data[data.building != 'building'] #kick multiple headers
data.index = range(len(data.index))
#Make class labels
from collections import Counter
counter_fines = Counter(zip(data.building,data.room,data.colleagues))
etiquette_fines = {key:i for i,key in enumerate(counter_fines.keys())}
counter_laches = Counter(zip(data.building,data.room))
etiquette_laches = {key:i for i,key in enumerate(counter_laches.keys())}
data['classes'] = [str(etiquette_fines[tpl])+"_"+str(stream) for tpl in zip(data.building,data.room,data.colleagues)]
#export
ouput_path = "wifi_stream.h5"
ouput_key = 'wifi'
data.to_hdf(ouput_path,ouput_key)
