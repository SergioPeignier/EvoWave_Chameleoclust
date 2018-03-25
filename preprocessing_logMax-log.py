#-*- coding: utf-8 -*-
# 12/02/2016 : This script import evowave aggregated data (csv), build class labels, apply a metric, and export to hdf5 pandas table.

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
random.seed(0)

datas = pd.read_hdf("wifi_stream.h5")

#convert data
data_dims = data.iloc[:,18:-2]
data_dims = data_dims.astype(np.float64)
med_max   = data_dims.max().median()
data_dims[data_dims>med_max] = med_max
data_dims = np.log10(med_max) - np.log10(data_dims.values)
data.iloc[:,18:-2] = data_dims
data.iloc[:,18:-2] = data.iloc[:,18:-2].fillna(0)

#plot histogram
"""
res = []
for i in range(256):
	res += [a for a in data.iloc[:,18+i] if not np.isnan(a)]
plt.hist(res)
plt.show()
"""
#export
ouput_path = 'wifi_preprocessed.h5'
ouput_key = 'logMax_logx'
data.to_hdf(ouput_path,ouput_key)
