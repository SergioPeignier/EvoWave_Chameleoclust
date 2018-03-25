import sys
sys.path.append("Chameleoclust_OSx")
from Chameleoclust import Chameleoclust
import pandas as pd
import numpy as np
from numpy import array
import matplotlib.pyplot as plt

def local_test(X,Y,cham,i,iterations,save_local_results =False):
    print i
    cham._fit_local(X = X.iloc[i-sliding_sample_size:i,:],nb_iterations = iterations)
    cham.compute_local_stats()
    y = cham._predict_local()
    y_t = array(Y[i-sliding_sample_size:i],dtype=object)
    phenotype = cham._get_best_individual_phenotype()
    confusion_matrix = build_confusion_table(y,y_t)
    if save_local_results:
        confusion_matrix_to_save = modify_confusion_table_for_plots(confusion_matrix, i)
        confusion_matrix_to_save.to_hdf(confusion_matrix_name,key(seed,i))
        phenotype.to_hdf(phenotype_file_name,key(seed,i))
        model = cham._get_best_individual_model()
        model.to_hdf(model_file_name,key(seed,i))    
    return confusion_matrix,phenotype

#load data
dataset = "datasets/pendigits.h5"
data = pd.read_hdf(dataset)
precision = 1000
vecfunc = np.vectorize(lambda x: int(x*precision) if not np.isnan(x) else np.nan)
X = data.iloc[:,0:-1]
X = (X - X.mean()) / X.std()
X = vecfunc(X)
Y = list(data.iloc[:,-1])
print X
#build Chameleoclust
seed = 0
cham = Chameleoclust(cmax=10, seed = seed)
#train and predict
cham.fit(pd.DataFrame(X))
results = pd.DataFrame()
results["clusters"] = np.array(cham.predict(pd.DataFrame(X)))
results["classes"]  = np.array(Y)
phenotype = cham._get_best_individual_phenotype()
print 
print "Number of clusters:"
print len(np.unique(results["clusters"]))

print "Number of core-points:"
print len(phenotype)

print "Core-points coordinates:"
print phenotype

print "Confusion matrix:"
print pd.crosstab(results["classes"],results["clusters"])
