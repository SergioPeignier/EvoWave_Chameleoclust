import sys
sys.path.append("Chameleoclust_OSx")
from Chameleoclust import Chameleoclust
from modevoevoPY.evaluationfunctions import Compute_only_entropy,Compute_only_accuracy,Compute_only_F1,Compute_only_CE,mapped,ValidityClusterChecking,egetavgdim,egetnbcorepoints
import pandas as pd
import numpy as np
from numpy import array
from scipy import stats
import random
import matplotlib.pyplot as plt

phenotype_file_name = "phenotype.h5"
model_file_name = "model.h5"
confusion_matrix_name = "confusion_matrix.h5"
aggregated_data_name  = "aggregated_data.h5"
phenotype_stats_name  = "clustering_stats.h5"
key = lambda seed, generation: "seed"+str(seed)+"/"+"generation"+str(generation)


def build_confusion_table(y_clusters,y_classes):
    confusion_table = pd.crosstab(y_clusters,y_classes)
    confusion_table.index.name = "cluster label"
    confusion_table.columns.name = "class label"
    return confusion_table

def modify_confusion_table_for_plots(confusion_table,generation=None):
    order = confusion_table.sum(1).sort_values(ascending=0).index
    confusion_table = confusion_table.loc[order]
    if generation is not None:
        confusion_table.index = [str(generation)+"-"+str(cluster_id) for cluster_id in confusion_table.index]
    return confusion_table

def compute_statistics(confusion_matrix,phenotype):
    confusion_matrix = confusion_matrix.T
    found_clusters_effective = confusion_matrix.sum(0)
    valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity=0)
    mapped_clusters = mapped(confusion_matrix)
    #clusters evaluation
    e = Compute_only_entropy(confusion_matrix, valid_clusters)
    a = Compute_only_accuracy(confusion_matrix, valid_clusters, found_clusters_effective)
    f1 = Compute_only_F1(confusion_matrix, valid_clusters, mapped_clusters)
    ce = Compute_only_CE(confusion_matrix, found_clusters_effective, valid_clusters)
    #phenotype description
    clusters = phenotype.loc[list(confusion_matrix.loc[:,valid_clusters].columns),:]
    avg_core_points_dim = egetavgdim(phenotype)
    nb_core_points      = egetnbcorepoints(phenotype)
    avg_clusters_dim = egetavgdim(clusters)
    nb_clusters      = egetnbcorepoints(clusters)   
    return [e,a,f1,ce,nb_clusters,avg_clusters_dim,nb_core_points,avg_core_points_dim]

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
    print model, phenotype 
    return confusion_matrix,phenotype

#load data
dataset = "datasets/wifi_filter_120_LogMax-LogX_final_classes.h5"
data = pd.read_hdf(dataset)
precision = 10000
vecfunc = np.vectorize(lambda x: int(x*precision) if not np.isnan(x) else np.nan)
data.iloc[:,0:-1]  = vecfunc(data.iloc[:,0:-1])
X = data.iloc[:,0:-1]
Y = list(data['classes'])
#build Chameleoclust
seed =2
sliding_sample_size = 100
initital_training   = 100
training_generations = 10
distances_stream_to_zero = []
stats_clustering = pd.DataFrame(columns=["entropy","accuracy","F1","CE","NumClusters","AvgClusterDim","NumCorePoints","AvgCorePointDim"])
cham = Chameleoclust(cmax=10, sliding_sample_size=sliding_sample_size, seed = seed)
#train and predict
confusion_matrix,phenotype = local_test(X = X,Y = Y,cham = cham,i = sliding_sample_size,iterations = initital_training, save_local_results = True)
stats_clustering.loc[100] = compute_statistics(confusion_matrix,phenotype)
distances_stream_to_zero.append(X.iloc[0:sliding_sample_size,:].abs().sum().sum()*1./sliding_sample_size)
for i in range(sliding_sample_size+1,110):#,len(X)):
    confusion_matrix,phenotype = local_test(X = X,Y = Y,cham = cham,i = i,iterations = training_generations, save_local_results = True)
    stats_clustering.loc[i] = compute_statistics(confusion_matrix,phenotype)
    distances_stream_to_zero.append(distances_stream_to_zero[-1] + (X.iloc[i,:].abs().sum() - X.iloc[i-sliding_sample_size,:].abs().sum()) * 1./sliding_sample_size)
    print stats_clustering.loc[i] 
#save results
full_stats = cham._get_stats()
full_stats["distance_to_0"] = distances_stream_to_zero
full_stats["distance_to_0"] = full_stats["distance_to_0"]*-1

full_stats.to_hdf(aggregated_data_name,"seed"+str(seed))
stats_clustering.to_hdf(phenotype_stats_name,"seed"+str(seed))
