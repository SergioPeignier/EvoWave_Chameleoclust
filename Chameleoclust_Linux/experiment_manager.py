import os
import fnmatch
import pandas as pd
import scipy.sparse
from multiprocessing import Pool
from Chameleoclust import *
random.seed(0)
np.random.seed(0)


def simulation(dic_params):
	chameleoclustparams = dic_params["chameleoclustparams"]
	runparams = dic_params["runparams"]
	pop_genome = dic_params["set_pop"]
	c = Chameleoclust(**chameleoclustparams)
	if len(pop_genome.keys()):
		c.simulation.esetpopulationgenome(pop_genome["genome"])
	c.erun(**runparams)
	return 1

def executeparallelsimulations(arguments,f = simulation):
    pool = Pool()
    pool.map(f,arguments)
    pool.close()
    pool.join()
    
def executeseriesofparallelexperiments(arguments,workers=4, f = simulation):
	i = 0
	while True:
		i_n = min(i+workers, len(arguments))
		executeparallelsimulations(arguments[i:i_n], f)
		i = i_n
		if i >= len(arguments): break
	
    
    

def find_related_files_in_folder(folder_path,pattern):
	matches = []
	for root, dirnames, filenames in os.walk(folder_path):
		for fold in fnmatch.filter(filenames, pattern):
			matches.append(os.path.join(root, fold))
	return matches

def computematrixensembleclustering(folder_path,pattern):
	files_to_use = find_related_files_in_folder(folder_path,pattern)
	if len(files_to_use) == 0:
		return []
	dataset = load_object_from_file(files_to_use[0])
	nb_objects = len(dataset.index)
	frequency  = [[0 for i in xrange(nb_objects)] for i in xrange(nb_objects)]
	for file in files_to_use:
		exp = load_object_from_file(file)
		object_clusters = list(exp["cluster_found"])
		for o_i in xrange(nb_objects):
			for o_j in xrange(nb_objects):
				if object_clusters[o_i] == object_clusters[o_j]:
					frequency[o_i][o_j] += 1
	for o_i in xrange(nb_objects):
		for o_j in xrange(nb_objects): 
			frequency[o_i][o_j] *= 1./len(files_to_use)
	return frequency

def computematrixdimensionsensembleclustering(folder_path,pattern_phenotype):
	files_phenotype = find_related_files_in_folder(folder_path,pattern_phenotype)
	if len(files_phenotype) == 0 :
		return []
	pheno = load_object_from_file(files_phenotype[0])
	nb_dims    = len(pheno.columns)
	matrix_subspaces = [[0 for k in xrange(nb_dims)] for j in xrange(nb_dims)]
	nb_clusters = 0
	for i in xrange(len(files_phenotype)):
		pheno = load_object_from_file(files_phenotype[i]).notnull()
		for row_index, row in pheno.iterrows():
			truedims = [i for i, x in enumerate(list(row)) if x]
			for td1 in truedims:
				for td2 in truedims:
					matrix_subspaces[td1][td2] += 1
			nb_clusters += 1
	return [[j * 1./nb_clusters for j in i] for i in matrix_subspaces]
			
	
def computematrixsubspacesensembleclustering(folder_path,pattern_clustering,pattern_phenotype):
	files_phenotype = find_related_files_in_folder(folder_path,pattern_phenotype)
	files_clustering = find_related_files_in_folder(folder_path,pattern_clustering)
	if len(files_phenotype) == 0 or len(files_clustering) == 0:
		return []
	pheno = load_object_from_file(files_phenotype[0])
	clus  = load_object_from_file(files_clustering[0])
	nb_objects = len(clus.index)
	nb_dims    = len(pheno.columns)
	matrix_subspaces = [[0 for k in xrange(nb_dims)] for j in xrange(nb_objects)]
	for i in xrange(len(files_phenotype)):
		clust = load_object_from_file(files_clustering[i])
		pheno = load_object_from_file(files_phenotype[i]).notnull()
		object_clusters = list(clust["cluster_found"])
		for o_i in xrange(nb_objects):
			subspace = list(pheno.ix[object_clusters[o_i]])
			for o_j in xrange(nb_dims):
				if subspace[o_j]:
					matrix_subspaces[o_i][o_j] += 1
	"""
	for o_i in xrange(nb_objects):
		for o_j in xrange(nb_dims):
			matrix_subspaces[o_i][o_j] *= 100.0 /  len(files_phenotype)
			matrix_subspaces[o_i][o_j] = int(matrix_subspaces[o_i][o_j])
	"""
	return matrix_subspaces

def floodandcomputenbconnectedcomponents(matrix):
	nb_connected_components = []
	for i in xrange(matrix[0][0]):
		loc_matrix = [[max(c - i,0) for c in r] for r in matrix]
		nb_connected_components.append(scipy.sparse.csgraph.connected_components(loc_matrix,directed=False)[0])
	return nb_connected_components
	
def computeweightedlistdims(folder_path,pattern_phenotype):
	files_phenotype = find_related_files_in_folder(folder_path,pattern_phenotype)
	if len(files_phenotype) == 0 :
		return {}
	res = {}
	for i in xrange(len(files_phenotype)):
		pheno = load_object_from_file(files_phenotype[i]).notnull()
		nb_clus = len(pheno.index)
		nb_dims = len(pheno.columns)
		for o_i in xrange(nb_clus):
			subspace = list(pheno.iloc[o_i])
			for o_j in xrange(nb_dims):
				if subspace[o_j] and o_j in res.keys():
					res[o_j] += 1
				if subspace[o_j] and o_j not in res.keys():
					res[o_j] = 1
	return res

def concatenatedffilesfromdifferentexperiments(folder_path,pattern):
	files_to_use = find_related_files_in_folder(folder_path,pattern)
	if len(files_to_use) > 0:
		dataset = load_object_from_file(files_to_use[0])
		i = 1
		while i < len(files_to_use):
			dataset = pd.concat((dataset,load_object_from_file(files_to_use[i])))
			i += 1
		return dataset
		
	
def computestatisticsfromdifferentexperiments(folder_path,pattern):
	files_to_use = find_related_files_in_folder(folder_path,pattern)
	if len(files_to_use) > 0:
		dataset = load_object_from_file(files_to_use[0])
		stats = {}
		i = 1
		if isinstance(dataset, dict):
			while i < len(files_to_use):
				for k in dataset:
					stats[k] = {}
					dataset_i = load_object_from_file(files_to_use[i])	
					if isinstance(dataset[k], dict):
						for j in dataset[k]:
							dataset[k][j] = pd.concat((dataset[k][j],dataset_i[k][j]))
					else:
						dataset_i = load_object_from_file(files_to_use[i])			
						for k in dataset:
							dataset[k] = pd.concat((dataset[k],dataset_i[k]))
				i += 1
		else:
			while i < len(files_to_use):
				dataset = pd.concat((dataset,load_object_from_file(files_to_use[i])))
				i += 1
			stats = dataset.groupby(dataset.index)
			return stats
		
		if not isinstance(dataset[dataset.keys()[0]], dict) :
			for k in dataset:
				stats[k] = dataset[k].groupby(dataset[k].index)
		else:
			for k in dataset:
				for j in dataset[k]:
					stats[k][j] = dataset[k][j].groupby(dataset[k][j].index)
		return stats
		
def computetimesfromdifferentexperiments(folder_path,pattern,sep="\t",header=None):
        files_to_use = find_related_files_in_folder(folder_path,pattern)
        if len(files_to_use) > 0:
                dataset = pd.DataFrame.from_csv(files_to_use[0],sep=sep,header=header)
                stats = {}
                i = 1
                while i < len(files_to_use):
                        dataset_i = pd.DataFrame.from_csv(files_to_use[i],sep=sep,header=header)
                        dataset_i.columns = ["UTIME","STIME","MAXRSS"]
                        dataset = pd.concat((dataset,dataset_i))
                        i += 1
                stats = dataset.groupby(dataset.index)
                return stats