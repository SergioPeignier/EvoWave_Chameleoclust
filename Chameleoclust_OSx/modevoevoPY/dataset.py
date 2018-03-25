from definitions		import *
import numpy as np
from scipy import stats
import random
import pylab
import pandas as pd
import numpy as np
import pdb

def clusters_column_from_clusters_membership(cluster_membership):
	nb_points = 0
	for clus in cluster_membership:
		nb_points += len(clus)
	clusters_column = range(nb_points)
	for i,clus in enumerate(cluster_membership):
		for j,point in enumerate(clus):
			clusters_column[point] = i
	return clusters_column
	
class true_file:
	def __init__(self,file_name,dim_i="",noise_clusters = ['-1']):
		self.file_name = file_name
		self.noise_clusters = noise_clusters
		f = file(file_name,"r")
		dim = f.readline().split(";")[0].split("=")
		if len(dim) > 1:
			dim = int(dim[1])
		else:
			dim = dim_i
		self.clusters_subspaces = []
		self.clusters_membership = []
		for line in f:
			l = map(int,line.split(" "))
			self.clusters_subspaces.append(l[0:dim])
			self.clusters_membership.append(l[dim+1:len(l)])
		f.close()
		self.clusters_subspaces = pd.DataFrame(self.clusters_subspaces)
		self.true_arff_dico = {}
		self.arff_true_dico = {}
		self.df_true_dico   = {}
	
	def generate_dico_arff_true_clusters(self, arff_membership):
		self.true_arff_dico = {}
		self.arff_true_dico = {}
		for i,hidden_cluster_true in enumerate(self.clusters_membership):
			self.true_arff_dico[i] = []
			for point in hidden_cluster_true:
				if arff_membership[point] not in self.true_arff_dico[i]:
					self.true_arff_dico[i].append(arff_membership[point])
		self.arff_true_dico = self.invert_dict(self.true_arff_dico)

		return self.check_arff_true_dico(arff_membership)
	
	
	def invert_dict(self,dico):
		inv_map = {}
		for k, v in dico.iteritems():
			for vi in v:
				inv_map[vi] = inv_map.get(vi, [])
				inv_map[vi].append(k)
		return inv_map
	
	def check_arff_true_dico(self,arff_membership):
		for i,arff_i in enumerate(arff_membership):
			if arff_i not in self.noise_clusters:
				for true_cluster in self.arff_true_dico[arff_i]:
					if i not in self.clusters_membership[true_cluster]:
						return 0
		return 1
	
	def compute_df_true_dico(self,arff_df_dico):
		self.df_true_dico = {}
		for arff_i in self.arff_true_dico.keys():
			self.df_true_dico[arff_df_dico[arff_i]] = self.arff_true_dico[arff_i]
	
	
class data:
	def __init__(self,file_name,
                 hidden_cluster=[],
                 found_cluster=[],
                 useless=[],
                 index=False,
                 header= False,
                 comments = "#",
                 extension = "csv",
                 sep = ",",
                 make_convertion_dico=True):
		if extension == "csv":
			self.df = pd.DataFrame().from_csv(file_name, header=header, index_col=index, sep= sep)
 		if extension == "arff":
 			self.df = pd.DataFrame(np.genfromtxt(file_name , skip_header=header,delimiter=sep,dtype=str,comments = comments))
		if extension == "hdf": 
			self.df = pd.read_hdf(file_name)
		self.useless     = useless
		self.hidden_cluster = hidden_cluster
		self.found_cluster = found_cluster
		self.features = [e for e in list(self.df.columns.values) if e not in self.hidden_cluster and e not in self.found_cluster and e not in self.useless]
		self.__cast()
		self.__drop_useless()
		self.conversion_dico = {}
		self.shuffled_index = list(self.df.index)
		if make_convertion_dico:
			if self.hidden_cluster != []:
				self.hidden_cluster_membership = list(self.df[self.hidden_cluster[0]])[:] 
			if self.found_cluster != []:
				self.convert_cluster_membership_to_numerical(self.found_cluster)
			if self.hidden_cluster != []:
				self.convert_cluster_membership_to_numerical(self.hidden_cluster)
	
	def __drop_useless(self):
		for feature in self.useless:
			self.df = self.df.drop(feature, 1)

	def __cast(self):
		for feature in self.features:
			self.df[feature] = self.df[feature].astype(float)
			
	def convert_cluster_membership_to_numerical(self,cluster_col):
		self.conversion_cluster_membership_dico(cluster_col)
		for i in self.df.index:
			self.df.loc[i,cluster_col[0]] = self.conversion_dico[cluster_col[0]][self.df.loc[i,cluster_col[0]]]
		self.df[cluster_col[0]] = self.df[cluster_col[0]].astype(int)
		
	def conversion_cluster_membership_dico(self,cluster_col):
		df_col = self.df[cluster_col[0]]
		self.conversion_dico[cluster_col[0]] = {k:i for i,k in enumerate(df_col.unique())}

	def standardize_table(self):
		self.df[self.features] = stats.zscore(self.df[self.features])
		self.df[self.features] = self.df[self.features].fillna(0)

 	def center_table(self):
		self.df[self.features] = pylab.demean(self.df[self.features], axis=1)

	def normalize_table_by_std(self):
		self.df[self.features] = self.df[self.features]/np.std(self.df[self.features],axis=0)
		self.df[self.features] = self.df[self.features].fillna(0)

	def multiply_table_by_precision_and_convert_to_int(self,precision=STD_PRECISION):
		vecfunc = np.vectorize(lambda x: int(x*precision) if not np.isnan(x) else np.nan)
		self.df[self.features] = vecfunc(self.df[self.features])

	def shift_table(self,shift):
		self.df[self.features] = self.df[self.features].values + shift

	def logscale_table(self):
		table = self.df[self.features]
		self.df[self.features]  = np.log(table*(table>0)+(table<=0))-np.log(-table*(table<0)+(table>=0))

	def table_to_data(self,features = ""):
		if features == "":
			features = self.features
		self.data = []
		i = 0
		while i < len(self.df):
			self.data.append([])
			if self.hidden_cluster == []:self.data[i].append(-1)
			else: self.data[i].append(self.df.loc[self.df.index[i],self.hidden_cluster[0]])
			if self.found_cluster == []:self.data[i].append(-1)
			else: self.data[i].append(self.df.loc[self.df.index[i],self.found_cluster[0]])
			self.data[i] += [[j,e] for j,e in enumerate(self.df.loc[self.df.index[i],features]) if not np.isnan(e)]
			i+=1
		return self.data

	def table_to_subobjects_stream(self):
		self.stream = []
		i = 0
		while i < len(self.df):
			for j,e in enumerate(self.df.loc[i,self.features]):
				self.stream.append([i,j,e])
			i+=1
		return self.stream

	def set_column_values(self,values,column):
		for i,e in enumerate(values):
			self.df[column].values[i] = e
	
	def shuffle(self,shuffled_index = []):
		if shuffled_index == []:
			self.shuffled_index = range(len(self.df.index))
			random.shuffle(self.shuffled_index)
		else:
			self.shuffled_index  = shuffled_index
		new_index = [list(self.df.index)[i] for i in self.shuffled_index]
		self.df = self.df.ix[new_index]
		self.df = self.df.reset_index(drop=True)
		self.data = [self.data[i] for i in self.shuffled_index]

	def reorder_dataframe(self,dataframe):
		dataframe.index = self.shuffled_index
		return dataframe.sort()