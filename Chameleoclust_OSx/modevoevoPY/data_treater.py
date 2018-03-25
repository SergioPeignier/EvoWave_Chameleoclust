from definitions		import *
import numpy as np
from scipy import stats
import random
import pylab
class data:
	def __init__(self):
		self.table 		 = []
		self.data			 = []
		self.header			= ""
		self.clusters		= []
		self.labels			= []
		self.clusters_dico = {}
		self.labels_dico	 = {}
	def open(self,file_name,sep=" ",header=0,missing_values="na",c_cluster="",c_label="",comments="#",c_not_to_use=[]):
		not_data_cols = []
		self.table = np.genfromtxt(file_name , skip_header=header,delimiter=sep, filling_values=(missing_values,0),dtype=str,comments = comments)
		if c_cluster != "" or c_label != "" or c_label == "last":
			if c_cluster != "":
				self.clusters = np.genfromtxt(file_name , skip_header=header,delimiter=sep,dtype=str,usecols=(c_cluster),comments = comments).tolist()
				not_data_cols.append(c_cluster)
			if c_label	 != "" and c_label	 != "last":
				self.labels = np.genfromtxt(file_name , skip_header=header,delimiter=sep,dtype=str,usecols=(c_label),comments = comments).tolist()
				not_data_cols.append(c_label)
			if c_label == "last":
				self.labels = self.table[:,len(self.table[1,])-1]
				not_data_cols.append(len(self.table[1,])-1)

		cols			 = list(set(range(len(self.table[0,:])))-set(not_data_cols+c_not_to_use))
		self.table = self.table[:,cols]
		vecfunc = np.vectorize(lambda x: float(x))
		self.table = vecfunc(self.table)
		return self.table.tolist()

	def denoise(self,noise_class):
		denoised_table = []
		denoised_labels = []
		for i,c in enumerate(self.labels):
			if c not in noise_class:
				denoised_labels.append(i)
		self.table  = self.table[denoised_labels,:]
		self.labels = self.labels[denoised_labels]
		
	def find_labels_or_cluster_dico(self,clus_or_label):
		dico = {}
		i = 0
		while i < len(clus_or_label):
			if clus_or_label[i] not	in dico.keys():
				dico[clus_or_label[i]]=len(dico.keys())
			i += 1
		return dico

	def generate_clus_or_label_numeric_vector(self,clus_or_label,dico):
		vec = []
		for c in clus_or_label:
			vec.append(dico[c])
		return vec

	def standardize_table(self):
		self.table = stats.zscore(self.table)
		self.table	= np.nan_to_num(self.table)
		return self.table.tolist()

 	def center_table(self):
		return pylab.demean(self.table, axis=1)

	def normalize_table_by_std(self):
		self.table = self.table/np.std(self.table,axis=0)
		self.table	= np.nan_to_num(self.table)
		return self.table.tolist()

	def multiply_table_by_precision_and_convert_to_int(self,precision=1000):
		vecfunc = np.vectorize(lambda x: int(x*precision))
		self.table = vecfunc(self.table)
		return self.table.tolist()

	def shift_table(self,shift):
		self.table = np.array(self.table) + shift
		return self.table.tolist()

	def logscale_table(self):
		self.table = np.array(self.table)
		self.table = np.log(self.table*(self.table>0)+(self.table<=0))-np.log(-self.table*(self.table<0)+(self.table>=0))

	def table_to_data(self):
		numeric_cluster = []
		numeric_label	 = []
		if self.clusters != []:
			self.clusters_dico = self.find_labels_or_cluster_dico(self.clusters)
			numeric_cluster		= self.generate_clus_or_label_numeric_vector(self.clusters,self.clusters_dico)
		if self.labels != []:
			self.labels_dico = self.find_labels_or_cluster_dico(self.labels)
			numeric_label		= self.generate_clus_or_label_numeric_vector(self.labels,self.labels_dico)
		self.data = []
		i = 0
		while i < len(self.table):
			point = self.table[i,:]
			self.data.append([])
			if numeric_label == []:self.data[i].append(-1)
			else: self.data[i].append(numeric_label[i])

			if numeric_cluster == []:self.data[i].append(-1)
			else: self.data[i].append(numeric_cluster[i])

			for j,element in enumerate(point):
				self.data[i].append([j,element])
			i+=1
		return self.data

	def str_labels_or_cluster_dico(self,dico=0):
		if dico == 0:
			dico = self.labels_dico
		txt = ""
		for i in dico:
			txt += str(i)+" "+str(dico[i])+"\n"
		return txt

	def data_to_table(self):#,label_index = "", cluster_index = ""):
		#ver si puedo darle la vuelta a la lista de clus y labels
		cluster_c = []
		label_c	 = []
		self.table = []
		for i,point in enumerate(self.data):
			j = 2
			self.table.append([])
			while j < len(point):
				self.table[i].append(point[j][1])
				j+=1
			label_c.append([point[0]])
			cluster_c.append([point[1]])
		self.table = np.array(self.table)
		label_c		= np.array(label_c)
		cluster_c	= np.array(cluster_c)
		self.table = np.hstack((label_c,cluster_c,self.table))
		return self.table.tolist()

	def eat_data(self,seen_data_size,p_food):
		nb_food_data = np.random.binomial(len(self.data) - seen_data_size, p_food)
		return self.eat_data_n_points(seen_data_size,nb_food_data)

	def eat_data_n_points(self,seen_data_size,nb_food_data):
		i = 0
		data_size = len(self.data)
		while i < nb_food_data:

			if seen_data_size < data_size :
				replacing = self.data.pop(random.randint(seen_data_size,data_size - 1))
				self.data[random.randint(0,seen_data_size - 1)] = replacing
				data_size -= 1
			i += 1
		return self.data

	def eat_data_one_point(self,seen_data_size):
		return self.eat_data_n_points(seen_data_size,1)

"""
a = data()
print a.open("dataa",sep=" ",header=1,c_cluster=5)
print a.standardize_table()
print a.multiply_table_by_precision_and_convert_to_int()
print a.table_to_data()
print a.data_to_table()
"""
