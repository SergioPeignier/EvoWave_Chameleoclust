#!/usr/bin/python
from definitions import *
import pandas as pd
import numpy as np
from sklearn.utils.linear_assignment_ import linear_assignment


def convert_contingence_table_arff_to_true(contingence_table, true_dataset):
	new_contingence_table = pd.DataFrame(np.zeros((len(true_dataset.true_arff_dico.keys()),len(contingence_table.columns))),index = true_dataset.true_arff_dico.keys(), columns = contingence_table.columns)
	for k in true_dataset.df_true_dico.keys():
		for t in true_dataset.df_true_dico[k]:
			if k in contingence_table.index:
				new_contingence_table.loc[t,:] += contingence_table.loc[k,:]
	#print true_dataset.df_true_dico
	#print new_contingence_table
	return new_contingence_table

def Subspaces_from_DF(DF):
	return DF.notnull()
	
def ecountdimpercluster(dfphenotype):
	subspaces = Subspaces_from_DF(dfphenotype)
	return subspaces.sum(1)
		
def egetavgdim(dfphenotype):
	dims = ecountdimpercluster(dfphenotype)
	return sum(dims)*1./len(dims)
	
def egetnbcorepoints(dfphenotype):
	return len(list(dfphenotype.index))

def ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity=0.0):
	percentage_points_per_cluster = found_clusters_effective * 1./sum(found_clusters_effective)
	return percentage_points_per_cluster >= threshold_cluster_validity

def Compute_only_entropy(contingence_table, valid_clusters):
	contingence_table = contingence_table.loc[:,valid_clusters]
	found_clusters_effective = contingence_table.sum(0)
	p_h_in_c = contingence_table*1./(found_clusters_effective+0.0000000001)
	log_p_h_in_c = np.log(p_h_in_c)
	pre_ec = -1. * p_h_in_c * log_p_h_in_c
	pre_ec = pre_ec.fillna(0)
	ec     = pre_ec.sum(0)
	num    = sum(ec*found_clusters_effective)
	denum  = sum(found_clusters_effective) * np.log(len(contingence_table.index))
	return 1. - num*1./denum
		
def Entropy(cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	valid_clusters = ValidityClusterChecking(contingence_table.sum(0),threshold_cluster_validity)
	return Compute_only_entropy(contingence_table, valid_clusters)

def Compute_only_accuracy(contingence_table, valid_clusters, found_clusters_effective):
	best_matching_hidden_cluster = contingence_table==contingence_table.max(0)
	best_matching_hidden_cluster_weight = 1./best_matching_hidden_cluster.sum(0)
	correctly_predicted_objects  = contingence_table * best_matching_hidden_cluster * best_matching_hidden_cluster_weight
	correctly_predicted_objects  *= valid_clusters
	accuracy = sum(correctly_predicted_objects.sum(0)) * 1./sum(found_clusters_effective )
	return accuracy

def Accuracy(cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	return Compute_only_accuracy(contingence_table, valid_clusters, found_clusters_effective)

def mapped(contingence_table):
	mapped_clusters = (contingence_table.T * 1./contingence_table.T.sum(0)).T
	return mapped_clusters == mapped_clusters.max(0)

def Compute_only_recall(contingence_table, found_clusters_effective, valid_clusters, mapped_clusters):
	num = mapped_clusters * contingence_table
	num = num.loc[:,valid_clusters]
	num = num.T.sum(0)
	denum = contingence_table.T.sum(0)
	ans = num * 1./(denum+0.0000000001)

def Recall(cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	mapped_clusters = mapped(contingence_table)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	return Compute_only_recall(contingence_table, found_clusters_effective, valid_clusters, mapped_clusters)
	
def Compute_only_precision(contingence_table, found_clusters_effective, valid_clusters, mapped_clusters):
	num = mapped_clusters * contingence_table
	num = num.loc[:,valid_clusters]
	num = num.T.sum(0)
	denum = mapped_clusters * contingence_table.sum(0) * valid_clusters
	denum = denum.T.sum(0)
	ans = num * 1./(denum+0.0000000001)
	return ans
	
def Precision(cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	mapped_clusters = mapped(contingence_table)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	return Compute_only_precision(contingence_table, found_clusters_effective, valid_clusters, mapped_clusters)

def Compute_only_F1(contingence_table, valid_clusters, mapped_clusters):
	num = mapped_clusters * contingence_table
	num = num.loc[:,valid_clusters]
	num = num.T.sum(0)
	#print num
	denum_recall = contingence_table.T.sum(0)
	#print denum_recall
	recall = num * 1./(denum_recall+0.0000000001)
	#print recall
	denum_precision =  mapped_clusters * contingence_table.sum(0) * valid_clusters
	denum_precision = denum_precision.T.sum(0)
	#print denum_precision
	precision = num * 1./(denum_precision+0.0000000001)
	#print precision
	denum = recall + precision
	num = 2 * recall * precision
	f1 = sum(num* 1.0/(denum+0.0000000001))* 1./(len(num)+0.0000000001)
	return f1
	
def F1(contingence_table = "", cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	mapped_clusters = mapped(contingence_table)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	return Compute_only_F1(contingence_table, valid_clusters, mapped_clusters)
	
def Compute_only_CE(contingence_table, found_clusters_effective, valid_clusters):
	valid_contingence_table =  contingence_table.loc[:,valid_clusters]
	best_hidden_found_couples = linear_assignment(-valid_contingence_table)
	return sum([valid_contingence_table.iloc[couple[0],couple[1]] for couple in best_hidden_found_couples])*1./sum(contingence_table.sum(0))


def CE(contingence_table = "", cluster_hidden = "", cluster_found = "", threshold_cluster_validity=0.0):
	if contingence_table == "":
		contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	return Compute_only_CE(contingence_table, found_clusters_effective, valid_clusters)
	
	
def Structural_Evaluation(phenotype, cluster_hidden = "", cluster_found = "", threshold_cluster_validity = 0.0, noise_hidden_cluster_id = [666]):
	not_noise      = [ ch not in noise_hidden_cluster_id for ch in cluster_hidden]
	cluster_hidden = cluster_hidden[not_noise]
	cluster_found  = cluster_found[not_noise]
	contingence_table = pd.crosstab(cluster_hidden,cluster_found)
	found_clusters_effective = contingence_table.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	for i in list(contingence_table.loc[:,valid_clusters].columns):
		if i not in phenotype.index:
			return [0,0,0,0]
	clusters = phenotype.loc[list(contingence_table.loc[:,valid_clusters].columns),:]
	avg_core_points_dim = egetavgdim(phenotype)
	nb_core_points      = egetnbcorepoints(phenotype)
	avg_clusters_dim = egetavgdim(clusters)
	nb_clusters      = egetnbcorepoints(clusters)	
	evaluation = [nb_core_points,avg_core_points_dim,nb_clusters,avg_clusters_dim]
	return evaluation 

def Pairwise_DF_Subspaces_Intersection(DF1,DF2):
	ans = pd.DataFrame([[sum(DF1.loc[i,:]&DF2.loc[j,:]) for j in DF2.index] for i in DF1.index],columns = list(DF2.index),index = list(DF1.index))
	return ans
	
def Pairwise_DF_Subspaces_Union(DF1,DF2):
	ans = pd.DataFrame([[sum(DF1.loc[i,:]|DF2.loc[j,:]) for j in DF2.index] for i in DF1.index],columns = list(DF2.index),index = list(DF1.index))
	return ans
	
def SubObjects_contingency_table(objects_contingence_table , subspaces_contigence_table):
	for i in list(objects_contingence_table.columns):
		if i not in subspaces_contigence_table.columns:
			return pd.DataFrame()
	subspaces_contigence_table = subspaces_contigence_table.loc[:,list(objects_contingence_table.columns)]
	return objects_contingence_table * subspaces_contigence_table
	
def Compute_only_RNIA(subobjects_intersection, subobjects_union,valid_clusters):
	subobjects_intersection = subobjects_intersection.loc[:,valid_clusters]
	subobjects_union = subobjects_union.loc[:,valid_clusters]
	I = sum(subobjects_intersection.sum(0))
	U = sum(subobjects_union.sum(0))
	return I *1./ (U+0.0000000001)
	
def Compute_only_SSCE(subobjects_intersection, subobjects_union,valid_clusters):
	subobjects_intersection = subobjects_intersection.loc[:,valid_clusters]
	#print subobjects_intersection
	subobjects_union = subobjects_union.loc[:,valid_clusters]
	#print subobjects_union
	best_hidden_found_couples = linear_assignment(-subobjects_intersection)
	return sum([subobjects_intersection.iloc[couple[0],couple[1]] for couple in best_hidden_found_couples])*1./(sum(subobjects_union.sum(0))+0.0000000001)

def Compute_coverage(found_clusters_effective, valid_clusters):
	return sum(found_clusters_effective[valid_clusters]) * 1.0/sum(found_clusters_effective)

def Functional_Evaluation(phenotype, cluster_hidden, cluster_found, true_dataset = [], threshold_cluster_validity = 0.0, noise_hidden_cluster_id = [666]):	
	not_noise      = [ ch not in noise_hidden_cluster_id for ch in cluster_hidden]
	cluster_hidden = cluster_hidden[not_noise]
	cluster_found  = cluster_found[not_noise]
	contingence_table_arff = pd.crosstab(cluster_hidden,cluster_found)
	found_clusters_effective = contingence_table_arff.sum(0)
	valid_clusters = ValidityClusterChecking(found_clusters_effective,threshold_cluster_validity)
	if not isinstance(true_dataset,list):
		contingence_table = convert_contingence_table_arff_to_true(contingence_table_arff, true_dataset)
		subspace_phenotype = Subspaces_from_DF(phenotype)
		I = Pairwise_DF_Subspaces_Intersection(true_dataset.clusters_subspaces,subspace_phenotype)	
		U = Pairwise_DF_Subspaces_Union(true_dataset.clusters_subspaces,subspace_phenotype)
		I = SubObjects_contingency_table(contingence_table ,I)
		U = SubObjects_contingency_table(contingence_table ,U)
		mapped_clusters = mapped(contingence_table)
		evaluation = [Compute_only_entropy(contingence_table, valid_clusters),
		Compute_only_accuracy(contingence_table, valid_clusters, found_clusters_effective),
		Compute_only_F1(contingence_table, valid_clusters, mapped_clusters),
		Compute_only_CE(contingence_table, found_clusters_effective, valid_clusters),
		Compute_only_RNIA(I,U,valid_clusters),
		Compute_only_SSCE(I,U,valid_clusters),
		Compute_coverage(found_clusters_effective, valid_clusters)]
	else:
		contingence_table = contingence_table_arff
		mapped_clusters = mapped(contingence_table)
		evaluation = [Compute_only_entropy(contingence_table, valid_clusters),
		Compute_only_accuracy(contingence_table, valid_clusters, found_clusters_effective),
		Compute_only_F1(contingence_table, valid_clusters, mapped_clusters),
		Compute_only_CE(contingence_table, found_clusters_effective, valid_clusters),
		None,
		None,
		Compute_coverage(found_clusters_effective, valid_clusters)]
	return evaluation
