#!/usr/bin/python
# -*- coding: utf-8 -*-
from definitions import *
import random

class eboudaryCondition:

	def __init__(self,
				 emin,
				 emax):
		''' 
			Creats an eboudaryCondition instance
			:param emin: Lower bound
			:type  emin: int
			:param emax: Upper bound
			:type  emax: int
			'''
		self.emin = emin
		self.emax = emax
		
	def get_this_as_tuple(self):
		''' 
			Returns eboudaryCondition as a tuple
			'''
		return (self.emin,self.emax)

	def get_this_from_tuple(self,tuple_in):
		''' 
			Makes an eboudaryCondition from a tuple
			'''
		self.emin = tuple_in[0]
		self.emax = tuple_in[1]


class estatisticalDistributionLaw:
	def __init__(self,
				 emin = 0,
				 emax = 0,
				 emean = 0,
				 estandard_deviation = 0,
				 law = "",
				 transition_matrix = ""):
		'''
			Creates an estatisticalDistributionLaw instance
			
			:param emin: Lower bound of uniform distribution
			:type  emin: int
			:param emin: Upper bound of uniform distribution
			:type  emin: int
			:param emean: Gaussian distribution mean
			:type  emean: float
			:param estandard_deviation: Standard deviation for gaussian distribution
			:type  estandard_deviation: float
			:param law: ID of statistical distribution
			:type  law: int
			:param transition_matrix: transition matrix probabilities
			:type  law: list
			'''
		self.law										 = law
		if law == GAUSSIAN:
			self.emin								= emin
			self.emax								= emax
			self.emean							 = emean
			self.estandard_deviation = estandard_deviation
			self.transition_matrix	 = []

		if law == UNIFORM:
			self.emin								= emin
			self.emax								= emax
			self.emean							 = emean
			self.estandard_deviation = estandard_deviation
			self.transition_matrix	 = []

		if law == TRANSITION_MATRIX:
			self.emin								= 0
			self.emax								= 0
			self.emean							 = 0
			self.estandard_deviation = 0
			self.transition_matrix	 = self.make_wheel_of_fortune(self.transitionProbability2stepProbability(transition_matrix))

		if law != GAUSSIAN and law != UNIFORM and law != TRANSITION_MATRIX:
			print "This law has not been implemented yet [0,1024[ uniform law will be used"
			self.emin								= 0
			self.emax								= 1025
			self.emean							 = 0
			self.estandard_deviation = 0
			self.transition_matrix	 = []
			
	def make_wheel_of_fortune(self,matrix):
		'''
			Makes a wheel of fortune from a transition matrix
			
			:param matrix: transition matrix probabilities
			:type  law: list
			'''
		new_matrix = [[0]+[sum(l[:i+1])*1./sum(l) for i in range(len(l))] for l in matrix]
		return new_matrix

	def get_this_as_tuple(self):
		''' 
			Return the object in a long tuple
			'''
		return (self.emin,self.emax,self.emean,self.estandard_deviation,self.law,self.transition_matrix)

	def get_this_from_tuple(self,tuple_to_load):
		''' 
			Load a eparameters objet from a set of parameters entered in a tuple
			:param tuple_to_load : Tuple of parameters):
			:type  tuple_to_load : tuple
			'''
		self.emin								    = tuple_to_load[0]
		self.emax								    = tuple_to_load[1]
		self.emean								    = tuple_to_load[2]
		self.estandard_deviation	                = tuple_to_load[3]
		self.law									= tuple_to_load[4]
		self.transition_matrix		                = self.stepProbability2transitionProbability(tuple_to_load[5])

	def circularity(self,x,len):
		return (len+(x%len))%len

	def circularitydrawrandomnumber(self,x):
		if self.emax-self.emin == 0 : return	self.emax
		else: return ((self.emax-self.emin)+(x%(self.emax-self.emin)))%(self.emax-self.emin)+self.emin

	def transitionProbability2stepProbability(self,transition_matrix):
		step_transition_matrix = [[0 for c in l] for l in enumerate(transition_matrix)]
		for i,l in enumerate(transition_matrix):
			for j,c in enumerate(transition_matrix[i]):
				step_transition_matrix[i][self.circularity(j-i,len(transition_matrix))] =	transition_matrix[i][j]
		return step_transition_matrix

	def stepProbability2transitionProbability(self,step_transition_matrix):
		transition_matrix = []
		for i,l in enumerate(step_transition_matrix):
			transition_matrix.append([])
			for j,c in enumerate(step_transition_matrix[i]):
				transition_matrix[i].append(0)
			for j,c in enumerate(step_transition_matrix[i]):
				transition_matrix[i][self.circularity(j+i,len(step_transition_matrix))] =	transition_matrix[i][j]
		return transition_matrix


	def drawrandomnumber(self):
		if self.law == GAUSSIAN:
			return(int(self.circularitydrawrandomnumber(int(random.gauss(self.emean,self.estandard_deviation)))))

		if self.law == UNIFORM:
			return(int(self.circularitydrawrandomnumber(random.uniform(self.emin,self.emax))))
		return 0
		

class egeneElementMutProbLaw:
	def __init__(self,mut_law_gene_elmnt = STD_GENE_MUTATION_PROBABILITIES):
		self.mut_law_gene_elmnt          = mut_law_gene_elmnt
		self.mut_law_wheel_of_fortune    = [0]+[sum(self.mut_law_gene_elmnt[:i+1])*1./sum(self.mut_law_gene_elmnt) for i in range(len(self.mut_law_gene_elmnt))]
	def get_this_as_tuple_of_arrays(self):
		return (self.mut_law_gene_elmnt,self.mut_law_wheel_of_fortune)

	def get_this_from_tuple_of_arrays(self,tuple_of_arrays):
		self.mut_law_gene_elmnt			 = tuple_of_arrays[0]
		self.mut_law_wheel_of_fortune    = tuple_of_arrays[1]

class eintergeniccut:
	def __init__(self,cut = 0,cut_order=STD_ORDER_FOR_INTERGENIC_CUT_OFF,emin=1,emax=STD_GENE_SIZE):
		self.cut = cut
		self.cut_order = cut_order
		self.emin = emin
		self.emax=emax
	def get_this_as_tuple(self):
		return tuple([self.cut,self.emin,self.emax]+self.cut_order)
	def get_this_from_tuple(self,tup):
		lis = list(tuple)
		self.cut = lis[0]
		self.emin = lis[1]
		self.emax = lis[2]
		self.cut_order = lis[3:len(lis)]


def generateeasymutationlawarray(dataset,nb_max_clusters, transition_matrix_c_nc = NONABSORBINGSTATE , init_c_mean = 0, init_c_std = 1, init_c_max=2, init_c_min = 0, cnc_law =  TRANSITION_MATRIX):
	position_lims = dataset.df[dataset.features].abs().max().max()
	cluster_max   = nb_max_clusters
	dim_max       = len(dataset.features)
	
	if cnc_law ==  TRANSITION_MATRIX:
		ans = [estatisticalDistributionLaw(law = TRANSITION_MATRIX, transition_matrix = transition_matrix_c_nc),
		estatisticalDistributionLaw(emin = -cluster_max/2 - 1,emax = cluster_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -dim_max/2 - 1,emax =  dim_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -position_lims-1,emax = position_lims+1,law = UNIFORM)]

	if cnc_law == GAUSSIAN:
		ans = [estatisticalDistributionLaw(emean = init_c_mean,estandard_deviation = init_c_std, law = GAUSSIAN,emax = init_c_max, emin = init_c_min ),
		estatisticalDistributionLaw(emin = -1,emax = cluster_max ,law = UNIFORM),
		estatisticalDistributionLaw(emin = -1,emax = dim_max ,law = UNIFORM),
		estatisticalDistributionLaw(emin = -position_lims-1,emax = position_lims+1,law = UNIFORM)]
	return ans
	