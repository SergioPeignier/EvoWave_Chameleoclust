#!/usr/bin/python
import sys
import time
import os
import pandas as pd
from modevoevo       import *
from genotypemanager import *
from definitions	 import *
import random
import numpy
import math
import copy
#import operator


def autocorr(x, t=1):
	return(numpy.corrcoef(numpy.array([x[0:len(x)-t], x[t:len(x)]]))[1,0])

class EvoEvoManager:
	def __init__(self,eparams_input="", elogger = "", save_info_in_log_file = 1):
		self.eparams_evoevo = eparams_input
		self.modevoevo_simulation = esimulation(self.eparams_evoevo, elogger,save_info_in_log_file)
		self.eparams_random_walk = copy.deepcopy(eparams_input)
		self.eparams_random_walk.population_size = 1
		self.modevoevo_random_walk_simulation = esimulation(self.eparams_random_walk, elogger, save_info_in_log_file)
		self.parents_stats = []
		self.children_stats = []
		self.aggregatedstats = pd.DataFrame()
	
	def ecomputeNchildrenfromparentalgenome(self,parentalgenome,dataset="",replace=""):
		self.setdata(dataset,replace)
		self.modevoevo_simulation.esetpopulationgenome(parentalgenome)
		self.modevoevo_simulation.ecomputefitnesses()
		self.modevoevo_simulation.egetstats()
		self.parents_stats = copy.deepcopy(self.modevoevo_simulation.egetstatsasdf())
		self.modevoevo_simulation.eiterate(1)
		self.modevoevo_simulation.egetstats()
		self.children_stats = self.modevoevo_simulation.egetstatsasdf()

	def ecomputeNchildrenfromparentalfile(self,parentalgenomefile,dataset="",replace="",extension="",dico=""):
		parentalgenome = eloadgenomefromfile(parentalgenomefile,extension)
		ecomputeNchildrenfromparentalgenome(dataset,parentalgenome)
		
		
	def ecomputeFv(self):
		parentalfitMIN = min(self.parents_stats["fitness"])
		parentalfitMAX = max(self.parents_stats["fitness"])
		children_fitnesses = self.children_stats["fitness"]
		fv = sum((children_fitnesses <= parentalfitMAX) & (children_fitnesses >= parentalfitMIN))
		return fv * 1.0/len(list(children_fitnesses))

	def ecomputeFw(self):
		parentalfitMIN = min(self.parents_stats["fitness"])
		children_fitnesses = self.children_stats["fitness"]
		fv = sum(children_fitnesses < parentalfitMIN)
		return fv * 1.0/len(list(children_fitnesses))

	def ecomputeFb(self):
		parentalfitMAX = max(self.parents_stats["fitness"])
		children_fitnesses = self.children_stats["fitness"]
		fv = sum(children_fitnesses > parentalfitMAX)
		return fv * 1.0/len(list(children_fitnesses))

	def ecomputeFwFvFb(self):
		fv = self.ecomputeFv()
		fw = self.ecomputeFw()
		fb = self.ecomputeFb()
		return pd.DataFrame([[fw,fv,fb]], columns =  STD_FW_FV_FB_COLUMNS, index = [self.modevoevo_simulation.current_generation])

	def ecomputestatsspeedvectors(self):
		parental_stat = list(self.parents_stats.loc[0,:])
		return pd.DataFrame([[self.children_stats.loc[indx,col] - parental_stat[i] for i,col in enumerate(self.children_stats.columns)] for indx in list(self.children_stats.index)],columns=list(self.children_stats.columns),index = list(self.children_stats.index))

	def setdata(self,dataset,replace):
		if dataset != "" and replace != "":
			self.modevoevo_random_walk_simulation.esetdata(dataset,replace)
			self.modevoevo_simulation.esetdata(dataset,replace)
	
	def emakerandomwalkandcomputerhotau(self,parentalgenome_i,funtion_to_save,columns_names,lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 0,index = 0,dataset = "",replace = ""):
		signal = []
		self.setdata(dataset,replace)
		parentalgenome   = parentalgenome_i
		self.modevoevo_random_walk_simulation.esetpopulationgenome(parentalgenome)
		self.modevoevo_random_walk_simulation.ecomputefitnesses()
		parentalstats	 = self.modevoevo_random_walk_simulation.egetstats()[0][:]
		t = 0
		while t < lengthwalk:
		 	trials = 0
	 		childstats = 0
	 		childgenome = []
			while (not childstats and trials < max_trials):
	 		 	self.modevoevo_random_walk_simulation.eiterate(1)
				childstats	 = random.sample(self.modevoevo_random_walk_simulation.egetstats(),1)[0][:]
				restart = (avoidnulllengthindividuals and childstats[STATS_GENOME_LENGTH] == 0)	or childstats[STATS_GENOME_LENGTH] >= self.modevoevo_random_walk_simulation.eparams.genom_limit_size
				if restart:
	 	 			childstats = 0
	 	 			simulation.esetpopulationgenome(parentalgenome)
	 			trials += 1
			#print(new_parental_organism[0])
			childgenome = self.modevoevo_random_walk_simulation.egetindividualgenome(int(childstats[STATS_INDEX]))[:]
			signal.append(funtion_to_save([parentalstats,parentalgenome,childstats,childgenome]))
			parentalstats = childstats[:]
			parentalgenome = childgenome[:]
			self.modevoevo_random_walk_simulation.esetpopulationgenome(parentalgenome)
			t += 1
		#print(signal)
		rho = autocorr(signal,t=1)
		tau = -1.0/(math.log(abs(rho)+0.000000001)+0.000000001)
		return (signal,pd.DataFrame([[rho,tau]],columns = columns_names, index = [self.modevoevo_random_walk_simulation.current_generation]))

	def emakerandomwalkandcomputerhotaufitness(self,parentalgenome,lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 0,index = 0,dataset = "",replace = ""):
		fitness_function = lambda prtst_prtg_chst_chg : prtst_prtg_chst_chg[2][STATS_FITNESS]
		return self.emakerandomwalkandcomputerhotau(parentalgenome,fitness_function, STD_FITNESS_RANDOM_WALK_STUDY,lengthwalk,max_trials,avoidnulllengthindividuals ,index, dataset, replace)
	
	def emakerandomwalkandcomputerhotaufv(self,parentalgenome,lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 0,index = 0,dataset = "",replace = ""):
		def fv_function(prtst_prtg_chst_chg):
			self.ecomputeNchildrenfromparentalgenome(prtst_prtg_chst_chg[3])
			return self.ecomputeFv()
		return self.emakerandomwalkandcomputerhotau(parentalgenome,fv_function, STD_Fv_RANDOM_WALK_STUDY, lengthwalk,max_trials,avoidnulllengthindividuals ,index, dataset, replace)

	def emakerandomwalkandcomputerhotaufw(self,parentalgenome,lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 0,index = 0,dataset = "",replace = ""):
		def fw_function(prtst_prtg_chst_chg):
			self.ecomputeNchildrenfromparentalgenome(prtst_prtg_chst_chg[3])
			return self.ecomputeFw()
		return self.emakerandomwalkandcomputerhotau(parentalgenome,fw_function, STD_Fw_RANDOM_WALK_STUDY, lengthwalk,max_trials,avoidnulllengthindividuals ,index, dataset, replace)

	def emakerandomwalkandcomputerhotaufb(self,parentalgenome,lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 0,index = 0,dataset = "",replace = ""):
		def fb_function(prtst_prtg_chst_chg):
			self.ecomputeNchildrenfromparentalgenome(prtst_prtg_chst_chg[3])
			return self.ecomputeFb()
		return self.emakerandomwalkandcomputerhotau(parentalgenome,fb_function,STD_Fb_RANDOM_WALK_STUDY, lengthwalk,max_trials,avoidnulllengthindividuals ,index , dataset, replace)

	def eaddtoaggregatedstats(self,generation = 0):
		self.modevoevo_random_walk_simulation.current_generation = generation
		self.modevoevo_simulation.current_generation = generation
		#agregate fvfwfb with rho tau afterwards  
		self.aggregatedstats = pd.concat([self.aggregatedstats,self.ecomputeFwFvFb()])
		

"""
self.save_info_in_log_file = save_info_in_log_file
if nb_individuals_for_evoevo > 0:
			self.eparams_evoevo = copy.deepcopy(eparams_input)
			self.eparams_evoevo.population_size = nb_individuals_for_evoevo
			self.evoevomanager = EvoEvoManager(self.eparams_evoevo, self.elogger , save_info_in_log_file)	

"""




"""
def ecomputeNchildrenfromparent(self,parameters,log,dataset,name_save_log):
		children_simu = esimulation(parameters, self.elogger)
		children_simu.esetdatafromdf(dataset,replace=1)
		children_simu.eiterate(1)
		children_simu.ecomputefitnesses()
		return p_children_simu_c

codear funccion que hace lo mismo pero que calcula la evolvabilidad
def emakerandomwalkandsavefitnesssignal(a,parentalgenomefile,lengthwalk=100,max_trials=100):
	signal = []
	parentalgenome		= a.eloadgenomefromfile(parentalgenomefile)
	t = 0
	new_parental_organism = 0
	while t < lengthwalk:
		a.esetpopulationgenome(parentalgenome)
		a.eiterate(1)
		parentalstats	 = a.egetstats()
		trials = 0
		while (not new_parental_organism and trials < max_trials):
	 	 new_parental_organism = random.sample(parentalstats, 1)[0]
	 	 if new_parental_organism[STATS_GENOME_LENGTH] == 0:
	 	 	new_parental_organism = 0
	 	 trials += 1
		signal.append(new_parental_organism[STATS_FITNESS])
		parentalgenome = a.egetindividualgenome(new_parental_organism[STATS_INDEX])
		t += 1
	return signal
"""

"""
def ecomputerhotaumanyindividualssameparams(parameters,dataset,stream_init_size,parentalgenomefiles,title,datasetName,parametersName,generationpatt = "generation=",sep="_",indivpatt="indiv=",computeeach=100,savename="RhoTauFile.txt",lengthwalk=100,max_trials=100,avoidnulllengthindividuals = 1):
	dataset.table_to_data()
	random.shuffle(dataset.data)
	log = logger(parent_folder_name = "./RESULTS_EVO_EVO_MESURES",folder_name="/"+title+"_DATASET_"+datasetName+"_PARAMETERS_"+parametersName)
	a = esimulation(parameters,log)
	a.esetdata(dataset.data[0:stream_init_size],1,dataset.str_labels_or_cluster_dico())
	ans = [["seed","generation","index","fitness","coding_ratio","length","rho","tau"]]
	t = len(parentalgenomefiles)
	for i,pg in enumerate(parentalgenomefiles):
		print t - i
		seed	= pg.split("seed_")[1]
		seed	= seed.split("_")[0]
		indiv_description = os.path.basename(pg).split(".")[0]
		generation				= indiv_description.split(generationpatt)[1]
		generation				= int(generation.split(sep)[0])
		indiv						 = indiv_description.split(indivpatt)[1]
		indiv						 = int(indiv.split(sep)[0])
		if not (generation % computeeach)-1:
			ans.append([int(seed),int(generation),int(indiv)]+emakerandomwalkandcomputerhotau(a,pg,lengthwalk,max_trials,avoidnulllengthindividuals))
	a.elogger.table_to_file(ans, savename,rewrite = True)
	return ans


def ecomputeFvmanyindividualssameparams(parameters,dataset,stream_init_size,parentalgenomefiles,title,datasetName,parametersName,generationpatt = "generation=",sep="_",indivpatt="indiv=",computeFveach=100,savename="FvFile.txt"):
	dataset.table_to_data()
	random.shuffle(dataset.data)
	log = logger(parent_folder_name = "./RESULTS_EVO_EVO_MESURES",folder_name="/"+title+"_DATASET_"+datasetName+"_PARAMETERS_"+parametersName)
	a = esimulation(parameters,log)
	a.esetdata(dataset.data[0:stream_init_size],1,dataset.str_labels_or_cluster_dico())
	ans = [["seed","generation","index","fitness","coding_ratio","length","Fv"]]
	t = len(parentalgenomefiles)
	for i,pg in enumerate(parentalgenomefiles):
		print t - i
		seed	= pg.split("seed_")[1]
		seed	= seed.split("_")[0]
		indiv_description = os.path.basename(pg).split(".")[0]
		generation				= indiv_description.split(generationpatt)[1]
		generation				= int(generation.split(sep)[0])
		indiv						 = indiv_description.split(indivpatt)[1]
		indiv						 = int(indiv.split(sep)[0])
		if not (generation % computeFveach)-1:
			ans.append([int(seed),int(generation),int(indiv)]+emakereplicationsandcomputeFv(a,pg))
	a.elogger.table_to_file(ans, savename,rewrite = True)
	return ans
def eclassifydatafromsavedindividual(parameter,dataset,parentalgenomefile,genotypesavename,phenotypesavename,classifysavename,logname="logger_classify_with_given_individual.txt",time_name="time_for_classification_given_individual"):
	dataset.table_to_data()
	random.shuffle(dataset.data)
	log = logger(parent_folder_name = "./",folder_name=os.path.dirname(parentalgenomefile),file_name=logname,time_name=time_name)
	a = esimulation(parameter,log)
	parentalgenome = a.eloadgenomefromfile(parentalgenomefile)
	a.esetpopulationgenome(egenome=parentalgenome)
	a.esetdata(dataset.data,1,dataset.str_labels_or_cluster_dico())
	a.ecomputefitnesses()
	indiv = a.egetbestindividualindex(a.egetstats())
	a.esavegenotype(indiv,genotypesavename)
	a.esavephenotype(indiv,phenotypesavename)
	a.esaveclassifieddata(indiv,classifysavename)

def eonecompleteiteration(a,i,nb_best_individuals_to_save,saveinfoeach,savegenotypes,savephenotypes,savetables,shufflegenome,pg_convertion):
	a.eiterate(1)
	a.esavestats("stats.dat")
	pop_stats		= a.egetstats()
	N_best_indiv = a.egetNbestindividualsindex(pop_stats,nb_best_individuals_to_save)
	a.emakesaves(i,saveinfoeach,N_best_indiv,savegenotypes,savephenotypes,savetables)
	if shufflegenome:
		a.eshuffleallpopulationgenome()
	if pg_convertion:
		a.emodifypseudogenesallpopulation(pg_convertion)
	a.elogger.write_time_mesures(i)


def erunexperiments(title,parameters,datasets,nb_best_individuals_to_save,savegenotypes,savephenotypes,savetables,pg_convertion,shufflegenome,saveinfoeach,proportion_2_update,nb_training_generations_per_update, nb_total_generations,stream_size):
	ans = []
	numpy.random.seed(0)
	random.seed(0) 
	for k in datasets:
		#numpy.random.shuffle(datasets[k].table)
		
		#random.shuffle(datasets[k].data)
		for j,p in enumerate(parameters):
			#initialization

			for exp in parameters[p]:
				random.seed(exp.prng_seed)
				numpy.random.seed(exp.prng_seed)
				datasets[k].table_to_data()
				random.shuffle(datasets[k].data)
				log = logger(parent_folder_name = "./RESULTS_"+title+"_"+str(k)+"_"+str(p),seed = exp.prng_seed)
				a = esimulation(exp,log)
				i = 0
				a.esetdata(datasets[k].data[0:stream_size],1,datasets[k].str_labels_or_cluster_dico())
				a.ecomputefitnesses()
				a.esavestats("stats.dat")
				if shufflegenome:
					a.eshuffleallpopulationgenome()
				if pg_convertion:
					a.emodifypseudogenesallpopulation(pg_convertion)
				nb_updates = 0
	 	 		while 1:
	 	 			j = 0
	 	 			while j < nb_training_generations_per_update and i < nb_total_generations:
						eonecompleteiteration(a,i,nb_best_individuals_to_save,saveinfoeach,savegenotypes,savephenotypes,savetables,shufflegenome,pg_convertion)
	 	 	 	 		j += 1
	 	 	 	 		i += 1
	 	 	 	 	if i == nb_total_generations:
	 	 	 	 		break		
	 	 	 	 	pos1 = (len(datasets[k].table)+((nb_updates*(int(proportion_2_update*stream_size))+stream_size)%len(datasets[k].table)))%len(datasets[k].table)
	 	 	 	 	pos2 = (len(datasets[k].table)+(((nb_updates+1)*(int(proportion_2_update*stream_size))+stream_size)%len(datasets[k].table)))%len(datasets[k].table)
	 	 	 	 	if pos1 < pos2:
	 	 	 	 		a.esetdata(datasets[k].data[pos1:pos2],0)
	 	 	 	 		#print str(pos1)+"  "+str(pos2)
	 	 	 	 	if pos2 < pos1:
	 	 	 	 		a.esetdata(datasets[k].data[pos1:len(datasets[k].data)]+datasets[k].data[0:pos2],0)
	 	 	 	 		#print str(pos1)+"  "+str(len(datasets[k].data))+" 0 "+str(pos2)
	 	 	 	 	nb_updates += 1
	 	 	 	 		
				pop_stats		= a.egetstats()
				N_best_indiv = a.egetNbestindividualsindex(pop_stats,nb_best_individuals_to_save)
				if i > 1:
					a.emakesaves(i,i-1,N_best_indiv,savegenotypes,savephenotypes,savetables)
				ans.append(a)
	return ans
"""