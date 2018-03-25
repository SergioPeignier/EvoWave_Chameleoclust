#!/usr/bin/python
import modevoevoC.modevoevo_c as modevoevo_c
import pandas as pd
from modevoevoPY.eparameters     import *
from modevoevoPY.logger          import *
from modevoevoPY.definitions     import *
from modevoevoPY.evaluationfunctions import *
import random
import numpy
import math
import time
import os

class esimulation:

	def __init__(self,eparams_input="", elogger = "", save_info_in_log_file = 0):
		self.save_info_in_log_file = save_info_in_log_file
		if elogger != "": self.elogger = elogger
		else: self.elogger             = logger()
		if eparams_input != "" : self.eparams = eparams_input
		else: self.eparams = eparams()
		self.ewritelogfile()
		self.p_simulation_c		= modevoevo_c.esimulation( *self.eparams.getastuple())
		self.nb_indiv_features  = modevoevo_c.egetindividualfeaturesnb(self.p_simulation_c)
		self.current_generation = 0
		self.stats				= [[0.0 for k in range(self.nb_indiv_features)] for l in range(self.eparams.population_size)]
		self.genome_getter		= [[0.0 for i in range(self.eparams.gene_size)] for j in range(self.eparams.genom_limit_size)]
		if self.eparams.fitness_mode == FITNESS_AEVOL:
			self.phenotype_getter	= [[1,1]+[[0 for k in range(2)] for i in range(self.eparams.boundary_conditions[2].emax-self.eparams.boundary_conditions[2].emin -1)] for j in range(self.eparams.boundary_conditions[2].emax-self.eparams.boundary_conditions[1].emin)]
			self.phenotype_dim_clus = [0 for j in range(self.eparams.boundary_conditions[2].emax-self.eparams.boundary_conditions[2].emin)]
		else:
			self.phenotype_getter	= [[1,1]+[[0 for k in range(2)] for i in range(self.eparams.boundary_conditions[2].emax-self.eparams.boundary_conditions[2].emin -1)] for j in range(self.eparams.boundary_conditions[1].emax-self.eparams.boundary_conditions[1].emin)]
			self.phenotype_dim_clus = [0 for j in range(self.eparams.boundary_conditions[1].emax-self.eparams.boundary_conditions[1].emin)]
		self.data = [[1,1]+[[0 for k in range(2)] for i in range(self.eparams.boundary_conditions[2].emax-self.eparams.boundary_conditions[2].emin -1)] for j in range(self.eparams.size_data_buffer)]
		self.aggregatedstats = {}

	def esetparameters(self,parameters):
		self.eparams = parameters
		modevoevo_c.esetparameters(self.p_simulation_c,*self.eparams.getastuple())
		
	def ewritelogfile(self,name=""):
		if self.save_info_in_log_file:
			if name == "":self.elogger.save_object_to_file(self.eparams,self.eparams.log_file)
			else         :self.elogger.save_object_to_file(self.eparams,name)
			if self.save_info_in_log_file:
				self.elogger.write_time_mesures(0)

	def eloadparametersfromsimulation(self):
		self.elogger.eloadfromtuple(modevoevo_c.egetparameters(self.p_simulation_c))

	def eiterate(self,eiter=1):
		if self.p_simulation_c == 0: return 0
		#X = dir(self)
		modevoevo_c.eiterate(self.p_simulation_c,eiter)
		#Y = dir(self)
		#for x in X:
		#	if x not in Y:
		#		print x 
		self.current_generation += eiter
		if self.save_info_in_log_file:
			self.elogger.write_time_mesures(self.current_generation)
		#self.egetstats()
		#self.egetstatsasdf()
		
	def emodifypseudogenesoneindividual(self,modification_laws,individual):
		if self.p_simulation_c == 0: return 0
		modification_laws_tuple = [modification_law.get_this_as_tuple() for modification_law in modification_laws]
		modevoevo_c.emodifypseudogenesoneindividual(self.p_simulation_c,modification_laws_tuple,individual)

	def eordergenomeindividual(self,individual):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eordergenomeindividual(self.p_simulation_c,individual)

	def eordergenomeallpopulation(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eordergenomeallpopulation(self.p_simulation_c)

	def eshuffleallpopulationgenome(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eshuffleallpopulationgenome(self.p_simulation_c)

	def eshuffleoneindividualgenome(self,individual):
		modevoevo_c.eshuffleoneindividualgenome(self.p_simulation_c,individual)

	def ecomputefitnesses(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.ecomputefitnesses(self.p_simulation_c)
		#self.egetstats()
		#self.egetstatsasdf()
	
	def esetsizedatabuffer(self,size_data_buffer = -1,size_data_point = -1):		
		if size_data_buffer == -1:
			size_data_buffer = self.eparams.size_data_buffer
		if size_data_point == -1:
			size_data_point = self.eparams.size_data_point
		self.eparams.size_data_point   = int(size_data_point)
		self.eparams.size_data_buffer  = int(size_data_buffer)
		return modevoevo_c.esetsizedatabuffer(self.eparams.size_data_buffer,self.eparams.size_data_point,self.p_simulation_c)

	def eresetdata(self,newsize):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eresetdata(self.p_simulation_c,newsize)
		
	def esetdata(self,data,replace=0,dico="",file_name="TRANSLATION_DICO.PICKLE"):
		if self.p_simulation_c == 0: return 0
		if dico != "": self.elogger.save_object_to_file(dico,file_name)
		modevoevo_c.esetdata(self.p_simulation_c,data,replace)	

	def esetdatafromdf(self,datadf,replace=0,dico="",file_name="TRANSLATION_DICO.PICKLE"):
		if self.p_simulation_c == 0: return 0
		if dico != "": self.elogger.save_object_to_file(dico,file_name)
		modevoevo_c.esetdata(self.p_simulation_c,datadf.data[0:self.eparams.size_data_buffer],replace)	
	
	def esaveobject(self,object,file_name):
		self.elogger.save_object_to_file(object,file_name)
	
	def egetindividualgenome(self,eindiv=0):
		if self.p_simulation_c == 0: return 0
		genome_size = modevoevo_c.egetgenomesize(self.p_simulation_c,eindiv)
		if genome_size > self.eparams.gene_size: self.genome_getter = [[0.0 for i in range(self.eparams.gene_size)] for j in range(genome_size)]
		modevoevo_c.egetindividualgenome(self.p_simulation_c,eindiv,self.genome_getter)
		return self.genome_getter[0:genome_size]
	
	def esetindividualgenome(self,eindiv=0,egenome=[]):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esetindividualgenome(self.p_simulation_c,eindiv,egenome)

	def esetbestindividualgenome(self,egenome=[]):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esetbestindividualgenome(self.p_simulation_c,egenome)

	def egetbestindividualgenome(self):
		return self.egetindividualgenome(self.egetbestindividualindex()[0])
		
	def esetpopulationgenome(self,egenome=[]):
		i = 0
		while i < self.eparams.population_size:
			self.esetindividualgenome(i,egenome)
			i += 1
		self.esetbestindividualgenome(egenome)
			
	def echangeseed(self,newSeed):
		if self.p_simulation_c == 0: return 0
		self.eparams.prng_seed = newSeed
		modevoevo_c.echangeseed(self.p_simulation_c,newSeed)
		
	def eprintgenomes(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eprintgenomes(self.p_simulation_c)
		
	def emodifypseudogenesallpopulation(self,modification_laws):
		if self.p_simulation_c == 0: return 0
		if len(modification_laws) > self.eparams.gene_size:
			modification_laws[:] =	modification_laws[0:self.eparams.gene_size]
		if len(modification_laws) < self.eparams.gene_size:
			raise "WARNING: Not enough information to replace pseudo-genes by random genes"
		modification_laws_tuple = [modification_law.get_this_as_tuple() for modification_law in modification_laws]
		modevoevo_c.emodifypseudogenesallpopulation(self.p_simulation_c,modification_laws_tuple)
		
	def egetbestindividualindex(self,stats = ""):
		if isinstance(stats,str):
			stats = self.stats
		stats.sort(reverse = True)
		return [int(i[STATS_INDEX]) for i in stats[0:1]]

	def egetNbestindividualsindex(self,N=0):
		stats.sort(reverse = True)
		return [int(i[STATS_INDEX]) for i in stats[0:N]]

	def egetindividual(self,individual_index):
		if self.p_simulation_c == 0: return 0
		return modevoevo_c.egetindividual(self.p_simulation_c,individual_index)

	def esetindividual(self,capsule_individual,individual_index):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esetindividual(self.p_simulation_c,capsule_individual,individual_index)

	def emakeindividualclassifydata(self,eindiv=0):
		if self.p_simulation_c == 0: return 0
		return modevoevo_c.eclassifydata(self.p_simulation_c,eindiv,self.data)
	
	def egetindividualclassifieddata(self,eindiv=0):
		self.emakeindividualclassifydata(eindiv)
		return self.egetdata()

	def egetindividualclassifieddataasdf(self,eindiv=0):
		self.emakeindividualclassifydata(eindiv)
		return self.egetdataasdf()	
	
	
	def egetindividualphenotype(self,eindiv=0):
		if self.p_simulation_c == 0: return 0
		nb_c = modevoevo_c.egetindividualphenotype(self.p_simulation_c,eindiv,self.phenotype_getter,self.phenotype_dim_clus)
		return [self.phenotype_getter[i][0:d] for i,d in enumerate(self.phenotype_dim_clus[0:nb_c])]
	
	def datatopandasdf(self,data,dim_max = "",cluster_columns = ["cluster_hidden","cluster_found"]):
		if isinstance(dim_max,str):
			dim_max = self.eparams.boundary_conditions[DIMENSION_GENE_TYPE].emax
		df = [[None for i in xrange(dim_max)] for point in data]
		#print data
		for i,p in enumerate(data):
			j = 2
			while j < len(p):
				df[i][int(p[j][0])] = p[j][1]
				j += 1
			df[i].append(p[0])
			df[i].append(p[1])
		df = pd.DataFrame(df,columns = [str(d) for d in range(dim_max)]+cluster_columns)
		return df
	
	def egetindividualphenotypeasdataframe(self,eindiv=0):
		if self.p_simulation_c == 0: return 0
		nb_c = modevoevo_c.egetindividualphenotype(self.p_simulation_c,eindiv,self.phenotype_getter,self.phenotype_dim_clus)
		ph = [self.phenotype_getter[i][0:d] for i,d in enumerate(self.phenotype_dim_clus[0:nb_c])]
		if ph != []:
			ph_df = self.datatopandasdf(ph,cluster_columns = ["cluster_id","_"])
			ph_df = ph_df.set_index("cluster_id",drop = 1)
			ph_df = ph_df.drop("_",1)
		else:
			ph_df = pd.DataFrame()
		return ph_df
		
	def egetstats(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.egetstats(self.p_simulation_c,self.stats)
		return self.stats
		
	def egetstatsasdf(self,stats=""):
		if isinstance(stats,str):
			stats = self.stats
		self.statsdf = pd.DataFrame(stats,columns=STATS_FEATURES)
		self.statsdf = self.statsdf.set_index("index",drop = 1)
		self.statsdf = self.statsdf.drop(STATS_TO_DROP,1)
		return self.statsdf
		
	def eaggregatebestindividualevaluation(self,true_dataframe = [],g = 1,f = 1,s = 1,threshold_cluster_validity= 0.0, noise_hidden_cluster_id = [666]):
		best_individual_index =  self.egetbestindividualindex()[0]
		best_individual_evaluation = self.egetindividualevaluation(index = best_individual_index,
																   true_dataframe = true_dataframe,
																   stats = self.statsdf,
																   g = 1,
																   f = 1,
																   s = 1,
																   threshold_cluster_validity= threshold_cluster_validity,
																   noise_hidden_cluster_id = noise_hidden_cluster_id)		
		if "best" not in self.aggregatedstats.keys():
			self.aggregatedstats["best"] = pd.DataFrame()
		self.aggregatedstats["best"] = pd.concat([self.aggregatedstats["best"],best_individual_evaluation])
		
	def eaddtoaggregatedstats(self,quantile=[1],bounds = [-0.1],col="fitness"):
		stats_quantile = self.statsdf[col].quantile(quantile)
		next_quantile  = [ bounds[i] + q for i,q in enumerate(quantile) ]
		stats_quantile_bound = self.statsdf[col].quantile(next_quantile)
		for i,quant in enumerate(quantile):
			aggregation = self.eaggregatestats(self.statsdf,stats_quantile[quant],stats_quantile_bound[next_quantile[i]],col)
			if col not in self.aggregatedstats.keys():
				self.aggregatedstats[col] = {}
			if quant not in self.aggregatedstats[col].keys():
				self.aggregatedstats[col][quant] = pd.DataFrame()
			self.aggregatedstats[col][quant] = pd.concat([self.aggregatedstats[col][quant],aggregation])
		return self.aggregatedstats
	
	def egetindividualevaluation(self,index,true_dataframe, stats = "",g = 1,f = 1,s = 1,threshold_cluster_validity= 0.0, noise_hidden_cluster_id = [666]):
		ev_g = pd.DataFrame()
		ev_s = pd.DataFrame()
		ev_f = pd.DataFrame()
		if isinstance(stats,str):
			stats = self.stats
		if g :
			ev_g = pd.DataFrame([list(stats.loc[index,STD_USEFULL_STATS_AGGREGATION])],columns = STD_USEFULL_STATS_AGGREGATION,index= [self.current_generation])
		best_clust_mod = self.egetindividualclassifieddataasdf(eindiv =index)
		best_phenotype = self.egetindividualphenotypeasdataframe(eindiv = index)
		#print best_clust_mod
		if f:
			ev_f = Functional_Evaluation(best_phenotype,cluster_hidden =best_clust_mod["cluster_hidden"],cluster_found = best_clust_mod["cluster_found"], true_dataset = true_dataframe,threshold_cluster_validity = threshold_cluster_validity, noise_hidden_cluster_id = noise_hidden_cluster_id)
			ev_f = pd.DataFrame([ev_f],columns = STD_FUNCTIONAL_EVALUATION_COLUMNS,index = [self.current_generation])
		if s: 
			ev_s = Structural_Evaluation(phenotype = best_phenotype, cluster_hidden =best_clust_mod["cluster_hidden"],cluster_found = best_clust_mod["cluster_found"], threshold_cluster_validity = threshold_cluster_validity,noise_hidden_cluster_id = noise_hidden_cluster_id)
			ev_s = pd.DataFrame([ev_s],columns = STD_STRUCTURAL_EVALUATION_COLUMNS,index= [self.current_generation])
		return pd.concat([ev_g,ev_s,ev_f],1)
	
		
	def eaggregatestats(self,stats,quantile,quantile_next_bound,col):
		quantile_min = min(quantile,quantile_next_bound)
		quantile_max = max(quantile,quantile_next_bound)
		to_aggregate = (pd.DataFrame(stats[col]>=quantile_min) & pd.DataFrame(stats[col]>=quantile_max))[col]
		ans = pd.DataFrame([list(stats.loc[to_aggregate,STD_USEFULL_STATS_AGGREGATION].mean())],columns =  STD_USEFULL_STATS_AGGREGATION, index = [self.current_generation])
		ans.loc[:,col] = quantile
		return ans

	
	def egetdata(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.egetdata(self.p_simulation_c,self.data)
		return self.data
	
	def egetdataasdf(self):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.egetdata(self.p_simulation_c,self.data)
		return self.datatopandasdf(self.data)

	def esavesimulation(self,name_backup_simulation, name_backup_prng , name_backup_params):
		modevoevo_c.esavesimulation(self.elogger.get_full_name(name_backup_simulation),self.elogger.get_full_name(name_backup_prng),self.elogger.get_full_name(name_backup_params), self.p_simulation_c)

	def eloadsimulation(self,name_backup_simulation, name_backup_prng , name_backup_params):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eloadsimulation(self.p_simulation_c,name_backup_simulation, name_backup_prng, name_backup_params)

#old

	def esavestats(self,name_file):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esavestats(self.p_simulation_c,self.elogger.get_full_name(name_file))

	def esaveindividualstat(self,individual,name_file):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esaveindividualstat(self.p_simulation_c,individual,self.elogger.get_full_name(name_file))

	def esavegenotype(self,individual,name_file):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esavegenotype(self.p_simulation_c,individual,self.elogger.get_full_name(name_file))

	def esavephenotype(self,individual,name_file):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.esavephenotype(self.p_simulation_c,individual,self.elogger.get_full_name(name_file))

	def esaveclassifieddata(self,individual,name_file):
		if self.p_simulation_c == 0: return 0
		modevoevo_c.eclassifydataset(self.p_simulation_c,individual,self.elogger.get_full_name(name_file))

	def emakesaves(self,title,ID_indivs,savegenotypes,savephenotypes,savetables,savestats):
		if savestats:
			self.esavestats("stats.dat")
		if savegenotypes:
			for indiv in ID_indivs:
				self.esavegenotype(indiv,title+"_"+"GENOTYPE"+"_generation="+str(self.current_generation)+"_indiv="+str(indiv)+".txt")
		if savephenotypes:
			for indiv in ID_indivs:
				self.esavephenotype(indiv,title+"_"+"PHENOTYPE"+"_generation="+str(self.current_generation)+"_indiv="+str(indiv)+".txt")
		if savetables:
			for indiv in ID_indivs:
				self.esaveclassifieddata(indiv,title+"_"+"CLASSIFIED_DATA"+"_generation="+str(self.current_generation)+"_indiv="+str(indiv)+".txt")

	def emakesavesnbestindividuals(self,title,nb_best_individuals_to_save,savegenotypes,savephenotypes,savetables,savestats):
		pop_stats = self.egetstats()
		ID_indivs = self.egetNbestindividualsindex(pop_stats,nb_best_individuals_to_save)
		self.emakesaves(title,ID_indivs,savegenotypes,savephenotypes,savetables,savestats)
		
	def ecomputeindividualfitness(self, individual_index):
		if self.p_simulation_c == 0:
			return 0
		return modevoevo_c.ecomputeindividualfitness(self.p_simulation_c, individual_index)

	def ecomputeclassification(self,eindiv=0):
		return self.ecomputeindividualfitness(eindiv),self.egetindividualclassifieddataasdf(eindiv)
	
	def esetgene(self,individualIndex,geneIndex,newGene):
		if self.p_simulation_c == 0:
			return 0
		return modevoevo_c.emodifygene(self.p_simulation_c,individualIndex,geneIndex,newGene)
