#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
@author: Sergio Peignier, Christophe Rigotti, Guillaume Beslon
'''
from definitions import *
from egenomeparameters import *

class eparams:
	def __init__(self,
				 mutation_rate                         = STD_MUT_RATE,
				 gene_size                             = STD_GENE_SIZE,
				 nb_init_genes                         = STD_INIT_GENES,
				 size_data_buffer                      = STD_DATA_BUFFER,
				 size_data_point                       = STD_SIZE_DATA_POINT,
				 genom_limit_size                      = STD_GENOME_LIM_SIZE,
				 population_size                       = STD_POP_SIZE ,
				 current_generation_index              = STD_CURRENT_GENERATION_INDEX,
				 selection_pressure                    = STD_SELECT_PRESSURE,
				 kmeans_iterations                     = STD_KMEAN_ITERATIONS,
				 norm_exponent                         = STD_NORM_EXPONENT,
				 fitness_mode                          = STD_FITNESS_MODE,
				 log_file                              = STD_LOG_FILE,
				 prng_seed                             = STD_PRNG_SEED,
				 sort_by_position                      = STD_SORT_BY_POSITION,
				 saveBestIndividual                    = STD_SAVE_BEST_INDIVIDUAL,
				 non_coding_genome_fitness             = STD_NON_CODING_GENOME_FITNESS,
				 shuffle_genome                        = STD_SHUFFLE_GENOME,
				 more_than_one_point_substitution_rate = STD_CASCADE_MUT,
				 intergenic_cut                        = "",
				 mutation_laws                         = "",
				 init_positions_laws                   = "",
				 boundary_conditions                   = "",
				 gene_Element_Mut_Prob_Law             = "",
				 large_duplication_rate                = "",
				 large_deletion_rate                   = "",
				 large_translocation_rate              = "",
				 large_invertion_rate                  = "",
				 point_substitution_rate               = "",
				 point_insertion_rate                  = "",
				 point_deletion_rate                   = "",
				 unknown_dimention_random_values_generator = ""):
		'''
				 Creates a eparams instance
				 
				 :param mutation_rate: Mutation rate
				 :type  mutation_rate: float
				 :param gene_size    : Tuple size :math: '|\gamma|'
				 :type  gene_size    : int
				 :param nb_init_genes : Initial genome size :math: '|\Gamma|_{t0}'
				 :type  nb_init_genes : int
				 :param size_data_buffer : Number of objets treated
				 :type size_data_buffer  : int
				 :param size_data_point  : Number of dimensions per object
				 :type size_data_point   : int
				 :param genom_limit_size : Maximal genome size authorized
				 :type genom_limit_size  : int
				 :param population_size  : Number of individuals in population
				 :type  population_size  : int
				 :param current_generation_index : Current generation index
				 :type  current_generation_index : int
				 :param selection_pressure                    : Selection pressure s :math: `(s-1)\frac{s^{N-r}}{s^N-1}`
				 :type selection_pressure                     : int
				 :param kmeans_iterations                     : Number of kmeans iterations
				 :type kmeans_iterations                      : int
				 :param norm_exponent                         : Norm of the lp metric to compute distances
				 :type norm_exponent                          : int
				 :param fitness_mode                          : Fitness evaluation function to use
				 :type fitness_mode                           : int
				 :param log_file                              : Log file name
				 :type  log_file                              : str
				 :param prng_seed                             : Pseudo random numbers generator seed
				 :type  prng_seed                             : float
				 :param sort_by_position                      : Should sort individuals by position at each generation
				 :type  sort_by_position                      : bool
				 :param saveBestIndividual                    : Elitist strategy
				 :type  saveBestIndividual                    : bool
				 :param non_coding_genome_fitness             : Fitness of individuals without phenotype
				 :type  non_coding_genome_fitness             : float
				 :param shuffle_genome                        : Should shuffle the individuals genomes at each generation
				 :type  shuffle_genome                        : bool
				 :param more_than_one_point_substitution_rate : Should produce a cascade of mutations for point substitutions
				 :type  more_than_one_point_substitution_rate : bool
				 :param intergenic_cut                        : Intergenic cut details
				 :type  intergenic_cut                        : eintergeniccut
				 :param mutation_laws                         : Tuple elements mutations distributions
				 :type  mutation_laws                         : float
				 :param init_positions_laws                   : Tuple elements initial values distributions
				 :type  init_positions_laws                   : list
				 :param boundary_conditions                   : Tuple elements bounds
				 :type  boundary_conditions                   : list
				 :param gene_Element_Mut_Prob_Law             : Tuple elements mutation probbilities
				 :type  gene_Element_Mut_Prob_Law             : egeneElementMutProbLaw
				 :param large_duplication_rate                : Duplication rate per tuple
				 :type  large_duplication_rate                : float
				 :param large_deletion_rate                   : Deletion rate per tuple
				 :type  large_deletion_rate                   : float
				 :param large_translocation_rate              : Translocation rate per tuple
				 :type  large_translocation_rate              : float
				 :param large_invertion_rate                  : Inversion rate per tuple
				 :type  large_invertion_rate                  : float
				 :param point_substitution_rate               : Point substitution rate per tuple
				 :type  point_substitution_rate               : float
				 :param point_insertion_rate                  : Point insertion rate per tuple (not implemented) 
				 :type  point_insertion_rate                  : float
				 :param point_deletion_rate                   : Point deletion rate per tuple (not implemented) 
				 :type  point_deletion_rate                   : float
				 :param unknown_dimention_random_values_generator : Distribution to produce positions for dimensions not present in the dataset points):
				 :type  unknown_dimention_random_values_generator : estatisticalDistributionLaw
				 
			'''
		self.shuffle_genome						   = shuffle_genome
		self.genom_limit_size					   = genom_limit_size
		self.gene_size							   = gene_size
		self.nb_init_genes						   = nb_init_genes
		self.size_data_buffer					   = size_data_buffer
		self.size_data_point					   = size_data_point
		self.mutation_rate						   = mutation_rate
		self.large_duplication_rate				   = self.mutation_rate
		self.large_deletion_rate 				   = self.mutation_rate
		self.large_translocation_rate			   = 0
		self.large_invertion_rate				   = 0
		self.point_substitution_rate			   = self.mutation_rate
		self.point_insertion_rate				   = 0
		self.point_deletion_rate				   = 0
		self.sort_by_position 					   = sort_by_position
		self.saveBestIndividual					   = saveBestIndividual
		self.population_size					   = population_size
		self.current_generation_index			   = current_generation_index
		self.selection_pressure					   = selection_pressure
		self.kmeans_iterations		 			   = kmeans_iterations
		self.norm_exponent			  			   = norm_exponent
		self.log_file							   = log_file
		self.fitness_mode			 			   = fitness_mode
		self.prng_seed                             = prng_seed
		self.more_than_one_point_substitution_rate = more_than_one_point_substitution_rate
		self.non_coding_genome_fitness             = non_coding_genome_fitness

		if gene_Element_Mut_Prob_Law != "": self.gene_Element_Mut_Prob_Law = gene_Element_Mut_Prob_Law
		else:                               self.gene_Element_Mut_Prob_Law = egeneElementMutProbLaw()

		if intergenic_cut == "":            self.intergenic_cut = eintergeniccut()
		else:                               self.intergenic_cut = intergenic_cut

		if mutation_laws != "":	            self.mutation_laws = mutation_laws
		else:                               self.mutation_laws = [estatisticalDistributionLaw(emin = 0,emax = 2,law = UNIFORM),
																  estatisticalDistributionLaw(emin = 0,emax = 17,law = UNIFORM),
																  estatisticalDistributionLaw(emin = 0,emax = 1025,law = UNIFORM),
																  estatisticalDistributionLaw(emin = -1024,emax = 1025,law = UNIFORM)]

		if init_positions_laws !="":		self.init_positions_laws = init_positions_laws
		else: 								self.init_positions_laws = [estatisticalDistributionLaw(emin = 0,emax = 2,law = UNIFORM),
											    						estatisticalDistributionLaw(emin = 0,emax = 17,law = UNIFORM),
											 							estatisticalDistributionLaw(emin = 0,emax = 1025,law = UNIFORM),
																		estatisticalDistributionLaw(emin = -1024,emax = 1025,law = UNIFORM)]

		if boundary_conditions != "":       self.boundary_conditions = boundary_conditions
		else:					            self.boundary_conditions = [eboudaryCondition(-1,2),eboudaryCondition(-1,17),eboudaryCondition(-1,1025),eboudaryCondition(-1025,1025)]

		if large_duplication_rate != ""   : self.large_duplication_rate=large_duplication_rate
		if large_deletion_rate != ""      : self.large_deletion_rate = large_deletion_rate
		if large_translocation_rate != "" : self.large_translocation_rate = large_translocation_rate
		if large_invertion_rate != ""     : self.large_invertion_rate=large_invertion_rate
		if point_substitution_rate != ""  : self.point_substitution_rate = point_substitution_rate
		if point_insertion_rate != ""     : self.point_insertion_rate = point_insertion_rate
		if point_deletion_rate != ""      : self.point_deletion_rate = point_deletion_rate

		if unknown_dimention_random_values_generator != "":	 self.unknown_dimention_random_values_generator = unknown_dimention_random_values_generator
		else:                                                self.unknown_dimention_random_values_generator = estatisticalDistributionLaw(emin = self.boundary_conditions[POSITION_GENE_TYPE].emin,emax = self.boundary_conditions[POSITION_GENE_TYPE].emax,law = UNIFORM)#NAN penalty
		#print self.unknown_dimention_random_values_generator.emin,self.unknown_dimention_random_values_generator.emax
		self.check_parameters()

	def check_parameters(self):
		''' 
			Run some basic tests to check the parameters entered
			'''
		self.max_gene_size = MAX_GENE_SIZE
		if self.gene_size > self.max_gene_size:
			raise "ERROR: The gene size is too big"
			return None
		if len(self.init_positions_laws) > self.gene_size:
			self.init_positions_laws[:] = self.init_positions_laws[0:self.gene_size]
		if len(self.mutation_laws) > self.gene_size:
			self.mutation_laws[:] = self.mutation_laws[0:self.gene_size]
		if len(self.boundary_conditions) > self.gene_size:
			self.boundary_conditions[:] = self.boundary_conditions[0:self.gene_size]
		if len(self.gene_Element_Mut_Prob_Law.mut_law_gene_elmnt) > self.gene_size:
			self.gene_Element_Mut_Prob_Law = egeneElementMutProbLaw(self.gene_Element_Mut_Prob_Law.mut_law_gene_elmnt[0:self.gene_size])

		if len(self.init_positions_laws) < self.gene_size:
			raise "ERROR: not enough information to define the initial positions laws"
		if len(self.mutation_laws) < self.gene_size:
			raise "ERROR: not enough information to define the mutations laws"
		if len(self.boundary_conditions) < self.gene_size:
			raise "ERROR: not enough information to define the boundary conditions"
		if len(self.gene_Element_Mut_Prob_Law.mut_law_gene_elmnt) < self.gene_size:
			raise "ERROR: not enough information to the point mutation probabilities"

	def getastuple(self):
		''' 
			Return the parameters in a long tuple
			'''
		eparameters = (self.gene_size,
		self.nb_init_genes,
		self.size_data_buffer,
		self.size_data_point,
		self.genom_limit_size,
		self.large_duplication_rate,
		self.large_deletion_rate,
		self.large_translocation_rate,
		self.large_invertion_rate,
		self.point_substitution_rate,
		self.point_insertion_rate,
		self.point_deletion_rate,
		self.more_than_one_point_substitution_rate,
		self.population_size,
		self.current_generation_index,
		self.selection_pressure,
		self.kmeans_iterations,
		self.norm_exponent,
		self.fitness_mode,
		self.log_file,
		self.prng_seed,
		self.sort_by_position,
		self.intergenic_cut.get_this_as_tuple(),
		self.saveBestIndividual,
		tuple([element.get_this_as_tuple() for element in self.mutation_laws]),
		tuple([element.get_this_as_tuple() for element in self.init_positions_laws]),
		tuple([element.get_this_as_tuple() for element in self.boundary_conditions]),
		self.gene_Element_Mut_Prob_Law.get_this_as_tuple_of_arrays(),
		self.non_coding_genome_fitness,
		self.unknown_dimention_random_values_generator.get_this_as_tuple(),
		self.shuffle_genome)
		return eparameters
	
	def eloadfromtuple(self,parameters_to_load):
		''' 
			Load a eparameters objet from a set of parameters entered in a tuple
			:param parameters_to_load : Tuple of parameters):
			:type  parameters_to_load : tuple
			'''
		self.eparams.gene_size								= parameters_to_load[0]
		self.eparams.nb_init_genes							= parameters_to_load[1]
		self.eparams.size_data_buffer						= parameters_to_load[2]
		self.eparams.size_data_point						= parameters_to_load[3]
		self.eparams.genom_limit_size						= parameters_to_load[4]
		self.eparams.large_duplication_rate					= parameters_to_load[5]
		self.eparams.large_deletion_rate					= parameters_to_load[6]
		self.eparams.large_translocation_rate				= parameters_to_load[7]
		self.eparams.large_invertion_rate					= parameters_to_load[8]
		self.eparams.point_substitution_rate				= parameters_to_load[9]
		self.eparams.point_insertion_rate					= parameters_to_load[10]
		self.eparams.point_deletion_rate					= parameters_to_load[11]
		self.eparams.more_than_one_point_substitution_rate	= parameters_to_load[12]
		self.eparams.population_size						= parameters_to_load[13]
		self.eparams.current_generation_index				= parameters_to_load[14]
		self.eparams.selection_pressure						= parameters_to_load[15]
		self.eparams.kmeans_iterations						= parameters_to_load[16]
		self.eparams.norm_exponent							= parameters_to_load[17]
		self.eparams.fitness_mode							= parameters_to_load[18]
		self.eparams.log_file								= parameters_to_load[19]
		self.eparams.prng_seed								= parameters_to_load[20]
		self.eparams.sort_by_position						= parameters_to_load[21]
		self.eparams.intergenic_cut				        	= self.eparams.intergenic_cut.get_this_from_tuple(parameters_to_load[22])
		self.eparams.saveBestIndividual						= parameters_to_load[23]
		self.cat_objects_tuple_conversion(self.eparams.mutation_laws,parameters_to_load[24])
		self.cat_objects_tuple_conversion(self.eparams.init_positions_laws,parameters_to_load[25]),
		self.cat_objects_tuple_conversion(self.eparams.boundary_conditions,parameters_to_load[26]),
		self.eparams.gene_Element_Mut_Prob_Law.get_this_from_tuple_of_arrays(parameters_to_load[27]),

		self.eparams.non_coding_genome_fitness				= parameters_to_load[28]
		self.eparams.unknown_dimention_random_values_generator.get_this_from_tuple(parameters_to_load[29]),

		self.eparams.shuffle_genome							= parameters_to_load[30]


def generatestandardparameter(dataset,
							  nb_max_clusters,
							  prng_seed               = STD_PRNG_SEED,
							  stream_size_p           = STD_WINDOW_SIZE ,
							  transition_matrix_c_nc  = NONABSORBINGSTATE,
							  transition_matrix_init  = BOTHABSORBING,
							  tuples_mutations_prob   = STD_GENE_MUTATION_PROBABILITIES,
							  nb_init_genes           = STD_INIT_GENES ,
							  mutation_rate           = STD_MUT_RATE,
							  population_size         = STD_POP_SIZE,
							  selection_pressure      = STD_SELECT_PRESSURE,
							  extraparams             = {}): 
	''' 
		Generates a standard parameter object, the main parameters can be tuned directly
		
		:param mutation_rate: Mutation rate
		:type  mutation_rate: float
		:param nb_init_genes : Initial genome size :math: '|\Gamma|_{t0}'
		:type  nb_init_genes : int
		:param stream_size_p: sliding window size
		:type size_data_buffer  : float
		:param population_size  : Number of individuals in population
		:type  population_size  : int
		:param selection_pressure                    : Selection pressure s :math: `(s-1)\frac{s^{N-r}}{s^N-1}`
		:type selection_pressure                     : int
		:param prng_seed                             : Pseudo random numbers generator seed
		:type  prng_seed                             : float
		:param transition_matrix_c_nc                : Functional_non-functional transition matrix
		:type  transition_matrix_c_nc                : list
		:param transition_matrix_init                : Initial functional_non-functional transition matrix
		:type  transition_matrix_init                : list
		:param tuples_mutations_prob                 : Mutation probabilities for each tuple element
		:type  tuples_mutations_prob                 : list
		:param extraparams                           : More parameters to change
		:type  extraparams							 : dict
		:param dataset                               : Dataset to deal with
		:type  extraparams							 : edataset
		:param nb_max_clusters                       : Cmax parameter
		:type  nb_max_clusters				     	 : int
	
		'''
							  
	position_lims = dataset.df[dataset.features].abs().max().max()
	cluster_max   = nb_max_clusters
	dim_max       = len(dataset.features)
	if (stream_size_p <= 1.0):
		stream_size   = int(stream_size_p*len(dataset.df))
	else:
		stream_size   = int(stream_size_p)
	stream_size = min(stream_size,len(dataset.df.index))
	ml = generateeasymutationlawarray(dataset,nb_max_clusters, transition_matrix_c_nc =  transition_matrix_c_nc,  cnc_law =  TRANSITION_MATRIX)
	il = generateeasymutationlawarray(dataset,nb_max_clusters, transition_matrix_c_nc =  transition_matrix_init , cnc_law =  TRANSITION_MATRIX)
	gmp = egeneElementMutProbLaw(tuples_mutations_prob)
	bc = [eboudaryCondition(-1,2),
		  eboudaryCondition(-1,cluster_max ),
		  eboudaryCondition(-1,dim_max ),
		  eboudaryCondition(-position_lims-1,position_lims+1)]

	params = {"size_data_buffer"           : stream_size,
	          "mutation_laws"              : ml,
	          "init_positions_laws"        : il,
	          "boundary_conditions"        : bc,
	          "gene_Element_Mut_Prob_Law"  : gmp,
	          "prng_seed"                  : prng_seed,
	          "nb_init_genes"              : nb_init_genes,
	          "mutation_rate"              : mutation_rate,
	          "population_size":population_size,
	          "selection_pressure":selection_pressure}
	          	
	params.update(extraparams)
	params = eparams(**params)
	return params


def generateestandardparametersfrompandasDF(dataframe,
											cmax,
											sliding_sample_size    = STD_WINDOW_SIZE ,
											selection_pressure     = STD_SELECT_PRESSURE,
											init_genome_size       = STD_INIT_GENES ,
											mutation_rate          = STD_MUT_RATE,
											population_size        = STD_POP_SIZE,
											transition_matrix_c_nc = NONABSORBINGSTATE,
											transition_matrix_init = BOTHABSORBING,
											tuples_mutations_prob  = STD_GENE_MUTATION_PROBABILITIES,
											prng_seed              = STD_PRNG_SEED,
											elitism                = STD_SAVE_BEST_INDIVIDUAL,
											extraparams            = {}):
	
	position_lims = dataframe.abs().max().max()
	cluster_max   = cmax
	dim_max       = len(dataframe.columns)
	if (sliding_sample_size <= 1.0):
		stream_size   = int(sliding_sample_size*len(dataframe))
	else:
		stream_size   = int(sliding_sample_size)
	stream_size = min(stream_size,len(dataframe.index))

	ml=[estatisticalDistributionLaw(law = TRANSITION_MATRIX, transition_matrix = transition_matrix_c_nc),
		estatisticalDistributionLaw(emin = -cluster_max/2 - 1,emax = cluster_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -dim_max/2 - 1,emax =  dim_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -position_lims-1,emax = position_lims+1,law = UNIFORM)]

	il=[estatisticalDistributionLaw(law = TRANSITION_MATRIX, transition_matrix = transition_matrix_init),
		estatisticalDistributionLaw(emin = -cluster_max/2 - 1,emax = cluster_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -dim_max/2 - 1,emax =  dim_max/2 + 1,law = UNIFORM),
		estatisticalDistributionLaw(emin = -position_lims-1,emax = position_lims+1,law = UNIFORM)]

	gmp = egeneElementMutProbLaw(tuples_mutations_prob)

	bc = [eboudaryCondition(-1,2),
		  eboudaryCondition(-1,cluster_max),
		  eboudaryCondition(-1,dim_max),
		  eboudaryCondition(-position_lims-1,position_lims+1)]

	params = {"size_data_buffer"           : int(stream_size),
	          "mutation_laws"              : ml,
	          "init_positions_laws"        : il,
	          "boundary_conditions"        : bc,
	          "gene_Element_Mut_Prob_Law"  : gmp,
	          "prng_seed"                  : float(prng_seed),
	          "nb_init_genes"              : int(init_genome_size),
	          "mutation_rate"              : float(mutation_rate),
	          "population_size"            : int(population_size),
	          "selection_pressure"         : float(selection_pressure),
	          "saveBestIndividual"         : int(elitism)}
	          	
	params.update(extraparams)
	params = eparams(**params)
	return params
	