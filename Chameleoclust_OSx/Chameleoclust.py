#!/usr/bin/python
from modevoevo import esimulation,generateestandardparametersfrompandasDF
from pandas import DataFrame
from numpy import array,isnan
from copy import deepcopy
import sys

class Chameleoclust:
    def __init__(self,
        sliding_sample_size = 0.1,
        selection_pressure = 0.5,
        init_genome_size = 200,
        mutation_rate = 0.00142,
        population_size = 300,
        number_of_generations = 5000,
        elitism = 1,
        cmax = 10,
        seed = 0,
        gene_pseudogene_transition_matrix = [[0,1],\
                                             [1,0]]
        ):
        """
        Creates a new Chameleoclust++ instance
        
        :param sliding_sample_size: Size of the sliding sample, if
        0<sliding_sample_size <= 1 the size is relative to the dataset size.
        Otherwise it is the absolute size.
        :param selection_pressure: Selection pressure, if close to 1
        the selection presssure will be low, if close to 0 il will be very strong.
        :param init_genome_size: Initial genome size, the algorithm is not very
        sensitive to this parameter, specially if it is taken high enough
        (e.g., >50).
        :param mutation_rate: Mutation rate, suitable mutation rates are
        usually in between 0.01 and 0.0001. If it is chosen to low, evolution
        is very slow, if it is taken to high, organisms cannot evolve properly.
        :population_size: Population size, higher population sizes usually lead
        to better results, however it also leads to higher runtimes. In between
        100 and 300 individuals usually lead to good results and reasonable
        runtimes.
        :param number_of_generations: Number of generations, during how many
        generations should we evolve individuals.
        :param elitism: Elistism, if this parameter is set to True, then one
        copy of the best individual is copied whithout any changes  and is
        inserted into the next generation.
        :param cmax: Cmax, maximal number of distinct core-points that
        induviduals can effectively generate.
        :param seed: Seed, fixed seed for the pseudo random numbers generator.
        :param gene_pseudogene_transition_matrix: Matrix for Gene <-> Pseudogenes
        transition probabilities. The first element in the matrix corresponds to
        genes (functional) and the second to Pseudogenes (non-functional)

        :type sliding_sample_size: int or float
        :type selection_pressure: float in [0,1[
        :type init_genome_size: int >= 1
        :type mutation_rate: float in [0,1[
        :type population_size: int >= 1
        :type number_of_generations: int >= 1
        :type elitism: boolean
        :type cmax: int >= 1
        :type seed: float
        :type gene_pseudogene_transition_matrix: list
        """
        self.sliding_sample_size = sliding_sample_size
        self.selection_pressure = selection_pressure
        self.init_genome_size = init_genome_size
        self.mutation_rate = mutation_rate
        self.population_size = population_size
        self.number_of_generations = number_of_generations
        self.elitism = elitism
        self.cmax = cmax
        self.seed = seed
        self.data = None
        self.eparams = None
        self.eparams_predictor = None
        self.simulation = None
        self.predictor = None
        self._population_stats_up_to_data = False
        self._best_individual_index = None
        self.gene_pseudogene_transition_matrix = gene_pseudogene_transition_matrix
        self.stats_evolution = DataFrame(columns = ["best_fitness", 
                                                        "coding_ratio_best",
                                                        "genome_length_best",
                                                        "mean_fitness_quantile", 
                                                        "coding_ratio_quantile",  
                                                        "mean_genome_length_quantile"])


    def _generate_chameleo(self,X): 
        self.eparams = generateestandardparametersfrompandasDF(dataframe              = X,
                                                               cmax                   = self.cmax,
                                                               sliding_sample_size    = self.sliding_sample_size,
                                                               selection_pressure     = self.selection_pressure,
                                                               init_genome_size       = self.init_genome_size,
                                                               mutation_rate          = self.mutation_rate,
                                                               population_size        = self.population_size,
                                                               prng_seed              = self.seed,
                                                               elitism                = self.elitism,
                                                               transition_matrix_c_nc = self.gene_pseudogene_transition_matrix)
        self.eparams_predictor = deepcopy(self.eparams)
        self.eparams_predictor.population_size = 1
        self.eparams_predictor.saveBestIndividual = 0
        self.eparams_predictor.size_data_buffer = 1
        self.simulation = esimulation(self.eparams, save_info_in_log_file = 0)
        self.predictor  = esimulation(self.eparams_predictor, save_info_in_log_file = 0)


    def set_params(self, **params):
        """
        Modifies the values of the parameters set as arguments.
        """
        for name in params:
            if not hasattr(self, name):
                raise ValueError('Invalid parameter %s for estimator %s.'
                                'Check the list of available parameters '
                                'with `KymerClust.get_params().keys()`.' %
                                (name, self))
            else:
                setattr(self, name, params[name])

    def get_params(self, deep=True):
        """
        Returns values of the desired parameters
        """
        params = {param: getattr(self, param) for param in self._parameters}
        return params

    def _table_to_data(self,X,y=None):
        features = X.columns
        self.data = []
        i = 0
        while i < len(X):
            self.data.append([])
            if y is None:
                self.data[i].append(-1)
            else:
                if not isinstance(y[0],int):
                    print("Warning: only integer labels are accepted. The label won't be used")
                    self.data[i].append(-1)
                else:
                    self.data[i].append(y[i])
            self.data[i].append(-1)
            self.data[i] += [[j,e] for j,e in enumerate(X.iloc[i,range(len(features))]) if not isnan(float(e))]
            i+=1
        return self.data

    def _compute_evolution_stats_from_df(self,df,quantile):
        sorted_df = df.sort_values("fitness", ascending = False)
        sorted_df = sorted_df[['fitness', 'coding_ratio', 'genome_length']]
        best_individual = sorted_df.iloc[0,:]
        best_quantile   = sorted_df.iloc[0:int(quantile*self.population_size),:]
        self.stats_evolution.loc[self.simulation.current_generation] = list(best_individual) + list(best_quantile.mean())

    def compute_local_stats(self, quantile=0.1):
        self.simulation.egetstats()
        stats = self.simulation.egetstatsasdf()
        self._compute_evolution_stats_from_df(stats, quantile)

    def _get_stats(self):
        self._update_population_stats()
        return self.stats_evolution

    def _get_best_individual_model(self):
        self._update_population_stats()
        return self.simulation.egetindividualclassifieddataasdf(self._best_individual_index)

    def _get_best_individual_phenotype(self):
        self._update_population_stats()
        return self.simulation.egetindividualphenotypeasdataframe(self._best_individual_index).astype(float)

    def fit(self, X, y=None, nb_training_generations_per_update = 1, proportion_2_update = 1.0):
        """
        Feeds Chameleoclust++ organisms with the dataset X and let them evolve
        :param X: Dataset.
        :param y: Classes.
        :param nb_training_generations_per_update: Number of generations to
        perform at each update of the data sample.
        :param proportion_2_update: Proportion of the data sample to update.
        :type X: Numpy array 
        :type y: Numpy array or list
        :type nb_training_generations_per_update: int >= 1
        :type proportion_2_update: 0 <= float <=1
        """
        if X is None: return None
        if self.simulation is None: self._generate_chameleo(X)
        self._table_to_data(X, y)
        dataset_length = len(self.data)
        stream_size = self.simulation.eparams.size_data_buffer        
        self.simulation.esetdata(self.data[0:stream_size], replace = 1)
        self.simulation.ecomputefitnesses()
        self._update_population_stats()
        if nb_training_generations_per_update <=0 : nb_training_generations_per_update = 1
        if proportion_2_update < 0: proportion_2_update = 1.0
        if proportion_2_update <= 1.0: update_size = int(proportion_2_update * stream_size)
        else: update_size = int(proportion_2_update)
        nb_updates = 0
        i = 0
        exit = 0
        while 1:
            sys.stdout.write('\r'+'generation: '+`nb_updates`)
            sys.stdout.flush()
            j = 0
            while j < nb_training_generations_per_update:
                self._eiterate()
                j += 1
                i += 1
                if i == self.number_of_generations:
                    exit = 1
                    break
            if exit:
                break
            pos1 = (dataset_length + ((nb_updates * update_size + stream_size) % dataset_length)) % dataset_length 
            pos2 = (dataset_length + (((nb_updates + 1) * update_size + stream_size) % dataset_length )) % dataset_length 
            if pos1 == pos2:
                data_local = self.data
            if pos1 < pos2:
                data_local = self.data[pos1:pos2]
            if pos2 < pos1:
                data_local = self.data[pos1:dataset_length] + self.data[0:pos2]
            self.simulation.esetdata(data_local, replace = 0)
            nb_updates += 1
            if exit:
                break
        self._update_population_stats()

    def _fit_local(self, X, y=None,nb_iterations = 1):
        """
        Feeds Chameleoclust++ organisms with the dataset X and let them evolve
        :param X: Dataset.
        :param y: Classes.
        :param nb_iterations: Number of generations to perform on this dataset
        :type X: Numpy array 
        :type y: Numpy array or list
        :type nb_iterations: int >= 1

        """
        if X is None : return None
        else: self._table_to_data(X, y)
        if self.simulation is None: 
            self._generate_chameleo(X)
        sliding_sample_size = self.simulation.eparams.size_data_buffer
        if sliding_sample_size < len(X):
            print("WARNING: len(X) > sliding_sample_size\nOnly the first "+str(sliding_sample_size)+" elements will be used")
            X = X.iloc[0:sliding_sample_size,:]
            if y is None: y = y[0:sliding_sample_size]
        self._table_to_data(X,y)
        self.simulation.esetdata(self.data[0:sliding_sample_size], replace = 0)
        for _ in xrange(nb_iterations):
            self._eiterate()
            
    def _eiterate(self):
        self.simulation.eiterate(1)
        self._population_stats_up_to_data = False
        self._best_individual_index = None

    def _update_population_stats(self):
        if not self._population_stats_up_to_data:
            self.simulation.egetstats()
            self._best_individual_index = self.simulation.egetbestindividualindex()[0]
            self._population_stats_up_to_data = True

    def predict(self, X):
        """
        Outputs the cluster-membership for each object in the dataset set as
        parameter
        :param X: Dataset
        :returns: Cluster-membership
        :type: Numpy array
        :rtype: Numpy array
        """
        self._update_population_stats()
        if X is None: return None
        else: self._table_to_data(X)
        if self.simulation is None: return None
        self.predictor.esetpopulationgenome(egenome=self.simulation.egetbestindividualgenome())
        y = []
        for i in xrange(len(X)):
            self.predictor.esetdata(self.data[i:i+1],replace = 1)
            y += list(self.predictor.egetindividualclassifieddataasdf(0)["cluster_found"])
        return array(y, dtype=object)

    def _predict_local(self):
        """
        Returns cluster-membership of the data sample
        :returns: Cluster-membership of the data sample
        :rtype: Numpy array
        """
        self._update_population_stats()
        if self.simulation is None: return None
        return array(self._get_best_individual_model()["cluster_found"],dtype=object)

    def _fit_predict_local(self,X, y=None,nb_iterations = 1):
        """
        Calls _fit_local and then _predict_local see these functions for
        more information
        """
        self._fit_local(X, y,nb_iterations)
        return self._predict_local()

    def fit_predict(self, X, y=None):
        """
        Calls fit and predict functions, seed both functions for more details
        """
        self.fit(X, y)
        return self.predict(X)

