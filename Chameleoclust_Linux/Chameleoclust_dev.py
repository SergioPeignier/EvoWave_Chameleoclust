#!/usr/bin/python
import sys
import time
import os
import pandas as pd
from modevoevo import *
from evoevomanager import *
from genotypemanager import *
from definitions import *
from usualdatasets import *
from pandas import DataFrame
import random
import numpy
import math
import copy

class Chameleoclust_:
    def __init__(self,
                 eparams_input="",
                 elogger="",
                 nb_individuals_for_evoevo_test=300,
                 static_evaluator=1,
                 evoevo_evaluator=1,
                 save_best_individual=0,
                 save_txt_files={0: {"savegenotypes": 0, "savephenotypes": 0, "savetables": 0, "savestats": 0,
                                     "nb_best_individuals_to_save": 0}},
                 stat_computation_period=100,
                 aggregation_dict={"fitness": {"quantile": [1], "bounds": [-0.25]}},
                 best_individual_evaluation_params={"true_dataframe": [], "g": 1, "f": 1, "s": 1,
                                                    "threshold_cluster_validity": 0.0,
                                                    "noise_hidden_cluster_id": [666]},
                 evaluator_stream_size=10):
        self.best_individual_evaluation_params = best_individual_evaluation_params
        self.aggregation_dict = aggregation_dict
        self.save_txt_files = save_txt_files
        self.save_best_individual = save_best_individual
        self.eparams = eparams_input
        self.evoevo_evaluator = evoevo_evaluator
        self.stat_computation_period = stat_computation_period
        self.simulation = esimulation(self.eparams, elogger, save_info_in_log_file = 0)
        if self.evoevo_evaluator:
            self.eparams_evoevo = copy.deepcopy(eparams_input)
            self.eparams_evoevo.saveBestIndividual = 0
            self.eparams_evoevo.population_size = nb_individuals_for_evoevo_test
            self.evoevo_manager = EvoEvoManager(self.eparams_evoevo, self.simulation.elogger, save_info_in_log_file = 0)

        self.static_evaluator = static_evaluator
        if self.static_evaluator:
            self.eparams_evaluator = copy.deepcopy(eparams_input)
            self.eparams_evaluator.population_size = 1
            self.eparams_evaluator.saveBestIndividual = 0
            self.eparams_evaluator.size_data_buffer = evaluator_stream_size
            self.simulation_evaluator = esimulation(self.eparams_evaluator, self.simulation.elogger, save_info_in_log_file = 0)

    def should_compute_stats(self):
        return not (self.simulation.current_generation % self.stat_computation_period) - 1 or self.stat_computation_period == 1

    def compute_local_stats(self):
        self.simulation.egetstats()
        self.simulation.egetstatsasdf()
        for k in self.aggregation_dict.keys():
            self.simulation.eaddtoaggregatedstats(col=k, **self.aggregation_dict[k])

    def compute_stats(self, force=0):
        if self.should_compute_stats() or force:

            self.simulation.egetstats()
            self.simulation.egetstatsasdf()
            print "local_evaluation",self.simulation.statsdf.mean()
            make_txt_saves = self.save_txt_files.keys()[0]

            if self.static_evaluator:
                self.simulation_evaluator.current_generation = self.simulation.current_generation
                self.simulation_evaluator.esetpopulationgenome(egenome=self.simulation.egetbestindividualgenome())
                self.simulation_evaluator.ecomputefitnesses()
                self.simulation_evaluator.egetstats()
                self.simulation_evaluator.egetstatsasdf()
                self.simulation_evaluator.eaggregatebestindividualevaluation(**self.best_individual_evaluation_params)
                print "overall_evaluation",self.simulation_evaluator.statsdf.mean()
                if make_txt_saves:
                    self.save_txt_files[make_txt_saves]["title"] = "STATIC"
                    self.simulation.emakesavesnbestindividuals(**self.save_txt_files[make_txt_saves])

            if self.evoevo_evaluator:
                self.evoevo_manager.ecomputeNchildrenfromparentalgenome(
                    parentalgenome=self.simulation.egetbestindividualgenome())
                self.evoevo_manager.eaddtoaggregatedstats(generation=self.simulation_evaluator.current_generation)

            if self.save_best_individual:
                self.save_best_individual_function()

            for k in self.aggregation_dict.keys():
                self.simulation.eaddtoaggregatedstats(col=k, **self.aggregation_dict[k])

            if make_txt_saves:
                self.save_txt_files[make_txt_saves]["title"] = "STREAM"
                self.simulation.emakesavesnbestindividuals(**self.save_txt_files[make_txt_saves])
                self.simulation.eaggregatebestindividualevaluation(**self.best_individual_evaluation_params)

    def eiterate(self):
        self.simulation.eiterate(1)
        #self.compute_stats()

    def setdata(self, data, replace):
        self.simulation.esetdata(data, replace)
        if self.evoevo_evaluator:
            self.evoevo_manager.setdata(data, replace)

    def erun(self, 
             data,
             nb_training_generations_per_update=1,
             nb_total_generations=5000,
             proportion_2_update=1.0):
        dataset_length = len(data)
        stream_size = self.simulation.eparams.size_data_buffer
        random.seed(self.simulation.eparams.prng_seed)
        numpy.random.seed(self.simulation.eparams.prng_seed)
        if self.static_evaluator:
            self.simulation_evaluator.esetdata(data[0:self.simulation_evaluator.eparams.size_data_buffer], 1)
        self.setdata(data[0:stream_size], 1)
        self.simulation.ecomputefitnesses()
        self.compute_stats(1)
        if proportion_2_update < 0:
            proportion_2_update = 1.0
        if proportion_2_update <= 1.0:
            update_size = int(proportion_2_update * stream_size)
        else:
            update_size = int(proportion_2_update)
        nb_updates = 0
        i = 0
        exit = 0
        while 1:
            j = 0
            while j < nb_training_generations_per_update:
                self.eiterate()
                j += 1
                i += 1
                if i == nb_total_generations:
                    exit = 1
                    break
            if exit:
                break
            pos1 = (dataset_length + ((nb_updates * update_size + stream_size) % dataset_length)) % dataset_length 
            pos2 = (dataset_length + (((nb_updates + 1) * update_size + stream_size) % dataset_length )) % dataset_length 
            if pos1 == pos2:
                data_local = data
            if pos1 < pos2:
                data_local = data[pos1:pos2]
            if pos2 < pos1:
                data_local = data[pos1:dataset_length] + data[0:pos2]
            self.setdata(data_local,0)
            nb_updates += 1
            if exit:
                break
        self.compute_stats(1)
        self.save_aggregations()
        if not self.save_best_individual:
            self.save_best_individual_function()

    def save_aggregations(self):
        self.simulation.elogger.save_object_to_file(self.simulation.aggregatedstats, "aggregated_stream_results_" + str(self.simulation.current_generation) + ".PICKLE")
        if self.evoevo_evaluator:
            self.simulation.elogger.save_object_to_file(self.evoevo_manager.aggregatedstats, "EvoEvo_results_" + str(self.simulation.current_generation) + ".PICKLE")
        if self.static_evaluator:
            self.simulation.elogger.save_object_to_file(self.simulation_evaluator.aggregatedstats, "aggregated_static_results_" + str(self.simulation.current_generation) + ".PICKLE")

    def save_best_individual_function(self):
        best_clust_mod = 0
        best_phenotype = 0
        best_genome = 0
        if self.static_evaluator:
            best_clust_mod = self.simulation_evaluator.egetindividualclassifieddataasdf(0)
            best_phenotype = self.simulation_evaluator.egetindividualphenotypeasdataframe(0)
            best_genome = self.simulation_evaluator.egetindividualgenome(0)
        else:
            best_clust_mod = self.simulation.egetindividualclassifieddataasdf(self.simulation.egetbestindividualindex()[0])
            best_phenotype = self.simulation.egetindividualphenotypeasdataframe(self.simulation.egetbestindividualindex()[0])
            best_genome = self.simulation.egetindividualgenome(self.simulation.egetbestindividualindex()[0])
        self.simulation.elogger.save_object_to_file(best_clust_mod, "BEST_INDIVIDUAL_SUBSPACE_MODEL_" + str(self.simulation.current_generation) + ".PICKLE")
        self.simulation.elogger.save_object_to_file(best_phenotype, "BEST_INDIVIDUAL_PHENOTYPE_" + str(self.simulation.current_generation) + ".PICKLE")
        self.simulation.elogger.save_object_to_file(best_genome, "BEST_INDIVIDUAL_GENOTYPE_" + str(self.simulation.current_generation) + ".PICKLE")
