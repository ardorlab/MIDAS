import gc
import os
import sys
import copy
import math
import time
import h5py
import yaml
import numpy
import pickle
import shutil
import random
from multiprocessing import Pool
from midas.utils import fitness
from midas.utils.solution_types import evaluate_function,Unique_Solution_Analyzer,test_evaluate_function
from midas.utils.metrics import Optimization_Metric_Toolbox

"""
This file is for storing all the classes and methods specifically related to
performing an optimization via the genetic algorithm.

NOTE: Mutation is currently stored in this file. However, it's possible that it
will be moved to its own unique file, since it's planned that the mutation
methods will form how other optimization methods alter their solutions.

Written by Brian Andersen. 9/1/2019
"""

class Random_Solution(object):
    """
    Class for performing random solution generation

    Parameters:
        solution: Class
            Contains the genome, fitness, and optimization
            objective scores for solutions to the optimization problem.
        population: Class
            Class that contains the population size and stores the current
            solutions in the parent and child populations.
        file_settings: Dictionary
            The settings file read into the optimization. Carried through because
            some information needed to be carried along, but pickling all of the
            information didn't seem like a good way to carrty it thorugh the optimization.
    """
    def __init__(self, solution,
                 population,
                 fitness,
                 num_procs,
                 file_settings):

        self.solution = solution
        self.population = population
        self.fitness = fitness
        self.num_procs = num_procs
        self.file_settings = file_settings
       

    def generate_initial_solutions(self,name):
        """
        Function for creating the initial solutions in the optimization.
        """
        foo = self.solution()
        foo.name = f"{name}"
        foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        foo.add_additional_information(self.file_settings)
        #scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        if foo.fixed_genome:
            foo.new_generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                       self.file_settings['optimization']['fixed_groups'])
        #   foo.genome = scrambler.reproduce(foo.genome)
        else:
            foo.generate_initial(self.file_settings['genome']['chromosomes'])

        return foo

    def main_in_parallel(self):
        """
        Performs optimization using a genetic algorithm in with
        parallel computations

        Parameters: None

        Written by Brian Andersen. 1/9/2020
        """
        opt = Optimization_Metric_Toolbox()

        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        loading_pattern_tracker = open("loading_patterns.txt", 'w')
        loading_pattern_tracker.close()

        all_value_count = 0
        all_values = open('all_value_tracker.txt','w')
        all_values.close()
        for i in range(self.population.size):
            foo = self.generate_initial_solutions(f'solution_{i}')
            self.population.parents.append(foo)

        pool = Pool(processes=self.num_procs)
        self.population.parents = pool.map(evaluate_function, self.population.parents)
        print('finished solutions...')
        all_values = open('all_value_tracker.txt','a')
        for sol in self.population.parents:
            #print(sol.name)
            if not all_value_count:
                for param in sol.parameters:
                    all_values.write(f"{param},    ")
                all_values.write(f"Fitness,    ")
                all_values.write('\n')
            all_values.write(f"{all_value_count},    {sol.name},    ")
            for param in sol.parameters:
                all_values.write(f"{sol.parameters[param]['value']},    ")
                #print(f"{sol.parameters[param]['value']},    ")
            solList = self.fitness.calculate([sol])
            all_values.write(f"{solList[0].fitness},    ")
            all_values.write('\n')
            all_value_count += 1


        #self.population.parents = pool.map(self.crud.evaluator,self.population.parents)
        #self.population.children = pool.map(self.crud.evaluator,self.population.children)
        # save all param before selection perform
        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()
        all_values.close()

    def main_in_serial(self):
        """
        Performs optimization using a genetic algorithm in serial.
        Done because for some reason the neural network stuff seems to be
        breaking with parallel.

        Parameters: None

        Written by Brian Andersen 1/9/2020
        """
        opt = Optimization_Metric_Toolbox()

        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        loading_pattern_tracker = open("loading patterns.txt", 'w')
        loading_pattern_tracker.close()

        for i in range(self.population.size):
            foo = self.generate_initial_solutions(f'solution_{i}')
            foo.evaluate()
            self.population.parents.append(foo)

        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()

        opt.plotter()

    def cleanup(self):
        """
        Deletes solution results that are no longer relevant to the optimization.

        Written by Brian Andersen 10/28/2020.
        """
        gc.collect()
        if self.perform_cleanup:
            if self.generation.current <= self.number_generations_post_cleanup:
                pass
            else:
                for cl in range(self.population.size+20):
                    file_name = f'initial_parent_{cl}'
                    if os.path.isdir(file_name):
                        delete_file = True
                        for parent,child in zip(self.population.parents,self.population.children):
                            if parent.name == file_name:
                                delete_file = False
                                break
                            elif child.name == file_name:
                                delete_file = False
                                break
                        if delete_file:
                            if os.path.isfile(f"{file_name}/{file_name}_sim.out"):
                                os.remove(f"{file_name}/{file_name}_sim.out")
                    file_name = f'initial_child_{cl}'
                    if os.path.isdir(file_name):
                        delete_file = True
                        for parent,child in zip(self.population.parents,self.population.children):
                            if parent.name == file_name:
                                delete_file = False
                                break
                            elif child.name == file_name:
                                delete_file = False
                                break
                        if delete_file:
                            if os.path.isfile(f"{file_name}/{file_name}_sim.out"):
                                os.remove(f"{file_name}/{file_name}_sim.out")
                for clgen in range(self.generation.current):
                    for clpop in range(self.population.size+20):
                        file_name = f"child_{clgen}_{clpop}"
                        if os.path.isdir(file_name):
                            delete_file = True
                            for parent,child in zip(self.population.parents,self.population.children):
                                if parent.name == file_name:
                                    delete_file = False
                                    break
                                elif child.name == file_name:
                                    delete_file = False
                                    break
                            if delete_file:
                                if os.path.isfile(f"{file_name}/{file_name}_sim.out"):
                                    os.remove(f"{file_name}/{file_name}_sim.out")
