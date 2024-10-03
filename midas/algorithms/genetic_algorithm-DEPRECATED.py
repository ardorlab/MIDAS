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
import numpy as np
from multiprocessing import Pool
from midas.utils import fitness
from midas.utils.solution_types import evaluate_function,Unique_Solution_Analyzer
from midas.utils.metrics import Optimization_Metric_Toolbox

"""
This file is for storing all the classes and methods specifically related to
performing an optimization via the genetic algorithm.

NOTE: Mutation is currently stored in this file. However, it's possible that it
will be moved to its own unique file, since it's planned that the mutation
methods will form how other optimization methods alter their solutions.

Written by Brian Andersen. 9/1/2019
"""

class Genetic_Algorithm(object):
    """
    Class for performing optimization through a genetic algorithm

    Parameters:
        solution: Class
            Contains the genome, fitness, and optimization
            objective scores for solutions to the optimization problem.
        population: Class
            Class that contains the population size and stores the current
            solutions in the parent and child populations.
        generation: Class
            Keeps track of the current and total number of generations that
        reproduction: Class
            Class for performing all actions related to producing new solutions
            for the optimization.
        selection: Class _evaluate
            different solutions in the optimization.
        file_settings: Dictionary
            The settings file read into the optimization. Carried through because
            some information needed to be carried along, but pickling all of the
            information didn't seem like a good way to carrty it thorugh the optimization.

    Written by Brian Andersen. 1/9/2020
    """
    def __init__(self, solution,
                 population,
                 generation,
                 reproduction,
                 selection,
                 num_procs,
                 inp_obj,
                 file_settings):

        ## Store parameters in object
        self.solution = solution
        self.population = population
        self.generation = generation
        self.reproduction = reproduction
        self.selection = selection
        self.num_procs = num_procs
        self.file_settings = file_settings
        self.number_generations_post_cleanup = 2 #Arbitrarily chosen default.
        self.input = inp_obj
        
        if 'cleanup' in file_settings['optimization']:
            if file_settings['optimization']['cleanup']['perform']:
                self.perform_cleanup = True
            else: 
                self.perform_cleanup = False
        else:
            self.perform_cleanup = False
        
        if 'neural_network' in file_settings:
            from crudworks import CRUD_Predictor
            self.crud = CRUD_Predictor(file_settings)

    def generate_initial_solutions(self,name):
        """
        Function for creating the initial solutions in the optimization.
        """
        foo = self.solution(self.input) # solution_type is assigned, but not called, before this step
        foo.name = f"{name}"
        foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        #foo.add_additional_information(self.file_settings) #!TODO: redundant, remove.
        #scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        if foo.fixed_genome:
            foo.new_generate_initial_fixed(self.file_settings['genome']['chromosomes']['fuel'],
                                       self.file_settings['optimization']['fixed_groups'])
        #   foo.genome = scrambler.reproduce(foo.genome)
        else:
            foo.generate_initial(self.file_settings['genome']['chromosomes']['fuel'])

        return foo

    def main_in_parallel(self):
        """
        Performs optimization using a genetic algorithm with
        parallel computations

        Parameters: None

        Written by Brian Andersen. 1/9/2020
        """
        opt = Optimization_Metric_Toolbox()

        # Generate file. Clears existing.
        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        # Generate file. Clears existing.
        loading_pattern_tracker = open("loading_patterns.txt", 'w')
        loading_pattern_tracker.close()

        # Generate file. Clears existing.
        all_value_count = 0
        all_values = open('all_value_tracker.txt','w')
        all_values.close()
        
        for i in range(self.population.size):
            foo = self.generate_initial_solutions(f'initial_parent_{i}')
            self.population.parents.append(foo)

        for i in range(self.population.size):
            foo = self.generate_initial_solutions(f'initial_child_{i}')
            self.population.children.append(foo)

        pool = Pool(processes=self.num_procs)
        self.population.parents = pool.map(evaluate_function, self.population.parents)
        print('finished parents...')
        self.population.children = pool.map(evaluate_function, self.population.children)
        print('finished children...')
        all_values = open('all_value_tracker.txt','a')
        for sol in self.population.parents:
            if not all_value_count:
                for param in sol.parameters:
                    all_values.write(f"{param},    ")
                all_values.write('\n')
            all_values.write(f"{all_value_count},    {sol.name},    ")
            for param in sol.parameters:
                all_values.write(f"{sol.parameters[param]['value']},    ")
            all_values.write('\n')
            all_value_count += 1
        for sol in self.population.children:
            all_values.write(f"{all_value_count},    {sol.name},")
            for param in sol.parameters:
                all_values.write(f"{sol.parameters[param]['value']},    ")
            all_values.write('\n')
            all_value_count += 1        
        all_values.close()


        #self.population.parents = pool.map(self.crud.evaluator,self.population.parents)
        #self.population.children = pool.map(self.crud.evaluator,self.population.children)
        # save all param before selection perform
        opt.record_all_param(self.population,self.generation, flag=True)
        self.population = self.selection.perform(self.population)
        opt.check_best_worst_average(self.population.parents)
        
        opt.write_track_file(self.population, self.generation)

        #scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        #uniqueness = Unique_Solution_Analyzer(scrambler)
        for self.generation.current in range(self.generation.total):
            self.population.children = self.reproduction.reproduce(self.input, self.population.parents, 
                                                                   self.solution)
            #self.population.children = uniqueness.analyze(self.population.children)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                #solution.add_additional_information(self.file_settings)  #!TODO: redundant, remove.


            self.population.children = pool.map(evaluate_function, self.population.children)
            print('finished children...')
            all_values = open('all_value_tracker.txt','a')
            for sol in self.population.children:
                all_values.write(f"{all_value_count},   {sol.name},   ")
                #print(sol.name)
                for param in sol.parameters:
                    #print(f"{sol.parameters[param]['value']},    ")
                    all_values.write(f"{sol.parameters[param]['value']},    ")
                all_values.write('\n')
                all_value_count += 1
            all_values.close()
            opt.record_all_param(self.population,self.generation, flag=False)

            #self.population.children = pool.map(self.crud.evaluator,self.population.children)
            self.cleanup()
            self.population = self.selection.perform(self.population)
            opt.check_best_worst_average(self.population.parents)
            opt.write_track_file(self.population, self.generation)
            opt.record_optimized_solutions(self.population)

        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()
        all_values.close()
        opt.plotter()

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
            foo = self.generate_initial_solutions(f'initial_parent_{i}')
            foo.evaluate()
            self.population.parents.append(foo)

        for i in range(self.population.size):
            foo = self.generate_initial_solutions(f'initial_child_{i}')
            foo.evaluate()
            self.population.children.append(foo)

        self.population = self.selection.perform(self.population)
        opt.check_best_worst_average(self.population.parents)
        opt.write_track_file(self.population, self.generation)

        #scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        #uniqueness = Unique_Solution_Analyzer(scrambler)
        for self.generation.current in range(self.generation.total):
            self.population.children = self.reproduction.reproduce(self.population.parents, 
                                                                     self.solution)
            #self.population.children = uniqueness.analyze(self.population.children)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                #solution.add_additional_information(self.file_settings) #!TODO: redundant, remove.
                solution.evaluate()

            self.population = self.selection.perform(self.population)
            opt.check_best_worst_average(self.population.parents)
            opt.write_track_file(self.population, self.generation)

        opt.record_optimized_solutions(self.population)

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

    def evaluate_neural_networks(self,eval_pop):
        """
        So, the original method of analyzing neural networks in parallel didn't work.
        There were way too many memory failures. So I am now pursuing this course, of
        single evaluation. It will likely take a long time, unfortunately, but at present 
        it represents the best solution to the memory error problems. It will be significantly 
        sped up hopefully by performing the power analysis and interpolation in parallel.

        WRitten by Brian Andersen. 11/11/20
        """
        print("Performing neural network evaluations")
        for solution in eval_pop:
            print(solution.name)
            time1 = time.time()
            crud_matrix = None      #
            boron_matrix = None
            no_crud_matrix = True   #Flag saying I haven't done the CRUD predictions yet
            no_boron_matrix = True  #Flag saying I haven't done the Boron predictions yet.
            axially_adjusted_crud = False #Flag saying I haven't axially adjusted CRUD yet.
            if 'max_assembly_boron' in solution.parameters:
                if 'value' in solution.parameters['max_assembly_boron']:  #Don't do anything since it already failed inspection
                    pass
                else:
                    if os.path.isfile(f"{solution.name}/crud_input.h5"): #check input file exists. Maybe some sort of fluke
                        if no_boron_matrix:
                            if no_crud_matrix: #If the crud matrix hasn't been calculated yet.
                                crud_matrix = self.crud.predict_crud_distribution(solution.name)
                                no_crud_matrix = False
                            boron_matrix = self.crud.predict_boron_distribution(crud_matrix,solution.name)
                            boron_matrix = self.crud.adjust_prediction_axially(boron_matrix)
                            no_boron_matrix = False
                        
                        #crud_matrix = self.crud.adjust_prediction_axially(crud_matrix)
                        solution.parameters['max_assembly_boron']['value'] = self.crud.calculate_max_assembly_sum(boron_matrix)
                    else:
                        solution.parameters['max_assembly_boron']['value'] = 10000000. #If the input file doesn't exist give failed value
            if 'max_core_boron' in solution.parameters:
                if 'value' in solution.parameters['max_core_boron']:  #Don't do anything since it already failed inspection
                    pass
                else:
                    if os.path.isfile(f"{solution.name}/crud_input.h5"): #check input file exists. Maybe some sort of fluke
                        if no_boron_matrix:
                            if no_crud_matrix: #If the crud matrix hasn't been calculated yet.
                                crud_matrix = self.crud.predict_crud_distribution(solution.name)
                                no_crud_matrix = False
                            boron_matrix = self.crud.predict_boron_distribution(crud_matrix,solution.name)
                            boron_matrix = self.crud.adjust_prediction_axially(boron_matrix)
                            no_boron_matrix = False
                        
                        #crud_matrix = self.crud.adjust_prediction_axially(crud_matrix)
                        solution.parameters['max_core_boron']['value'] = self.crud.calculate_core_max_sum(boron_matrix)
                    else:
                        solution.parameters['max_core_boron']['value'] = 10000000. #If the input file doesn't exist give failed value
            if 'max_core_crud' in solution.parameters:
                if 'value' in solution.parameters['max_core_crud']:  #Don't do anything since it already failed inspection
                    pass
                else:
                    if os.path.isfile(f"{solution.name}/crud_input.h5"): #check input file exists. Maybe some sort of fluke
                        if no_crud_matrix: #If the crud matrix hasn't been calculated yet.
                            crud_matrix = self.crud.predict_crud_distribution(solution.name)
                            no_crud_matrix = False
                        if not axially_adjusted_crud:
                            crud_matrix = self.crud.adjust_prediction_axially(crud_matrix)
                            axially_adjusted_crud = True
                        solution.parameters['max_core_crud']['value'] = self.crud.calculate_core_max_sum(crud_matrix)
                    else:
                        solution.parameters['max_core_crud']['value'] = 10000000. #If the input file doesn't exist give failed value
            
            if 'max_assembly_crud' in solution.parameters:
                if 'value' in solution.parameters['max_core_crud']:  #Don't do anything since it already failed inspection
                    pass
                else:
                    if os.path.isfile(f"{solution.name}/crud_input.h5"): #check input file exists. Maybe some sort of fluke
                        if no_crud_matrix: #If the crud matrix hasn't been calculated yet.
                            crud_matrix = self.crud.predict_crud_distribution(solution.name)
                            no_crud_matrix = False
                        if not axially_adjusted_crud:
                            crud_matrix = self.crud.adjust_prediction_axially(crud_matrix)
                            axially_adjusted_crud = True
                        solution.parameters['max_assembly_crud']['value'] = self.crud.calculate_max_assembly_sum(crud_matrix)
                    else:
                        solution.parameters['max_assembly_crud']['value'] = 10000000. #If the input file doesn't exist give failed value
            print(time.time()-time1)
        return eval_pop

    def test_main_in_parallel(self):
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
            foo = self.solution()
            foo.name = "initial_parent_{}".format(i)
            foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
            #foo.add_additional_information(self.file_settings) #!TODO: redundant, remove.
            #scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
            if foo.fixed_genome:
                foo.new_generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                           self.file_settings['optimization']['fixed_groups'])
            #   foo.genome = scrambler.reproduce(foo.genome)
            else:
                foo.generate_initial(self.file_settings['genome']['chromosomes'])
            self.population.parents.append(foo)

        for i in range(self.population.size):
            foo = self.solution()
            foo.name = "initial_child_{}".format(i)
            foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
            #foo.add_additional_information(self.file_settings) #!TODO: redundant, remove.
            if foo.fixed_genome:
                foo.new_generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                           self.file_settings['optimization']['fixed_groups'])
            #    foo.genome = scrambler.reproduce(foo.genome)
            else:
                foo.generate_initial(self.file_settings['genome']['chromosomes'])
            self.population.children.append(foo)

        pool = Pool(processes=self.num_procs)
        self.population.parents = pool.map(test_evaluate_function, self.population.parents)
        self.population.children = pool.map(test_evaluate_function, self.population.children)
        all_values = open('all_value_tracker.txt','a')
        for sol in self.population.parents:
            #print(sol.name)
            if not all_value_count:
                for param in sol.parameters:
                    all_values.write(f"{param},    ")
                all_values.write('\n')
            all_values.write(f"{all_value_count},    {sol.genome}    ")
            for param in sol.parameters:
                all_values.write(f"{sol.parameters[param]['value']},    ")
                print(f"{sol.parameters[param]['value']},    ")
            all_values.write('\n')
            all_value_count += 1
        for sol in self.population.children:
            all_values.write(f"{all_value_count},   {sol.genome}   ")
            for param in sol.parameters:
                all_values.write(f"{sol.parameters[param]['value']},    ")
            all_values.write('\n')
            all_value_count += 1        
        all_values.close()


        #self.population.parents = pool.map(self.crud.evaluator,self.population.parents)
        #self.population.children = pool.map(self.crud.evaluator,self.population.children)
        self.population = self.selection.perform(self.population)
        opt.check_best_worst_average(self.population.parents)
        opt.write_track_file(self.population, self.generation)

        scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        uniqueness = Unique_Solution_Analyzer(scrambler)
        for self.generation.current in range(self.generation.total):
            self.population.children = self.reproduction.reproduce(self.population.parents, 
                                                                     self.solution)
            self.population.children = uniqueness.analyze(self.population.children)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                #solution.add_additional_information(self.file_settings) #!TODO: redundant, remove.


            self.population.children = pool.map(test_evaluate_function, self.population.children)
            all_values = open('all_value_tracker.txt','a')
            for sol in self.population.children:
                all_values.write(f"{all_value_count},    {sol.genome}    ")
                print(sol.name)
                for param in sol.parameters:
                    print(f"{sol.parameters[param]['value']},    ")
                    all_values.write(f"{sol.parameters[param]['value']},    ")
                all_values.write('\n')
                all_value_count += 1
            all_values.close()
            #self.population.children = pool.map(self.crud.evaluator,self.population.children)
            self.cleanup()
            self.population = self.selection.perform(self.population)
            opt.check_best_worst_average(self.population.parents)
            opt.write_track_file(self.population, self.generation)

        opt.record_optimized_solutions(self.population)

        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()
        all_values.close()
        opt.plotter()

    def test_main_in_serial(self):
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
            foo = self.solution()
            foo.name = "initial_parent_{}".format(i)
            foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
            #foo.add_additional_information(self.file_settings) #!TODO: redundant, remove.
            scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
            if foo.fixed_genome:
                foo.generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                           self.file_settings['optimization']['fixed_groups'])
                foo.genome = scrambler.reproduce(foo.genome)
            else:
                foo.generate_initial(self.file_settings['genome']['chromosomes'])
            test_evaluate_function(foo)
            self.population.parents.append(foo)

        for i in range(self.population.size):
            foo = self.solution()
            foo.name = "initial_child_{}".format(i)
            foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
            #foo.add_additional_information(self.file_settings) #!TODO: redundant, remove.
            if foo.fixed_genome:
                foo.generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                           self.file_settings['optimization']['fixed_groups'])
                foo.genome = scrambler.reproduce(foo.genome)
            else:
                foo.generate_initial(self.file_settings['genome']['chromosomes'])
            test_evaluate_function(foo)
            self.population.children.append(foo)

      #  self.population.parents = map(evaluate_function, self.population.parents)
      #  self.population.children = map(evaluate_function, self.population.children)
        self.population = self.selection.perform(self.population)
        opt.check_best_worst_average(self.population.parents)
        opt.write_track_file(self.population, self.generation)

        uniqueness = Unique_Solution_Analyzer(scrambler)
        for self.generation.current in range(self.generation.total):
            self.population.children = self.reproduction.reproduce(self.population.parents, 
                                                                     self.solution)
            self.population.children = uniqueness.analyze(self.population.children)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                #solution.add_additional_information(self.file_settings) #!TODO: redundant, remove.
                test_evaluate_function(solution)

           # self.population.children = map(evaluate_function, self.population.children)

            self.population = self.selection.perform(self.population)
            opt.check_best_worst_average(self.population.parents)
            opt.write_track_file(self.population, self.generation)

        opt.record_optimized_solutions(self.population)

        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()

        opt.plotter()

    def restart_main_in_parallel(self):
        """
        Performs optimization using a genetic algorithm in with
        parallel computations

        Parameters: None

        Written by Brian Andersen. 1/9/2020
        """
        opt = Optimization_Metric_Toolbox()

        track_file = open('optimization_track_file.txt', 'a')
        track_file.write("Restarting Optimization \n")
        track_file.close()

        all_value_count = 0
        all_values = open('all_value_tracker.txt','a')
        all_values.write("Restarting\n")
        all_values.close()
        with open('optimized_solutions.yaml') as current_solutions:
            solutions = yaml.safe_load(current_solutions)

        track_file = open('optimization_track_file.txt','r')
        track_lines = track_file.readlines()
        track_file.close()

        for line in track_lines:
            if "Current generation of the optimization" in line:
                elems = line.strip().split()
                current_gen = int(elems[-1]) + 1
        track_lines = None

        for solut in solutions:
            foo = self.solution()
            foo.name = solut
            foo.genome = solutions[solut]['genome']
            foo.parameters = solutions[solut]['parameters']
            foo.fitness = float(solutions[solut]['fitness'])
            self.population.parents.append(foo)

        pool = Pool(processes=self.num_procs)

        #self.population.parents = pool.map(self.crud.evaluator,self.population.parents)
        #self.population.children = pool.map(self.crud.evaluator,self.population.children)
        
        scrambler = Fixed_Genome_Mutator(1,1,200,self.file_settings)
        uniqueness = Unique_Solution_Analyzer(scrambler)
        for self.generation.current in range(current_gen,self.generation.total):
            self.population.children = self.reproduction.reproduce(self.population.parents, 
                                                                     self.solution)
            self.population.children = uniqueness.analyze(self.population.children)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                #solution.add_additional_information(self.file_settings) #!TODO: redundant, remove.


            self.population.children = pool.map(evaluate_function, self.population.children)
            all_values = open('all_value_tracker.txt','a')
            for sol in self.population.children:
                all_values.write(f"{all_value_count},   {sol.name},   ")
                print(sol.name)
                for param in sol.parameters:
                    print(f"{sol.parameters[param]['value']},    ")
                    all_values.write(f"{sol.parameters[param]['value']},    ")
                all_values.write('\n')
                all_value_count += 1
            all_values.close()
            #self.population.children = pool.map(self.crud.evaluator,self.population.children)
            self.cleanup()
            self.population = self.selection.perform(self.population)
            opt.check_best_worst_average(self.population.parents)
            opt.write_track_file(self.population, self.generation)
            opt.record_optimized_solutions(self.population)

        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()
        all_values.close()
        opt.plotter()


class Population(object):
    """
    Population Size of the genetic algorithm.
    """
    def __init__(self, population=None):
        self.size = population
        self.parents = []
        self.children = []
        self.solution_front = None

    def calculate_size(self, number_genes):
        """
        Calculates the population size based on total number genes.
        """
        self.size = int(10*math.sqrt(number_genes))


class Generation(object):
    """
    Number of generations over the course of the optimization.
    """
    def __init__(self, total=None):
        self.total = total
        self.current = 0

    def calculate_total_generations(self, number_changes):
        """
        Calculates the total number of generations in an optimization
        based on the number of changes.
        """
        self.total = int(5*math.sqrt(number_changes))


class GA_Selection(object):
    """Class to determine which solutions make it onto the solution front and
    which solutions become the parents of the next generation.

    Parameters
    ----------
    fitness: class
        Function for calculating the fitness of solutions in the optimization
    solution_front: class
        solutions on the solution front of the genetic algorithm.
    method: function
        Function for determining which solutions will be used to determine
        the next generation of parents.

    Written by Brian Andersen. 1/9/2020
    """
    def __init__(self, fitness=None,
                 method=None,
                 solution_front=None):

        self.fitness = fitness
        self.method = method
        self.solution_front = solution_front

    @staticmethod
    def tournament(solution_list, desired_number_solutions):
        """Performs a tournament selection of the solutions. Can be used for determining
        solution front or parents for the next generation or whatever.

        If the desired number of tournament winners is greater than the desired
        number of solutions, a second tournament is held to select the additional solutions
        needed. If the desired number of solutions is less than the number of solutions in
        solution list, else half the solutions are returned.
        Parameters
        -----------
        solution_list: list
            The solutions that will compete in the tournament
        desired_number_solutions: int
            The number of solutions chosen at the end of the tournament.

        Written by Brian Andersen. 1/9/2020
        """
        unused_solutions = solution_list[:]
        used_solutions = []
        winners = []
        for i in range(desired_number_solutions):
            one = random.choice(unused_solutions)
            two = random.choice(unused_solutions)
            while one == two:
                two = random.choice(unused_solutions)
            if one.fitness > two.fitness:
                winners.append(one)
                used_solutions.append(two)
            else:
                winners.append(two)
                used_solutions.append(one)

            if len(unused_solutions) < 2:
                if len(used_solutions) > 2:
                    unused_solutions = used_solutions
                    used_solutions = []
                else:
                    unused_solutions = copy.deepcopy(solution_list)
                    used_solutions = []

            unused_solutions.remove(one)
            unused_solutions.remove(two)

        return winners

    @staticmethod
    def roulette(solution_list, desired_number_solutions):
        """Performs a roulette selection of the solutions. Can be used for determining
        solution front or parents for the next generation or whatever."""
        unused_solutions = copy.deepcopy(solution_list)
        winners = []
        for i in range(desired_number_solutions):
            probability_sum = 0
            selection_probability = {}
            selection_probability['low_bound'] = []
            selection_probability['up_bound']  = []
            for solution in unused_solutions:
                selection_probability['low_bound'].append(probability_sum)
                probability_sum += solution.fitness
                selection_probability['up_bound'].append(probability_sum)

            value = random.random()
            value = value*probability_sum
            for j, solution in enumerate(unused_solutions):
                if(selection_probability['low_bound'][j] <= value
                   and value <= selection_probability['up_bound'][j]):
                    winners.append(solution)
                    unused_solutions.remove(solution)

            if not unused_solutions:
                unused_solutions = copy.deepcopy(solution_list)

        return winners

    def perform(self, population_class):
        solution_list = copy.deepcopy(population_class.parents)
        solution_list.extend(population_class.children)

        if not population_class.solution_front:
            pass
        else:
            solution_list.extend(population_class.solution_front)
        solution_list = self.fitness.calculate(solution_list)
        population_class.solution_front = self.solution_front.calculate(solution_list)
        if not population_class.solution_front:
            population_class.parents = self.method(solution_list, population_class.size)
        else:
            population_class.parents = self.method(population_class.solution_front,
                                                   population_class.size)
        population_class.children = []

        return population_class


class MOOGLE(GA_Selection):
    """
    MOOGLE selection and fitness class.

    MOOGLE uses binned solution space, and bases fitness off of
    binning as well.

    Parameters
    -------------
    bin_ : see bin_list.

    Variables:
    -------------
    solution_front: dict
        Dictionary of solutions that have made it onto the binned solution space used in the
        MOOGLE optimization.
    bin_list : dict
        Organized dictionary taken from the input yaml file.
        Dictionary is organized in the form:
            self.bin[parameter]['minimum']
            self.bin[parameter]['maximum']
            self.bin[parameter]['bin_size']
            self.bin[parameter]['objective']

    Written by Brian Andersen. 1/9/2020
    """
    def __init__(self):
        GA_Selection.__init__(self)
        temp = fitness.Binned_Fitness()
        self.fitness = fitness.Fitness_Maxizer(temp)
        temp = fitness.Ranked_Fitness()
        self.no_bin_fitness = fitness.Fitness_Maxizer(temp)

    def perform(self, population_class):
        """
        Performs solution selection using the MOOGLE methodology.

        Returns the population class instance with an updated parent list and
        solution front list.

        Parameters
        -----------
        population_class: Class
            The solution class used for the optimization.
        """
        if not population_class.solution_front:
            solution_list = population_class.parents[:]
            solution_list.extend(population_class.children)
        else:
            solution_list = population_class.solution_front[:]
            solution_list.extend(population_class.children)
        solution_list = self.calculate_bin_values(solution_list)
        population_class.solution_front = self.bin_solutions(solution_list)
        if not population_class.solution_front:
            solution_list = self.no_bin_fitness.calculate(solution_list)
            population_class.parents = self.tournament(solution_list, population_class.size)
        else:
            solution_list = self.fitness.calculate(solution_list)
            population_class.parents = self.roulette(population_class.solution_front,
                                                     population_class.size)

        return population_class

    def calculate_bin_values(self, solution_list):
        """
        Calculates the appropriate bin for the optimization solution.

        If the solution is outside of the desired solution space
        the function returns the string "solution_not_binned".
        Otherwise it returns a tuple of the appropriate bin
        values.

        Parameters
        ------------
        solution : class
            Solution to the optimization value with parameters
            already determined.
        """
        for solution in solution_list:
            for param in solution.parameters:
                if 'value' in solution.parameters[param]:
                    pass
                else:
                    solution.raise_value_error()
                goal = solution.parameters[param]['goal'].lower()
                #bin_ = solution.parameters[param]['bin']
                value = solution.parameters[param]['value']
                if goal == 'maximize':
                    if value < solution.parameters[param]['maximum'] and value > solution.parameters[param]['minimum']:
                        bin_ = value
                        bin_ -= solution.parameters[param]['minimum']
                        bin_ /= solution.parameters[param]['bin_size']
                        solution.parameters[param]['bin'] = bin_
                    else:
                        solution.parameters[param]['bin'] = 'solution_not_binned'
                elif goal == 'minimize':
                    if value < solution.parameters[param]['maximum'] and value > solution.parameters[param]['minimum']:
                        bin_ = solution.parameters[param]['maximum']
                        bin_ -= value
                        bin_ /= solution.parameters[param]['bin_size']
                        solution.parameters[param]['bin'] = bin_
                    else:
                        solution.parameters[param]['bin'] = 'solution_not_binned'
                elif goal == 'meet_target':
                    bin_ = solution.parameters[param]['target']
                    bin_ -= value
                    bin_ /= solution.parameters[param]['bin_size']
                    solution.parameters[param]['bin'] = bin_
                elif goal == 'less_than_target':
                    if value < solution.parameters[param]['target']:
                        solution.parameters[param]['bin'] = 0
                    else:
                        bin_ = 'solution_not_binned'
                elif goal == 'greater_than_target':
                    if value > solution.parameters[param]['target']:
                        solution.parameters[param]['bin'] = 0
                    else:
                        solution.parameters[param]['bin'] = 'solution_not_binned'
                else:
                    raise NotImplementedError

        return solution_list

    def bin_solutions(self, solution_list):
        """
        Determines which solutions make it onto the binned solution front based on the calculated
        bin values.

        Solution space binning is performed using a dictionary. Dictionary keys are formed by
        the bin numbers calculated. If multiple solutions occupy the same bin, the first solution
        to be placed into the bin is kept, and all other solutions are discarded.

        Once binning has occured the dictionary is converted to a list and returned.
        """
        bin_dictionary = {}
        for solution in solution_list:
            solution_binned = True
            temp = []
            for param in solution.parameters:
                if solution.parameters[param]['bin'] == 'solution_not_binned':
                    solution_binned = False
                    break
                else:
                    temp.append(solution.parameters[param]['bin'])
            if solution_binned:
                key = tuple(temp)
                if key not in bin_dictionary:
                    bin_dictionary[key] = solution

        bin_list = []
        for key in bin_dictionary:
            bin_list.append(bin_dictionary[key])

        if bin_list:
            return bin_list
        else:
            return None


class Reproduction(object):
    """
    Class for choosing which solutions undergo crossover and mutation, as well
    as performing crossover. Class has been updated as of 1/20/2020 to allow for
    fixed genome groups, e.g. the three assembly enrichments that must be
    maintained for the NE 412/512 class project. The entire class was modified,
    rather than creating a new class that inherits because the main crossover
    function is only changed by three lines, and the other added functions may
    prove useful for future classes to inherit these functions as well.

    Parameters
    -----------
    mutation: class
        The method in which solutions are mutated into new solutions in the
        genetic algorithm.

    Written by Brian Andersen. 1/9/2020
    """
    def __init__(self, mutation):
        self.mutation = mutation
        self.mutation_list = None
        self.crossover_list = None

    def select_reproduction_method(self, solution_list):
        """
        Function to determine which assemblies undergo crossover and which assemblies
        undergo mutation.
        """
        self.mutation_list  = []
        self.crossover_list = []
        for solution in solution_list:

            if random.random() < self.mutation.rate:
                self.mutation_list.append(solution.genome)
            else:
                self.crossover_list.append(solution.genome)

        if len(self.crossover_list) %2 == 1:
            swap_solution = random.choice(self.mutation_list)
            self.crossover_list.append(swap_solution)

    def mate_crossover_solutions(self):
        """
        Function for mating the solutions selected for crossover.
        """
        first_mate_list  = []
        second_mate_list = []
        while len(self.crossover_list) > 0:
            parent_one = random.choice(self.crossover_list)

            partner_similarities = -1
            for suitor in self.crossover_list:
                if suitor is parent_one:
                    pass
                else:
                    new_suitor_similarities = self.check_intersection(parent_one, suitor)
                    if new_suitor_similarities > partner_similarities:
                        parent_two = suitor[:]
                        partner_similarities = new_suitor_similarities

            self.crossover_list.remove(parent_one)
            if parent_two in self.crossover_list:
                self.crossover_list.remove(parent_two)

            first_mate_list.append(parent_one)
            second_mate_list.append(parent_two)
        return first_mate_list, second_mate_list

    def crossover(self, genome_one, genome_two, mutation_rate):
        """
        Function for performing crossover of the mated solutions by swapping
        differing genes between the genomes. Edited as of 1/20/2020 to account
        for fixed gene counts in problems.

        Parameters:
            genome_one: list
                The first genome that is to undergo crossover.
            genome_two: list
                The second geneome that is to undergo crossover.
            mutation_rate: float
                The percent of solutions that undergo mutation. Used here to
                determine the number of genes that are swapped.

        Written by Brian Andersen. 8/29/2019
        """
        difference_positions = self.return_different_positions(genome_one,
                                                               genome_two)
        child_one = []
        child_two = []
        position_count = 0
        for i, j in zip(genome_one, genome_two):
            if position_count in difference_positions:
                if random.random() < mutation_rate:
                    child_one.append(j)
                    child_two.append(i)
                else:
                    child_one.append(i)
                    child_two.append(j)
            else:
                child_one.append(i)
                child_two.append(j)
            position_count += 1

        return child_one, child_two

    def reproduce(self, inp_obj, solution_list, solution_class):
        """
        Performs all functions related to reproduction of the solution.
        """
        self.input = inp_obj
        self.select_reproduction_method(solution_list)

        first_mate_list, second_mate_list = self.mate_crossover_solutions()

        child_genome_list = []
        for mate_one, mate_two in zip(first_mate_list, second_mate_list):
            child_one, child_two = self.crossover(mate_one, mate_two, self.mutation.rate)
            child_genome_list.extend([child_one, child_two])

        for genome in self.mutation_list:
            child = self.mutation.reproduce(genome)
            child_genome_list.append(child)

        new_solution_list = []
        for child in child_genome_list:

            foo = solution_class(self.input)#!TODO:this needs to be narrowed to just self.input
            foo.genome = child

            new_solution_list.append(foo)

        return new_solution_list

    @staticmethod
    def check_intersection(list_one, list_two):
        """
        Counts the number of matches between two lists.
        Written by Brian Andersen. 11/24/2019
        """
        match_count = 0
        for i, j in zip(list_one, list_two):
            if i == j:
                match_count += 1

        return match_count

    @staticmethod
    def return_different_positions(list_one, list_two):
        """
        Returns the positions in the list where the two lists do not match.
        Written by Brian Andersen. 11/24/2019
        """
        diff_position_list = []
        position_count = 0
        for i, j in zip(list_one, list_two):
            if i == j:
                pass
            else:
                diff_position_list.append(position_count)
            position_count += 1

        return diff_position_list


class Dictionary_Genome_Reproducer(Reproduction):
    """
    Reproduction class used when genomes are designed by a dictionary, rather
    than the traditional list, this reproduction class is implemented. The
    The key difference is that the entire solution is passed in this class,
    not just the solution genome.

    Written by Brian Andersen, 11/15/2019.
    """
    def __init__(self, mutator_):
        Reproduction.__init__(self, mutator_)

    def select_reproduction_method(self, solution_list):
        """
        Class for choosing which solutions undergo crossover and mutation.
        """
        self.mutation_list  = []
        self.crossover_list = []
        for solution in solution_list:

            if random.random() < self.mutation.rate:
                self.mutation_list.append(solution)
            else:
                self.crossover_list.append(solution)

        if len(self.crossover_list) %2 == 1:
            swap_solution = random.choice(self.mutation_list)
            self.crossover_list.append(swap_solution)

    def mate_crossover_solutions(self):
        """
        Function for mating the solutions selected for crossover.
        """
        first_mate_list  = []
        second_mate_list = []
        while len(self.crossover_list) > 0:
            parent_one = random.choice(self.crossover_list)

            partner_similarities = 0
            for suitor in self.crossover_list:
                if suitor is parent_one:
                    pass
                else:
                    new_suitor_similarities = 0
                    for key in parent_one.genome:
                        new_suitor_similarities += self.check_intersection(parent_one.genome[key],
                                                                           suitor.genome[key])
                        if new_suitor_similarities > partner_similarities:
                            parent_two = suitor
                            partner_similarities = new_suitor_similarities

            self.crossover_list.remove(parent_one)

            self.crossover_list.remove(parent_two)

            first_mate_list.append(parent_one)
            second_mate_list.append(parent_two)
        return first_mate_list, second_mate_list

    def reproduce(self, solution_list, solution_class):
        """
        Performs all functions related to reproduction of the solution.
        """
        self.select_reproduction_method(solution_list)
        first_mate_list, second_mate_list = self.mate_crossover_solutions()
        child_genome_list = []
        for mate_one, mate_two in zip(first_mate_list, second_mate_list):
            child_one = {}
            child_two = {}
            for par in mate_one.genome:
                child_one[par], child_two[par] = self.crossover(mate_one.genome[par],
                                                                mate_two.genome[par],
                                                                self.mutation.rate)
            child_genome_list.extend([child_one, child_two])

        for solution in self.mutation_list:
            child = self.mutation.reproduce(solution)
            child_genome_list.append(child)

        new_solution_list = []
        for child in child_genome_list:
            foo = solution_class()
            foo.genome = child
            new_solution_list.append(foo)

        return new_solution_list


class Fixed_Gene_Reproducer(Reproduction):
    """
    Gene for reproducing when the genes intended to be used are held fixed.
    I.E. Every solution is going to be made from the same combination of
    genes.
    """
    def __init__(self, mutator_):
        Reproduction.__init__(self, mutator_)

    def crossover(self, genome_one, genome_two, mutation_rate):
        """
        Function for performing crossover of the mated solutions
        when the genes are held fixed.
        Written by Brian Andersen 12/4/2019
        """
        difference_positions = self.return_different_positions(genome_one,
                                                               genome_two)
        if len(difference_positions) <= 4:
            child_one = self.mutation.reproduce(genome_one)
            child_two = self.mutation.reproduce(genome_two)
        else:
            swap_list = random.sample(difference_positions,
                                      int(len(difference_positions)/2))
            match_list = []
            while len(swap_list) > 1:
                one = swap_list.pop(0)
                for two in swap_list[1:]:
                    if one in self.mutation.symmetry_list:
                        if two in self.mutation.symmetry_list:
                            if one == two:
                                pass
                            else:
                                is1ok = self.mutation.genome_map[genome_one[two]][one]
                                is2ok = self.mutation.genome_map[genome_one[one]][two]
                                is3ok = self.mutation.genome_map[genome_two[two]][one]
                                is4ok = self.mutation.genome_map[genome_two[one]][two]
                                if is1ok == 1 and is2ok == 1:
                                    if is3ok == 1 and is4ok == 1:
                                        match_list.append([one,two])
                                        swap_list.remove(two)
                                        break
                    else:
                        if two in self.mutation.symmetry_list:
                            pass
                        else:
                            if one == two:
                                pass
                            else:
                                is1ok = self.mutation.genome_map[genome_one[two]][one]
                                is2ok = self.mutation.genome_map[genome_one[one]][two]
                                is3ok = self.mutation.genome_map[genome_two[two]][one]
                                is4ok = self.mutation.genome_map[genome_two[one]][two]
                                if is1ok == 1 and is2ok == 1:
                                    if is3ok == 1 and is4ok == 1:
                                        match_list.append([one,two])
                                        swap_list.remove(two)
                                        break
            child_one = copy.deepcopy(genome_one)

            for match in match_list:
                temp_one = child_one[match[0]]
                temp_two = child_one[match[1]]
                child_one[match[0]] = temp_two
                child_one[match[1]] = temp_one

            child_two = copy.deepcopy(genome_two)
            for match in match_list:
                temp_one = child_two[match[0]]
                temp_two = child_two[match[1]]
                child_two[match[0]] = temp_two
                child_two[match[1]] = temp_one

        return child_one, child_two


class Unique_Genome_Reproducer(Reproduction):
    """
    Class for producing crossover when Unique Genes are present in the genome, but 
    pure genome fixity isn't necessary.
    """
    def __init__(self,mutator_,settings):
        Reproduction.__init__(self,mutator_)

        self.common_chromosomes = settings['optimization']['mutation']['common_chromosomes']
        self.fixed_groups = settings['optimization']['fixed_groups']

    def return_operation_positions(self,genome_one,genome_two,genome_key):
        """
        Returns the positions of genes that can be swapped directly via crossover and those
        that will be switched.

        Written by Brian Andersen. 11/20/2020.
        """
        crossover_position_list = []
        swap_position_list = []
        position_count = 0 #Position count is used because enumerate becomes complicated
                           #while using zip.
        for gene_one,gene_two in zip(genome_one,genome_two): #For every gene in both genomes
            if gene_one == gene_two:                         #Do nothing if they match
                pass
            else:                                            #Otherwise
                perform_swap = False
                if 'unique' in genome_key[gene_one]: #If gene one is unique and used in genome
                    if genome_key[gene_one]['unique']: #two the position will be marked for 
                        if gene_one in genome_two:     #swap
                            perform_swap = True
                if 'unique' in genome_key[gene_two]: #Same as above only for gene two and 
                    if genome_key[gene_two]['unique']: #genome one.
                        if gene_two in genome_one:
                            perform_swap = True
                if perform_swap:                              #If designated for swap, add to            
                    swap_position_list.append(position_count) #swap position list
                else:                                         #otherwise
                    for group in self.common_chromosomes:
                        if gene_one in self.common_chromosomes[group]: #determine which common
                            group_one = group                          #gene list each genome is
                        if gene_two in self.common_chromosomes[group]: #in. If they are in the
                            group_two = group                          #same group, perform 
                    if group_one == group_two:                         #crossover, otherwise add
                        crossover_position_list.append(position_count) #to swap list. Done to 
                    else:                                              #preserve number of genes
                        swap_position_list.append(position_count)      #in each gene group.

            position_count += 1

        return crossover_position_list,swap_position_list

    def gene_can_swap(self,gene_one,gene_two,genome_key):
        """
        Returns True if the gene can be swapped between the genome.
        Returns True if both genes are not unique and in the same common
        gene group.
        """
        gene_can_swap = True
        if 'unique' in genome_key[gene_one]:
            if genome_key[gene_one]['unique']:
                gene_can_swap = False
        if 'unique' in genome_key[gene_two]:
            if genome_key[gene_two]['unique']:
                gene_can_swap = False
        for group in self.common_chromosomes:
            if gene_one in self.common_chromosomes[group]:
                if gene_two in self.common_chromosomes[group]:
                    pass
                else:
                    gene_can_swap = False

        return gene_can_swap

    def crossover(self,genome_one,genome_two,mutation_rate):
        """
        Function for performing crossover of the mated solutions when unique genes are 
        present in the optimization. Crossover is performed normally, however, if one of 
        the genes being swapped is a unique gene already present in the genome, these positions
        are swapped.
        """
        file_ = open("genome_key",'rb')
        genome_key = pickle.load(file_)
        file_.close()

        
        #crossover_rate = 0.5 #I want roughly 50% of all differing genomes to 
        #                     #undergo crossover. This seems cleaner than random lists.


        difference_positions = self.return_different_positions(genome_one,
                                                               genome_two)
        if len(difference_positions) <= 4:
            for i in range(8):
                child_one = self.mutation.reproduce(genome_one)
                child_two = self.mutation.reproduce(genome_two)
            return child_one,child_two
        
        first_child_one = []
        first_child_two = []
        position_count = 0
        for one,two in zip(genome_one,genome_two):
            if position_count in difference_positions:
                if self.gene_can_swap(one,two,genome_key):
                    first_child_one.append(two)
                    first_child_two.append(one)
                    difference_positions.remove(position_count)
                else:
                    first_child_one.append(one)
                    first_child_two.append(two)
            else:
                first_child_one.append(one)
                first_child_two.append(two)
            position_count += 1
        
        swap_list = random.sample(difference_positions,
                                int(len(difference_positions)/2))
        match_list = []
        while len(swap_list) > 1:
            one = swap_list.pop(0)
            for two in swap_list[1:]:
                if one in self.mutation.symmetry_list:
                    if two in self.mutation.symmetry_list:
                        if one == two:
                            pass
                        else:
                            is1ok = self.mutation.genome_map[first_child_one[two]][one]
                            is2ok = self.mutation.genome_map[first_child_one[one]][two]
                            is3ok = self.mutation.genome_map[first_child_two[two]][one]
                            is4ok = self.mutation.genome_map[first_child_two[one]][two]
                            if is1ok == 1 and is2ok == 1:
                                if is3ok == 1 and is4ok == 1:
                                    match_list.append([one,two])
                                    swap_list.remove(two)
                                    break
                else:
                    if two in self.mutation.symmetry_list:
                        pass
                    else:
                        if one == two:
                            pass
                        else:
                            is1ok = self.mutation.genome_map[first_child_one[two]][one]
                            is2ok = self.mutation.genome_map[first_child_one[one]][two]
                            is3ok = self.mutation.genome_map[first_child_two[two]][one]
                            is4ok = self.mutation.genome_map[first_child_two[one]][two]
                            if is1ok == 1 and is2ok == 1:
                                if is3ok == 1 and is4ok == 1:
                                    match_list.append([one,two])
                                    swap_list.remove(two)
                                    break
        child_one = copy.deepcopy(first_child_one)
        for match in match_list:
            temp_one = child_one[match[0]]
            temp_two = child_one[match[1]]
            child_one[match[0]] = temp_two
            child_one[match[1]] = temp_one

        child_two = copy.deepcopy(first_child_two)
        for match in match_list:
            temp_one = child_two[match[0]]
            temp_two = child_two[match[1]]
            child_two[match[0]] = temp_two
            child_two[match[1]] = temp_one

        return child_one, child_two

    def first_crossover(self,genome_one,genome_two,mutation_rate):
        """
        Function for performing crossover of the mated solutions when unique genes are 
        present in the optimization. Crossover is performed normally, however, if one of 
        the genes being swapped is a unique gene already present in the genome, these positions
        are swapped.
        """
        file_ = open("genome_key",'rb')
        genome_key = pickle.load(file_)
        file_.close()

        crossover_rate = 0.5 #I want roughly 50% of all differing genomes to 
                             #undergo crossover. This seems cleaner than random lists.
            
        difference_positions = self.return_different_positions(genome_one,
                                                               genome_two)
        if len(difference_positions) < 4:
            child_one = self.mutation.reproduce(genome_one)
            child_two = self.mutation.reproduce(genome_two)
        else:
            position_count = 0
            child_one = [None]*len(genome_one)
            child_two = [None]*len(genome_one)
            for i,j in zip(genome_one,genome_two):
                if position_count in difference_positions:
                    if random.random() < crossover_rate:
                        if not child_one[position_count]:
                            if 'unique' in genome_key[i]:
                                if genome_key[i]['unique']:
                                    if i in genome_two:
                                        other_pos = genome_two.index(i)
                                        child_two[other_pos] = genome_two[position_count]
                                        child_two[position_count] = i
                                    else:
                                        child_two[position_count] = i
                                else:
                                        child_two[position_count] = i
                            else:
                                        child_two[position_count] = i
                        if not child_two[position_count]:
                            if 'unique' in genome_key[j]:
                                if genome_key[j]['unique']:
                                    if j in genome_one:
                                        other_pos = genome_one.index(j)
                                        child_one[other_pos] = genome_one[position_count]
                                        child_one[position_count] = j
                                    else:
                                        child_one[position_count] = j
                                else:
                                        child_one[position_count] = j
                            else:
                                        child_one[position_count] = j 
                    else:
                        if not child_one[position_count]:
                            child_one[position_count] = i
                        if not child_two[position_count]:
                            child_two[position_count] = j
                else:
                    if not child_one[position_count]:
                        child_one[position_count] = i
                    if not child_two[position_count]:
                        child_two[position_count] = j
                position_count += 1

        return child_one,child_two    


class Mcycle_Reproducer(Reproduction):
    """
    Class for reproducing genomes relevant to multicycle problems in the presence of both
    fresh and burned fuel.
    Written by G. K. Delipei 10/17/2023
    """
    def __init__(self, mutator_, settings):
        super().__init__(mutator_)
        self.unique_gene_list = []
        self.nonunique_gene_list = []
        self._fill_unique_list(settings)
        self.gene_list =  self.unique_gene_list + self.nonunique_gene_list
        self.ncycles = int(settings['genome']['parcs_data']['ncycles'])
    
    def _fill_unique_list(self,settings):
        """
        This is an idea for genes of a unique type, meaning that only one is allowed in the optimization solution.
        It is thought of as a way to include burnt fuel assemblies. Right now, each burnt fuel assembly is its own 
        group. But I was thinking what if I put them into a group. Then, in case you inclue a bad assembly, such as low 
        fuel enrichment or whatever, it could be easily selected out.
        """
        GENES = settings['genome']['chromosomes']
        for key in GENES:
            if 'unique' in GENES[key]:
                if GENES[key]['unique']:
                    self.unique_gene_list.append(key)
                else:
                    self.nonunique_gene_list.append(key)
            else:
                self.nonunique_gene_list.append(key)


    def crossover(self,genome_one,genome_two,mutation_rate):
        """
        Function for performing crossover of the mated solutions when unique genes are 
        present in the optimization. Crossover is performed normally, however, if one of 
        the genes being swapped is a unique gene already present in the genome, these positions
        are swapped.
        """
        file_ = open("genome_key",'rb')
        genome_key = pickle.load(file_)
        file_.close()

        set_all_genes = set(self.gene_list)
       
        child_one = []
        child_two = []
        for one,two in zip(genome_one,genome_two):
            if one == two:
                child_one.append(one)
                child_two.append(two)
            else:
                if one in self.nonunique_gene_list and two in self.nonunique_gene_list:
                    if random.random() < 0.5:
                        child_one.append(two)
                        child_two.append(one)
                    else:
                        child_one.append(one)
                        child_two.append(two)
                elif one in self.unique_gene_list and two in self.nonunique_gene_list:
                    if random.random() < 0.5:
                        child_one.append(two)
                        if one not in child_two:
                            child_two.append(one)
                        else:
                            non_unique_child_two =[]
                            non_unique_child_two.append(two)
                            for gene in child_two:
                                if gene in self.nonunique_gene_list:
                                    non_unique_child_two.append(gene)
                            set_avail_genes= set_all_genes-set(non_unique_child_two)
                            avail_genes = list(set_avail_genes)
                            one_new = random.choice(avail_genes)
                            child_two.append(one_new)
                    else:
                        child_one.append(one)
                        child_two.append(two)      
                elif one in self.nonunique_gene_list and two in self.unique_gene_list:
                    if random.random() < 0.5:
                        child_two.append(one)
                        if two not in child_one:
                            child_one.append(two)
                        else:
                            non_unique_child_one =[]
                            non_unique_child_one.append(one)
                            for gene in child_two:
                                if gene in self.nonunique_gene_list:
                                    non_unique_child_one.append(gene)
                            set_avail_genes= set_all_genes-set(non_unique_child_one)
                            avail_genes = list(set_avail_genes)
                            two_new = random.choice(avail_genes)
                            child_one.append(two_new)
                    else:
                        child_one.append(one)
                        child_two.append(two)  
                elif one in self.unique_gene_list and two in self.unique_gene_list:
                    if random.random() < 0.5:
                        if one not in child_two:
                            child_two.append(one)
                        else:
                            non_unique_child_two =[]
                            non_unique_child_two.append(two)
                            for gene in child_two:
                                if gene in self.nonunique_gene_list:
                                    non_unique_child_two.append(gene)
                            set_avail_genes= set_all_genes-set(non_unique_child_two)
                            avail_genes = list(set_avail_genes)
                            one_new = random.choice(avail_genes)
                            child_two.append(one_new)
                        if two not in child_one:
                            child_one.append(two)
                        else:
                            non_unique_child_one =[]
                            non_unique_child_one.append(one)
                            for gene in child_two:
                                if gene in self.nonunique_gene_list:
                                    non_unique_child_one.append(gene)
                            set_avail_genes= set_all_genes-set(non_unique_child_one)
                            avail_genes = list(set_avail_genes)
                            two_new = random.choice(avail_genes)
                            child_one.append(two_new)
                    else:
                        child_one.append(one)
                        child_two.append(two)  

        return child_one, child_two

    def reproduce(self, solution_list, solution_class):
        """
        Performs all functions related to reproduction of the solution.
        """
        ngen=len(solution_list[0].genome)
        cycle_genes = int(ngen/self.ncycles)
        start_cycle_id = []
        end_cycle_id = []
        for i in range(self.ncycles):
            istart = i*cycle_genes
            iend = i*cycle_genes + cycle_genes
            start_cycle_id.append(istart)
            end_cycle_id.append(iend)
    
        self.select_reproduction_method(solution_list)

        first_mate_list, second_mate_list = self.mate_crossover_solutions()

        child_genome_list = []
        for mate_one, mate_two in zip(first_mate_list, second_mate_list):
            child_one = []
            child_two = []
            for i in range(self.ncycles):
                ci_mate_one = mate_one[start_cycle_id[i]:end_cycle_id[i]]
                ci_mate_two = mate_two[start_cycle_id[i]:end_cycle_id[i]]
                ci_child_one, ci_child_two = self.crossover(ci_mate_one, ci_mate_two, self.mutation.rate)
                child_one = child_one + ci_child_one
                child_two = child_two + ci_child_two
            child_genome_list.extend([child_one, child_two])

        for genome in self.mutation_list:
            child = []
            for i in range(self.ncycles):
                cigenome = genome[start_cycle_id[i]:end_cycle_id[i]]
                cichild = self.mutation.reproduce(cigenome)
                child = child + cichild
            child_genome_list.append(child)

        new_solution_list = []
        for child in child_genome_list:

            foo = solution_class()
            foo.genome = child

            new_solution_list.append(foo)

        return new_solution_list


class Mutation(object):
    """
    Class for Genetic Algorithm Mutation.
    """
    def __init__(self, rate,
                 final_rate,
                 number,
                 settings
                 ):
        self.rate = rate
        self.number = number
        self.final_rate = final_rate
        if self.rate == self.final_rate:
            self.constant_mutation = True
        else:
            self.constant_mutation = False
        self.rate_incease = None
        self.genome_map = {}
        self.symmetry_list = []
        self._fill_genome_map(settings)

    def _fill_genome_map(self, settings):
        """
        Fills in the genome map used to determine which genes may be placed where
        when mutation is performed.

        Parameters:
            settings: Dictionary
                The input yaml file, as a dictionary. 

        Written by Brian Andersen. 1/21/2020.
        """
        if not settings:
            pass
        else:
            chrom_settings = settings['genome']['chromosomes']['fuel']
            for chrom in chrom_settings:
                if chrom == "symmetry_list":
                    self.symmetry_list = chrom_settings[chrom]
                else:
                    self.genome_map[chrom] = copy.deepcopy(chrom_settings[chrom]['map'])

    def calculate_rate_increase(self, number_generations):
        """
        Calculates the rate at which the mutation rate increases.
        """       
        if self.constant_mutation:
            pass
        else:
            temp = (1-self.final_rate)/(1-self.rate)
            temp = math.log(temp)
            temp /= number_generations
            self.rate_incease = math.exp(temp)

    def update_mutation_rate(self):
        """
        Updates the mutation rate of the genetic algorithm.
        """
        if self.constant_mutation:
            pass
        else:
            self.rate =  1-self.rate_incease*(1-self.rate)

    def reproduce(self, genome):
        """
        The function used to perform the specific mutation method.
        """
        raise NotImplementedError


class Mutate_By_Type(Mutation):
    """
    Performs mutation using common rod types. 
    Inherits from class Mutation.
    """
    def __init__(self, rate, final_rate, number, settings=None):
        Mutation.__init__(self, rate, final_rate, number, settings)

        self._fill_common_gene_list(settings)
        self.unique_gene_list = []
        self._fill_unique_list(settings)

    def reproduce(self, genome):
        """
        Mutates solution in the optimization.
        """
        child_genome = copy.deepcopy(genome)
        iteration = 0
        while child_genome == genome:
            for i in range(self.number):
                mutate = random.randint(0,len(child_genome)-1)
                old_gene = child_genome[mutate]
                new_gene = None
                if type(self.common_gene_list) == dict:
                    for common in self.common_gene_list:
                        lisst = self.common_gene_list[common]
                        if len(lisst) > 0:
                            if old_gene in lisst: 
                                temp = random.choice(lisst)
                                if temp in self.unique_gene_list:
                                    if temp in child_genome:
                                        pass
                                    else:
                                        new_gene = temp
                                else:
                                    new_gene = temp
                    if not new_gene:
                        new_gene = old_gene
                elif type(self.common_gene_list) == list:
                    if old_gene in self.common_gene_list:
                        new_gene = random.choice(self.common_gene_list)
                    else:
                        new_gene = copy.copy(old_gene)
                else:
                    raise ValueError("Unsupported data type for common gene list.")
                if self.genome_map[new_gene][mutate] == 1:
                    child_genome[mutate] = new_gene

            iteration += 1

        return child_genome

    def _fill_common_gene_list(self, settings):
        """
        Fills in the parameter common_gene_list. If the common gene list is provided it
        directly takes this value. Otherwise, it determines the common gene list from
        the provided gene maps of the chromosomes.

        Parameters:
            settings: Dictionary
                The input yaml file, as a dictionary.

        Written by Brian Andersen. 1/21/2020.
        """
        MUTE_SETTING = settings['optimization']['mutation']['common_chromosomes']
        if MUTE_SETTING == 'determine_from_maps':
            chromosomes = settings['genome']['chromosomes']
            examine_dictionary = {}
            current_map_count = 0
            for gene in chromosomes:
                if not examine_dictionary:
                    examine_dictionary[current_map_count] = {}
                    examine_dictionary[current_map_count]['map'] = chromosomes[gene]['map']
                    examine_dictionary[current_map_count]['gene_list'] = [gene]
                    current_map_count += 1
                else:
                    no_match_found = True
                    for i in examine_dictionary:
                        if examine_dictionary[i]['map'] == chromosomes[gene]['map']:
                            examine_dictionary[i]['gene_list'].append(gene)
                            no_match_found = False
                            break
                    if no_match_found:
                        examine_dictionary[current_map_count] = {}
                        examine_dictionary[current_map_count]['map'] = chromosomes[gene]['map']
                        examine_dictionary[current_map_count]['gene_list'] = [gene]
                        current_map_count += 1

            common_list = []
            for i in examine_dictionary:
                if len(examine_dictionary[i]['gene_list']) > 1:
                    common_list[i] = examine_dictionary[i]['gene_list']
            self.common_gene_list = common_list
        elif MUTE_SETTING == 'determine_from_fixity':
            chromosomes = settings['genome']['chromosomes']
            self.common_gene_list = {}
            for gene in chromosomes:
                group = chromosomes['gene']['gene_group']
                if group in self.common_gene_list:
                    self.common_gene_list[group].append(gene)
                else:
                    self.common_gene_list[group] = [gene]
        else:
            self.common_gene_list = copy.deepcopy(MUTE_SETTING)

    def _fill_unique_list(self,settings):
        """
        This is an idea for genes of a unique type, meaning that only one is allowed in the optimization solution.
        It is thought of as a way to include burnt fuel assemblies. Right now, each burnt fuel assembly is its own 
        group. But I was thinking what if I put them into a group. Then, in case you include a bad assembly, such as low 
        fuel enrichment or whatever, it could be easily selected out.
        """
        GENES = settings['genome']['chromosomes']
        for key in GENES:
            if 'unique' in GENES[key]:
                if GENES[key]['unique']:
                    self.unique_gene_list.append(key)


class Mutate_By_Unique(Mutation):
    """
    Performs mutation preserving genes uniqueness. 
    Inherits from class Mutation.
    """
    def __init__(self, rate, final_rate, number, settings=None):
        Mutation.__init__(self, rate, final_rate, number, settings)

        self.unique_gene_list = []
        self.nonunique_gene_list = []
        self._fill_unique_list(settings)
        self.gene_list =  self.unique_gene_list + self.nonunique_gene_list

    def reproduce(self, genome):
        """
        Mutates solution in the optimization.
        """
        child_genome = copy.deepcopy(genome)
        iteration = 0
        SWAP = False
        while child_genome == genome:
            for i in range(self.number):
                mutate = random.randint(0,len(child_genome)-1)
                old_gene = child_genome[mutate]
                new_gene = None
                temp = random.choice(self.gene_list)
                if new_gene in self.nonunique_gene_list:
                    new_gene = temp
                else:
                    if temp not in child_genome:
                        new_gene = temp
                    else:
                        SWAP = True
                        mutate1 = mutate
                        gene1 = copy.deepcopy(old_gene)
                        mutate2 = np.where(np.array(child_genome)==temp)[0][0]
                        new_gene = temp
                        gene2 = copy.deepcopy(temp)  
                if self.genome_map[new_gene][mutate] == 1:
                    if SWAP:
                        child_genome[mutate1] = gene2
                        child_genome[mutate2] = gene1
                    else:
                        child_genome[mutate] = new_gene
            iteration += 1

        return child_genome
    
    def _fill_unique_list(self,settings):
        """
        This is an idea for genes of a unique type, meaning that only one is allowed in the optimization solution.
        It is thought of as a way to include burnt fuel assemblies. Right now, each burnt fuel assembly is its own 
        group. But I was thinking what if I put them into a group. Then, in case you inclue a bad assembly, such as low 
        fuel enrichment or whatever, it could be easily selected out.
        """
        GENES = settings['genome']['chromosomes']
        for key in GENES:
            if 'unique' in GENES[key]:
                if GENES[key]['unique']:
                    self.unique_gene_list.append(key)
                else:
                    self.nonunique_gene_list.append(key)
            else:
                self.nonunique_gene_list.append(key)


class Mutate_By_Genome(Mutation):
    """
    Mutates Genomes based upon the solution space of the genomes.
    Inherits from class Mutation.
    """
    def __init__(self, rate, final_rate, number, settings):
        Mutation.__init__(self, rate, final_rate, number, settings)

    def reproduce(self, genome):
        """
        Generates new solution through mutation in the optimization.
        """

        child_genome = copy.deepcopy(genome)

        while child_genome == genome:
            for i in range(self.number):
                mutate = random.randint(0, len(child_genome)-1)
                old_gene = child_genome[mutate]

                key_list = list(self.genome_map.keys())
                new_gene = random.choice(key_list)
                if new_gene == old_gene:
                    pass
                else:
                    if self.genome_map[new_gene][mutate] == 1:
                        child_genome[mutate] = new_gene

        return child_genome


class Mutate_By_Common(Mutation):
    """
    Genes will only mutate into other genes already present in the solution to the
    optimization problem.
    Inherits from class Mutation.
    """
    def __init__(self, rate, final_rate, number, settings):
        Mutation.__init__(self, rate, final_rate, number, settings)

    def reproduce(self, genome):
        """
        Generates new solution through mutation in the optimization.
        """
        child_genome = copy.deepcopy(genome)

        gene_list = []
        for gene in child_genome:
            if gene in gene_list:
                pass
            else:
                gene_list.append(gene)

        while child_genome == genome:
            for i in range(self.number):
                mutate = random.randint(0, len(child_genome)-1)
                old_gene = child_genome[mutate]

                new_gene = random.choice(gene_list)

                if new_gene == old_gene:
                    pass
                else:
                    if self.genome_map[new_gene][mutate] == 1:
                        child_genome[mutate] = new_gene

        return child_genome


class Submutation_Mutator(Mutation):
    """
    Mutation type used when genes have subgenes, i.e. how used fuel assemblies
    are implemented.

    All genomes must be defined through a dictionary in this case.
    """
    def __init__(self, rate, final_rate, number, sub_dict,settings):
        Mutation.__init__(self, rate, final_rate, number,settings)

        self.mutation_key = {}
        for key in sub_dict:
            if 'mutate_by_type' == sub_dict[key]:
                self.mutation_key[key] = Mutate_By_Type(rate, final_rate, number,
                                                        sub_dict['common_gene_list'])
            elif 'mutate_by_genome' ==  sub_dict[key]:
                self.mutation_key[key] = Mutate_By_Genome(rate, final_rate, number,
                                                          sub_dict['map_genome'])
            elif 'mutate_by_common' == sub_dict[key]:
                self.mutation_key[key] = Mutate_By_Common(rate, final_rate, number,
                                                          sub_dict['map_genome'])
            elif 'mutate_fixed' == sub_dict[key]:
                self.mutation_key[key] = Fixed_Genome_Mutator(rate, final_rate, number,
                                                              sub_dict['map_genome'])

    def reproduce(self, solution):

        key_list = list(solution.genome.keys())
        child_genome = copy.deepcopy(solution.genome)
        key = random.choice(key_list)
        encode = child_genome[key]
        encode = self.mutation_key[key].reproduce(encode)
        child_genome[key] = encode
        return child_genome


class Fixed_Genome_Mutator(Mutation):
    def __init__(self, rate, final_rate, number, settings):
        Mutation.__init__(self, rate, final_rate, number,settings)

    def reproduce(self, genome):
        """
        Rearranges the genes in the genome.
        """
        non_symmetries = []
        for i in range(len(genome)):
            if i in self.symmetry_list:
                pass
            else:
                non_symmetries.append(i)
        child_genome = copy.deepcopy(genome)
        while child_genome == genome:
            for i in range(self.number):
                mutate1 = random.randint(0, len(child_genome)-1)
                if mutate1 in self.symmetry_list: 
                    mutate2 = random.choice(self.symmetry_list)
                    gene_one = child_genome[mutate1]
                    gene_two = child_genome[mutate2]
                    if self.genome_map[gene_one][mutate2] == 1:
                        if self.genome_map[gene_two][mutate1] == 1:
                            child_genome[mutate2] = gene_one
                            child_genome[mutate1] = gene_two
                else:
                    mutate2 = random.choice(non_symmetries)
                    gene_one = child_genome[mutate1]
                    gene_two = child_genome[mutate2]
                    if self.genome_map[gene_one][mutate2] == 1:
                        if self.genome_map[gene_two][mutate1] == 1:
                            child_genome[mutate2] = gene_one
                            child_genome[mutate1] = gene_two
        return child_genome


class Double_Mutator(Mutation):
    """
    Two different types of mutation are possible using this method. 
    """
    def __init__(self, rate, final_rate, number, first_method, second_method, settings=None):
        Mutation.__init__(self, rate, final_rate, number, settings)
        self.first_method = first_method(rate, final_rate, number, settings=settings)
        self.first_method.genome_map = self.genome_map
        self.second_method = second_method(rate, final_rate, number, settings=settings)
        self.second_method.genome_map = self.genome_map

    def reproduce(self, genome):
        """
        Originally the specific method that performs mutation. Has now been
        altered to choose which mutation method is performed. For now, written so
        that there is a 50/50 chance between the two methods.

        Parameters:
            genome: list
                The genome that is to be mutated.

        Written by Brian Andersen. 1/21/2020
        """
        chance = random.choice([True, False])
        if chance:
            child_genome = self.first_method.reproduce(genome)
        else:
            child_genome = self.second_method.reproduce(genome)

        return child_genome