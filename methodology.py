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
import fitness
from multiprocessing import Pool
from solution_types import evaluate_function,Unique_Solution_Analyzer,test_evaluate_function
from metrics import Optimization_Metric_Toolbox

class Solution_Generator(object):
    '''
    Class for generating solutions to the optimization problem.
    Note, only the solution portion of the optimization is generated. A class
    instance is not created.
    E.g., the genome for initial_parent_0 is created, but the solution initial_parent_0 
    is not.

    Written by Brian Andersen 4/8/2021
    '''
    def __init__(self):
        self.solution = None

    def return_solution_length(self,decision_variables):
        """
        Checks that all decision variable maps are of the same size. If they are, returns the 
        decision variable list and length.

        Written by Brian Andersen. 4/8/2021
        """
        solution_length = None
        self.variable_list = list(decision_variables.keys())
        if 'symmetry_list' in self.variable_list:
            self.variable_list.remove('symmetry_list')
        
        for variable in self.variable_list:
            if not solution_length:
                self.solution_length = len(decision_variables[variable]['map'])
            elif len(decision_variables[variable]['map']) == self.solution_length:
                pass
            else:
                error_message = 'Number of decisions for the Decision'
                error_message += f' Variable Map for {variable} is not '
                error_message += f'consistent with previous maps.'
                raise ValueError(error_message)

    def variables_in_group(self,decision_variables,group_name):
        """
        Returns a list of the genes in the chosen group
        """
        variable_list = []
        for variable in decision_variables:
            if variable == 'symmetry_list':
                pass
            else:
                if group_name == decision_variables[variable]['gene_group']:
                    variable_list.append(variable)

        return variable_list

    def is_variable_ok(self,decision_variables,variable,space):
        """
        Checks if the variable is allowed in the desired location
        """
        variable_is_ok = True
        if not decision_variables[variable]['map'][space]:
            variable_is_ok = False
        if space in decision_variables['symmetry_list']:
            if self.solution_groups[decision_variables[variable]['gene_group']] <= 1:
                variable_is_ok = False
        else:
            if not self.solution_groups[decision_variables[variable]['gene_group']]:
                variable_is_ok = False
        if 'unique' in decision_variables[variable]:
            if variable in self.solution:
                variable_is_ok = False

        return variable_is_ok

    def generate_non_fixed_solution(self,decision_variables):
        """
        Generates a random solution to the optimization problem without any
        restrictions on decision variable placement in the optimizaiton problem.
        
        Written by Brian Andersen. 4/8/2021.
        """
        self.return_solution_length(decision_variables)

        self.solution = []
        for i in range(self.solution_length):
            no_variable_found = True
            while no_variable_found:
                variable = random.choice(self.variable_list)
                if decision_variables[variable]['map']:
                    self.solution.append(variable)
                    no_variable_found = False

    def generate_fixed_solution(self,decision_variables,variable_groups):
        """
        Generates a random solution to the optimizaiton problem with restrictions on decision variable  
        placement in the optimization problem.

        Written by Brian Andersen. 4/8/2021
        """
        self.return_solution_length(decision_variables)

        no_solution_found = True
        while no_solution_found:
            attempts = 0
            self.solution_groups = copy.deepcopy(variable_groups)
            self.solution = [None]*self.solution_length
            unfilled_spaces = list(range(self.solution_length))
            while unfilled_spaces:  
                space_number = random.randint(0,len(unfilled_spaces)-1)
                group_name = None
                while not group_name:
                    random_group = random.choice(list(self.solution_groups.keys()))
                    if self.solution_groups[random_group] > 0:
                        group_name = random_group
                available_variable_list = self.variables_in_group(decision_variables,group_name)
                space = unfilled_spaces[space_number]
                variable = random.choice(available_variable_list)
                variable_is_ok = self.is_variable_ok(decision_variables,variable,space)
                if variable_is_ok:
                    self.solution[space] = variable
                    unfilled_spaces.remove(space)
                    if space in decision_variables['symmetry_list']:
                        self.solution_groups[decision_variables[variable]['gene_group']] -= 2
                    else:
                        self.solution_groups[decision_variables[variable]['gene_group']] -= 1             
                else:
                    attempts += 1
                if attempts == 100:
                    break

            bad_variable_list = []
            for i,variable in enumerate(self.solution):
                if not variable:
                    bad_variable_list.append(i)

            if not bad_variable_list:
                no_solution_found = False  

class Optimization(object):
    """
    Generic class for optimizations for containing shared methods such as generating initial
    solutions and deleting old solutions.

    Written by Brian Andersen. 4/8/2021.
    """
    def __init__(self,solution,file_settings):
        self.solution = solution
        
        self.possible_decision_variables = file_settings['genome']['chromosomes']
        self.decision_variable_groups = file_settings['optimization']['fixed_groups']
        self.objectives = file_settings['optimization']['objectives']
        self.file_settings = file_settings
        if 'fixed_geometry' in file_settings['optimization']:
            self.fixed_geometry = file_settings['optimization']['fixed_geometry']
        else:
            self.fixed_geometry = False
        
        self.delete_output = False
        self.delete_cax = False
        self.delete_log = False
        self.delete_cross_sections = False
        self.cs_library_name = None
        if 'cleanup' in file_settings['optimization']:
            if 'output' in file_settings['optimization']['cleanup']:
                self.delete_output = file_settings['optimization']['cleanup']['output']
            if 'cax_file' in file_settings['optimization']['cleanup']:
                self.delete_cax = file_settings['optimization']['cleanup']['cax_file']
            if 'log_file' in file_settings['optimization']['cleanup']:
                self.delete_log = file_settings['optimization']['cleanup']['log_file']
            if 'cross_sections' in file_settings['optimization']['cleanup']:
                self.delete_cross_sections = file_settings['optimization']['cleanup']['cross_sections']
                if self.delete_cross_sections:
                    self.cs_library_name = file_settings['genome']['assembly_data']['cs_library']
    
    def generate_initial_solutions(self,name):
        """
        Function for creating an initial, random solution to the optimization
        problem
        
        Written by Brian Andersen. 4/8/2021
        """
        generator = Solution_Generator()
        solution = self.solution()
        solution.name = name
        if self.fixed_geometry:
            generator.generate_fixed_solution(self.possible_decision_variables,
                                       self.decision_variable_groups)
        else:
            generator.generate_non_fixed_solution(self.possible_decision_variables)
        solution.genome = generator.solution
        solution.parameters = copy.deepcopy(self.objectives)
        solution.add_additional_information(self.file_settings)

        return solution

    def cleanup(self,consideration_list,do_not_delete_list):
        """
        Deletes solution results that are no longer relevant to the optimization.

        consideration_list: The list of all solutions that are to be examined for deletion.
        do_not_delete_list: The list of solutions that are still being used in the optimization
                            and should not be deleted. Used because I found it easier to not delete
                            the output files of my optimized solutions. I usually wanted to look at
                            them quickly, not rerun the evaluation code and then read through them.

        Written by Brian Andersen 4/8/2021.
        """
        for solution in consideration_list:
            if solution in do_not_delete_list:
                if self.delete_output:
                    if os.path.isfile(f"{solution}/{solution}.out"): #The way casmo output files 
                        os.remove(f"{solution}/{solution}.out")      #are named.
                    if os.path.isfile(f"{solution}/{solution}_sim.out"): #The way simulate output files 
                        os.remove(f"{solution}/{solution}_sim.out")      #are named.
                if self.delete_cax:
                    if os.path.isfile(f"{solution}/{solution}.cax"): #The way casmo output files 
                        os.remove(f"{solution}/{solution}.cax")      #are named.
                if self.delete_log:
                    if os.path.isfile(f"{solution}/{solution}.log"): #The way casmo output files 
                        os.remove(f"{solution}/{solution}.log")      #are named.
                if self.delete_cross_sections:
                    if os.path.isfile(f"{solution}/{self.cs_library_name}"):
                        os.remove(f"{solution}/{self.cs_library_name}")

    def main(self,num_procs,restart,test):
        """
        num_procs is the number of processors for the optimization to utilize.

        restart determines if the optimization is restarting from a previous run that failed before 
        completion. 
        
        The test command is used to test new functionality. It will not perform solution evaluation, 
        but will return random numbers for each parameter. Used for simple debugging tests for new 
        functionality.
        """
        raise NotImplementedError

class Genetic_Algorithm(Optimization):
    """
    Class for performing optimization through the genetic algorithm methodology.

    Written by Brian Andersen 4/12/2021.
    """
    def __init__(self,solution,file_settings,population,generation,reproduction,selection):
        Optimization.__init__(self,solution,file_settings)
        self.population = population
        self.generation = generation
        self.reproduction = reproduction
        self.selection = selection
    
    def main(self,num_procs,restart,test):
        """
        Main loop of the genetic algorithm.
        
        Originally a main in parallel and main in serial were used. These have been deprecated,
        because I don't think the benefits they bring are worth the pain of maintaining two 
        almost identical functions.

        Written by Brian Andersen. 4/12/2021.
        """
        opt = Optimization_Metric_Toolbox(restart,self.objectives) 

        if restart: #Optimization is restarting from a previous case
            with open('optimized_solutions.yaml') as current_solutions:
                solutions = yaml.safe_load(current_solutions)

                track_file = open(opt.human_track_file_name,'r')
                track_lines = track_file.readlines()
                track_file.close()

                for line in track_lines:
                    if "Current Optimization Genereration" in line:
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
        else: #Optimization is beginning from a fresh start
            for i in range(self.population.size):
                foo = self.generate_initial_solutions(f'initial_parent_{i}')
                self.population.parents.append(foo)

            for i in range(self.population.size):
                foo = self.generate_initial_solutions(f'initial_child_{i}')
                self.population.children.append(foo)

            pool = Pool(processes=num_procs)
            if test:
                self.population.parents = pool.map(test_evaluate_function, self.population.parents)
                self.population.children = pool.map(test_evaluate_function, self.population.children)
            else:
                self.population.parents = pool.map(evaluate_function, self.population.parents)
                self.population.children = pool.map(evaluate_function, self.population.children)
            opt.update_all_value_tracker(self.population.parents)
            opt.update_all_value_tracker(self.population.children)

            combined_population_list = []
            combined_population_list.extend(self.population.parents)
            combined_population_list.extend(self.population.children)
            self.population = self.selection.perform(self.population)
            opt.update_generation_track_file(self.population.parents)
            opt.update_human_track_file(self.population.parents)
            opt.generation_counter += 1
            self.cleanup(combined_population_list,self.population.parents)
            current_gen = 0

        for self.generation.current in range(current_gen,self.generation.total):
            self.population.children = self.reproduction.reproduce(self.population.parents, 
                                                                     self.solution)
            for i,solution in enumerate(self.population.children):
                solution.name = "child_{}_{}".format(self.generation.current, i)
                solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                solution.add_additional_information(self.file_settings)

            if test:
                self.population.children = pool.map(test_evaluate_function, self.population.children)
            else:
                self.population.children = pool.map(evaluate_function, self.population.children)
            opt.update_all_value_tracker(self.population.children)
            
            combined_population_list = []
            combined_population_list.extend(self.population.parents)
            combined_population_list.extend(self.population.children)
            self.population = self.selection.perform(self.population)
            opt.update_generation_track_file(self.population.parents)
            opt.update_human_track_file(self.population.parents)
            opt.generation_counter += 1
            self.cleanup(combined_population_list,self.population.parents)
            
            opt.record_optimized_solutions(self.population.parents)

        opt.end_of_optimization()

class Serial_Simulated_Annealing(Optimization):
    """
    Class for performing optimization through the Simulated Annealing methodology.
    Serial because this method is designed to analyze solutions one at a time. I
    do not know what differences exist in parallel simulated annealing so that 
    can be included as a class later.

    Written by Brian Andersen 4/12/2021.
    """
    def __init__(self,solution,
                      population,
                      generation,
                      cooling_schedule,
                      fitness,
                      mutation,
                      file_settings
                ):
        Optimization.__init__(self,solution,file_settings)
        self.population = population
        self.generation = generation
        self.cooling_schedule = cooling_schedule
        self.fitness = fitness
        self.mutation = mutation
        self.no_improvement_counter = 0

    def evaluate_best_solution_over_time(self,active):
        """
        I've noticed that the simulated annealing algorithm performs quite well in 
        finding desirable solutions, but that incorrectly set cooling schedules 
        cause the algorithm to occasionally wander from these desirable solutions.
        In order to fix this I am adding the best solution over time. It will be the
        best objective fitness ever encountered. After X number of generations, if 
        the solution hasn't improved beyond the current best solution over time the
        active solution will be reset to this solution.

        Written by Brian Andersen 4/13/2021.
        """
        if active.fitness < self.best_solution.fitness:
            self.best_solution = active
            self.no_improvement_counter = 0
        else:
            self.no_improvement_counter +=1
            if self.no_improvement_counter > self.cooling_schedule.no_improvement_limit:
                active = copy.deepcopy(self.best_solution)
                self.no_improvement_counter = 0

        return active

    def main(self,num_procs,restart,test):
        """
        Main function for the Serial Simulated Annealing Algorithm.

        Written by Brian Andersen 4/12/2021.
        """
        opt = Optimization_Metric_Toolbox(restart,self.objectives)

        if restart:     
            track_file = open('optimization_track_file.txt','r')
            track_lines = track_file.readlines()
            track_file.close()

            for line in track_lines:
                if "Current generation of the optimization" in line:
                    elems = line.strip().split()
                    current_gen = int(elems[-1]) + 1
            track_lines = None

            current_gen = current_gen / self.population.size
            with open('optimized_solutions.yaml') as current_solutions:
                solutions = yaml.safe_load(current_solutions)

            for solut in solutions:
                active = self.solution()
                active.name = solut
                active.genome = solutions[solut]['genome']
                active.parameters = solutions[solut]['parameters']
                active.fitness = float(solutions[solut]['fitness'])
            for i in range(0,current_gen):
                self.cooling_schedule.update()
        else:
            # defining the active solution 
            current_gen = 0
            active = self.generate_initial_solutions("intial_solution")
            if test:
                active.test_evaluate()
            else:
                active.evaluate()
            # will store the fit for active solution
            one = []
            one.append(active)
            one = self.fitness.calculate(one)
            self.best_solution = copy.deepcopy(active)
            opt.update_all_value_tracker(one)
            opt.update_generation_track_file(one)
            opt.update_human_track_file(one)

        for self.generation.current in range(current_gen,self.generation.total):
            for number in range(self.population.size):
                challenge = self.solution()
                challenge.genome = self.mutation.reproduce(active.genome)
                challenge.name = "solution_{}_{}".format(self.generation.current,number)
                challenge.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                challenge.add_additional_information(self.file_settings)
                if test:
                    challenge.test_evaluate()
                else:
                    challenge.evaluate()

                test = [challenge]
                test = self.fitness.calculate(test)
                
                opt.update_all_value_tracker(test)
                opt.update_generation_track_file([active,challenge])
                opt.update_human_track_file([active,challenge])
                #determining which solution makes the next generation 
                
                acceptance = numpy.exp(-1*(challenge.fitness-active.fitness)/self.cooling_schedule.temperature)
                if challenge.fitness < active.fitness:
                    active = copy.deepcopy(challenge)
                elif random.uniform(0,1) < acceptance:
                    active = copy.deepcopy(challenge)
                self.cleanup([active,challenge],[active,self.best_solution])
                self.evaluate_best_solution_over_time(active)
            self.cooling_schedule.update()