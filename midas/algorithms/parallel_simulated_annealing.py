## Import Block ##
import logging
from copy import deepcopy
import random
import numpy as np
import os
from midas.utils import optimizer_tools as optools
from midas.algorithms.simulated_annealing import Cooling_Schedule as SA_Cooling_Schedules
from midas.algorithms.simulated_annealing import SA_reproduction 

from itertools import repeat
from multiprocessing import Pool
from midas.codes import parcs342, parcs343
from midas.codes import nuscale_lut
from midas.codes import trace50p5
from midas.codes import polaris624

import statistics

class Parallel_Simulated_Annealing():
    """
    Class for performing optimization using parallel simulated annealing.
    parallel simulated annealing works by running multiple SAs in parallel.
    It was decided to make a unique file for PSA to keep the code clean and organized.

    Within MIDAS each iterations is refered to as a generation. This is not necessarily the correct nomenclature 
    for SA but it is named in this way for consistencey.

    *********** NOTE ************
    PSA borrows a large amount of functions from SA.
    Currently these functions are:
    perturb_by_genes in SA_reproduction
    all cooling schedules in SA_Cooling_Schedules class

    If a future developer wishes to add a new perturbation type please do so in SA_reproduction in 
    simulatted_annealing.py and ensure that reproduction is updated accordingly.

    If a future devloper wishes to add a new cooling schedule that can be utlized by both PSA and SA
    please do so in SA_Cooling_Schedules in simulatted_annealing.py and update Temperature_update_methods accordingly.

    If a future devloper wishes to add a new cooling schedule that can only be utlized by PSA, please do so in
    PSA_Cooling_Schedules in parallel_simulated_annealiing.py. Currently the LAM cooling schedule is the only
    CS unique to PSA.    

    *********** END NOTE ************

    Created by Jake Mikouchi. 09/26/2025
    """

    def __init__(self, input, eval_func):
        self.input = input   
        self.global_temperature = input.initial_temperature
        self.local_temperatures = [input.initial_temperature for i in range(self.input.num_procs)]
        self.buffer = []
        self.active_solutions = [[] for i in range(self.input.num_procs)]
        self.selected_solutions = [0 for i in range(self.input.num_procs)]
        self.best_in_gen = [0 for i in range(self.input.num_procs)]
        self.total_moves = 0
        self.holder = optools.Population(2, None)
        self.eval_func = eval_func
        self.fitness    = optools.Fitness()

    def reproduction(self, pop_list, gen_obj): #, population_step, proc):
        """
        Execution pathway for PSA. In this function, each of the SA threads are initialized and the data for each SA thread is stored/processed.
        The self object holds information for each of the different SA threads. 
        This returns a fully evaluated individual rather than just the chromosome.
        
        Created by Jake Mikouchi. 09/26/2025
        """
        logger = logging.getLogger("MIDAS_logger")

        # update / populate buffer
        self.update_buffer(pop_list, gen_obj.current)

        # calculate initial temperature only if lam
        if gen_obj.current == 1:
            if self.input.cooling_schedule == 'lam':
                self.global_temperature = PSA_Cooling_Schedules.lam_set_intial(self)
            logger.info(f"Initial Temperature: {self.global_temperature}")

        else: 
            # update global temperatures
            self.global_temperature = self.Temperature_update_methods(self.global_temperature, self.input.cooling_schedule, gen_obj.current, Global=True)
            logger.info(f"Updated Temperature: {self.global_temperature}")

        # distribute global temperature to local temperatures
        for proc in range(len(self.local_temperatures)):
            self.local_temperatures[proc] = self.global_temperature

        # store statistics
        self.total_moves = 0

        # select starting point from buffer
        self.active_solutions =  [[] for i in range(self.input.num_procs)]
        for proc in range(self.input.num_procs): 
            self.active_solutions[proc].append(PSA_reproduction.select_from_buffer(self, proc))

        # run PSA threads
        with Pool(self.input.num_procs) as p:
            holder = p.starmap(self.optimizer_thread,  zip(range(self.input.num_procs), repeat(gen_obj)))

        # retrieve information
        for i, (best, moves) in enumerate(holder):
            self.best_in_gen[i] = best
            self.total_moves += moves

        return self.best_in_gen


    def optimizer_thread(self, proc, gen_obj):
        """
        This is largely a copy of main in the Optimizer class. PSA needs its own version
        for each of the threads. Each SA thread computes a 'mini' optimization which
        provides new solutions and statistics to the greater PSA algorithm.

        NOTE any updates to main in the Optimizer class need to be reflected here. 
        
        Created by Jake Mikouchi. 09/26/2025
        """

        # initalize
        from midas.optimizer import  Optimizer # this has to be here
        best_in_current_gen = 0
        proc_total_moves = 0
        os.system(f'rm -rf ./{self.input.results_dir_name}/tmp_Proc_{proc}_*')

        # iterate for population size
        for iter_pops in range(self.input.population_size):
            proc_total_moves += self.SA_threads(self.active_solutions[proc], proc)
            self.holder.current = []
            for i in range(len(self.active_solutions[proc])):
                self.holder.current.append(Optimizer.generate_solution(self, f'tmp_Proc_{proc}_Gen_{gen_obj.current}_Indv_{iter_pops}', self.active_solutions[proc][i]))
            ## Evaluate fitness
            ## If chromosome exists in previous generations, skip call to external code.
            inactive_solutions = []
            for soln in self.holder.current:
                try:
                    soln_index = self.holder.archive['solutions'].index(soln.chromosome)
                    soln.fitness_value = self.holder.archive['fitnesses'][soln_index]
                    soln.parameters = self.holder.archive['parameters'][soln_index]
                    inactive_solutions.append(soln)
                    self.holder.current.remove(soln)
                except ValueError:
                    continue #chromosome is unique, do nothing.
            
            ## Execute and parse objective/constraint values
            self.holder.current = [self.eval_func(soln, self.input) for soln in self.holder.current]
            if 'cost_fuelcycle' in self.input.objectives.keys():
                for soln in self.holder.current:
                    soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, self.input)
            if 'av_fuelenrichment' in self.input.objectives.keys():
                for soln in self.holder.current:
                    soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, self.input)
            
            ## Calculate fitness from objective/constriant values
            for soln in self.holder.current:
                soln.fitness_value = self.fitness.calculate(soln.parameters)
            
            ## Recombine active and inactive solutions.
            for soln in inactive_solutions:
                self.holder.current.append(soln)

                    ## Archive results
            for soln in self.holder.current:
                self.holder.archive['solutions'].append(soln.chromosome)
                self.holder.archive['fitnesses'].append(soln.fitness_value)
                self.holder.archive['parameters'].append(soln.parameters)

            # save the best fitness from each thread
            self.active_solutions[proc] = self.holder.current
            for soln in self.active_solutions[proc]: 
                try:
                    if soln.fitness_value > best_in_current_gen.fitness_value:
                        best_in_current_gen =  soln
                except:
                    best_in_current_gen =  soln

        return best_in_current_gen, proc_total_moves

    def SA_threads(self, pop_list, proc):
        """
        Generates a new individual by perturbing the origional solution. The perturbation 
        operates in a similar manner to mutations in genetic algorithm. 

        This is largely a copy of reproduction in Simulated_Annealing class. PSA needs its own version
        because different arguments are provided and different outputs are given.
        
        Created by Jake Mikouchi. 09/26/2025
        """
        ## Container for holding new list of child chromosomes
        holder = PSA_reproduction.selection(self, proc, pop_list)
        primary_individual = [holder[0].chromosome]
        moved = holder[1]
        individual_pairs = deepcopy(primary_individual)
        ## preserve core parameters 
        core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]
         ## Perform perturbation
        if self.input.perturbation_type['method'] == "perturb_by_gene":
            individual_pairs.append(SA_reproduction.perturb_by_gene(self.input, primary_individual[0]))
        else:
            raise ValueError("Requested perturbation type not recognized.")
        # update temperature
        self.local_temperatures[proc] = SA_reproduction.Temperature_update_methods(self, self.local_temperatures[proc], self.input.secondary_cooling_schedule)

        self.active_solutions[proc] = individual_pairs
        
        # returns if the solution moved or not
        return moved


    def update_buffer(self, pop_list, current_generation):
        """
        Updates solution buffer by filling it with solutions or replacing old solutions with new ones

        Written by Jake Mikouchi. 04/25/2025
        """
        # populate if nothing is there
        if not self.buffer:
            for i in range(self.input.buffer_size):
                self.buffer.append(pop_list[i])

        else: 
            potential_candidates = []
            for candidate in self.best_in_gen:
                if candidate not in self.buffer:
                    potential_candidates.append(candidate)

            for candidate in potential_candidates:
                # find smallest fitness value in buffer
                smallest_indv = self.buffer[0]
                for buf_indv in self.buffer:
                    if buf_indv.fitness_value < smallest_indv.fitness_value:
                        smallest_indv = buf_indv

                # replace smallest fitness value if candidate has a higher one
                if candidate.fitness_value > smallest_indv.fitness_value:
                    self.buffer.pop(self.buffer.index(smallest_indv))
                    self.buffer.append(candidate) 

        for i in range(len(self.local_temperatures)):
            self.local_temperatures[i] = self.global_temperature

    def Temperature_update_methods(self, temperature, cooling_schedule, current_step, Global):
        """
        Method for distributing to the requested SA cooling schedule

        NOTE: If a future developer adds a new cooling schedule in SA that can be utilized in
        PSA, it has to be explicitly added to this function in order for PSA to utilize it. 

        Created by Jake Mikouchi. 09/26/2025
        """

        if Global:
            # updates the global temperatures
            if cooling_schedule == 'exponential_decrease':
                temperature = SA_Cooling_Schedules.exponential_decrease(self.input.update_factor, temperature)
            if cooling_schedule == 'linear_update':
                temperature = SA_Cooling_Schedules.linear_update(self.input.initial_temperature, current_step-2, self.input.num_generations-1)
            if cooling_schedule == 'log_update':
                temperature = SA_Cooling_Schedules.logarithmic_update(self.input.initial_temperature, current_step)
            if cooling_schedule == 'lam':
                temperature = PSA_Cooling_Schedules.lam(self, current_step)  
        else: 
            # updates the local temperatures
            if cooling_schedule == 'exponential_decrease':
                temperature = SA_Cooling_Schedules.exponential_decrease(self.input.update_factor, temperature)
            if cooling_schedule == 'linear_update':
                temperature = SA_Cooling_Schedules.linear_update(self.global_temperature, current_step, self.input.population_size)
            if cooling_schedule == 'log_update':
                temperature = SA_Cooling_Schedules.logarithmic_update(self.global_temperature, current_step)
    
        return temperature 
    

class PSA_reproduction():
    """
    Functions for performing reproduction of chromosomes using PSA
    methodologies.
     
    Written by Jake Mikouchi. 04/25/25
    """

    def select_from_buffer(self, proc):
        """
        selects a solution from the buffer to act as a starting point for the serial SAs
        
        Created by Jake Mikouchi. 09/26/2025
        """

        scalingfactor = 1000
        # finds the sum of the probabilities of each individual being selected in the Buffer
        sumprob = 0
        for i in range(self.input.buffer_size):
            sumprob += np.exp((self.buffer[i].fitness_value) / scalingfactor)

        # finds the unique probability of each individual being selected in the Buffer
        positional_prob = []
        for i in range(self.input.buffer_size):
            currprob = np.exp((self.buffer[i].fitness_value) / scalingfactor)
            totalprob = (currprob / sumprob)
            positional_prob.append(totalprob)

       # uses the probabilities to select which indivdual to use as a solution
        sum_prob_list = []
        for i in range(self.input.buffer_size):
            if i >= 1:
                sum_prob_list.append(sum_prob_list[i - 1] + positional_prob[i])
            if i == 0:
                sum_prob_list.append(positional_prob[i])

        random_num = random.uniform(0, 1)

        for i in range(self.input.buffer_size):
            if i == 0:
                if random_num < sum_prob_list[i]:
                    selected_indv = self.buffer[i]
                    break
            if i >= 1:
                if random_num < sum_prob_list[i] and random_num >= sum_prob_list[i - 1]:
                    selected_indv  = self.buffer[i]
                    break

        self.best_in_gen[proc] = selected_indv

        self.selected_solutions[proc] = selected_indv


        return selected_indv

    def selection(self, proc, pop_list):
        """
        Selects the current indivdiual in the SA optimization.
        PSA must have a unique function because each thread has a unique set of data in the self object
        
        Created by Jake Mikouchi. 09/26/2025
        """
        moved = 0
        try:
            if self.selected_solutions[proc].chromosome == self.active_solutions[proc][0].chromosome:
                primary = pop_list[0]
                challenger = pop_list[1]
            if self.selected_solutions[proc].chromosome == self.active_solutions[proc][1].chromosome:
                primary = pop_list[1]
                challenger = pop_list[0]
        except: 
            primary = pop_list[0]
            challenger = pop_list[0]

        selected = pop_list[0]

        if challenger.fitness_value >= primary.fitness_value:
            selected = challenger
        else: 
            acceptance_prob = np.exp(-1 * (primary.fitness_value - challenger.fitness_value) / self.local_temperatures[proc]) 
            chance = random.random()
            if chance < acceptance_prob:
                selected = challenger
                moved = 1
            else: 
                selected = primary

        self.selected_solutions[proc] = selected

        return selected, moved


class PSA_Cooling_Schedules(object):
    """
    Class for Parallel Simulated Annealing cooling schedules.

    The cooling schedule sets the tolerance for accepting new solutions.
    This dictates the "randomness" of the optimization and balances the exploration vs exploitation 
    of the optimization. Generally, it is best for the cooling schedule to
    start at a high temperature and gradually decrease throughout the optimization.

    The PSA algorithm can access all cooling schedules available in the SA script. 
    This class is for PSA specific schedules which SA cannot utilize.

    If a future developer wants to add cooling schedules that can be universally applied to both 
    SA and PSA please add it to the simulated_annealing.py file.

    Updated by Jake Mikouchi 09/26/2025
    """

    def __init__(self, generation):
        self.generation = generation

    def lam(self, current_generation):
        """
        lam cooling schedule. adaptively updates the temperature acroding to statistics gathered
        during the optimization. This is unique to PSA and is not available for SA.
        Look to https://nstopenresearch.org/articles/2-5 for more information on how it works

        updated by Jake Mikouchi. 09/26/2025
        """

        Buffer_fitnesses = [soln.fitness_value for soln in self.buffer]
        try:
            deviation = statistics.stdev(Buffer_fitnesses)
        except:
            deviation = 1
        
        if deviation == 0: 
            deviation = 0.01

        # calculate move acceptanceratio
        p = self.total_moves / (self.input.num_procs * self.input.population_size)
        if p == 1:
            p = 0.9
        
        Gp = (4 * p * ((1 - p) ** 2)) / ((2 - p) ** 2)

        sk = (1 / self.global_temperature)
        sk1 = sk + (self.input.quality_factor * (1 / deviation) * (1 / ((sk ** 2) * (deviation ** 2)))) * Gp
        temperature = 1 / sk1

        return temperature 

    def lam_set_intial(self):
        """
        calculates the intial temperature using statistics gathered from the buffer
        updated by Jake Mikouchi. 09/26/2025
        """

        Buffer_fitnesses = [soln.fitness_value for soln in self.buffer]
        try:
            deviation = statistics.stdev(Buffer_fitnesses)
        except:
            deviation = 1
        
        if deviation == 0: 
            deviation = 0.01
        
        temperature = deviation * self.input.scaling_factor 

        return temperature

  