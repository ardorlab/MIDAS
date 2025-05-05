## Import Block ##
import logging
from copy import deepcopy
import random
import numpy as np
from midas.utils import optimizer_tools as optools
from midas.algorithms.simulated_annealing import Cooling_Schedule as SA_Cooling_Schedules
from midas.algorithms.simulated_annealing import SA_reproduction 
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

    Written by Jake Mikouchi. 04/25/2025
    """

    def __init__(self, input):
        self.input = input   
        self.global_temperature = input.initial_temperature
        self.local_temperatures = [input.initial_temperature for i in range(self.input.num_procs)]
        self.buffer = []
        self.active_solutions = [[] for i in range(self.input.num_procs)]
        self.selected_solutions = [0 for i in range(self.input.num_procs)]
        self.best_in_gen = [0 for i in range(self.input.num_procs)]
        self.total_moves = 0

    def reproduction(self, population_step, proc):
        """
        Generates a new individual by perturbing the origional solution. The perturbation 
        operates in a similar manner to mutations in genetic algorithm. 
        Specific for PSA
        
        Updated by Jake Mikouchi. 04/26/2025
        """
        if population_step == 0:
            # takes solution from buffer if begining of serial SA
            primary_individual = [PSA_reproduction.select_from_buffer(self, proc).chromosome]
            self.total_moves = 0

        else: 
            # if not beginning of algorithm then use stored active solution
            primary_individual = [PSA_reproduction.selection(self, proc).chromosome]

        individual_pairs = deepcopy(primary_individual)
        ## preserve core parameters 
        LWR_core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]
        ## Perform perturbation
        if self.input.perturbation_type == "perturb_by_gene":
            individual_pairs.append(SA_reproduction.perturb_by_gene(self.input, primary_individual[0]))
        else:
            raise ValueError("Requested perturbation type not recognized.")

        self.active_solutions[proc] = individual_pairs

        self.local_temperatures[proc] = self.Temperature_update_methods(self.local_temperatures[proc], self.input.secondary_cooling_schedule, population_step, Global=False)

        return individual_pairs

    def update_active(self, full_active_solutions):
        """
        updates active_solutions with the final values

        Written by Jake Mikouchi. 04/29/2025
        """      
        self.active_solutions = full_active_solutions

    def update_buffer(self, pop_list, current_generation):
        """
        Updates solution buffer by filling it with solutions or replacing old solutions with new ones

        Written by Jake Mikouchi. 04/25/2025
        """

        if not self.buffer:
            for i in range(self.input.buffer_size):
                self.buffer.append(pop_list[i])

            self.global_temperature = PSA_Cooling_Schedules.lam_set_intial(self)
            logger = logging.getLogger("MIDAS_logger")
            logger.info(f"Initial Temperature: {self.global_temperature}")

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

            # calculate temperature using LAM cooling schedule
            logger = logging.getLogger("MIDAS_logger")
            self.global_temperature = self.Temperature_update_methods(self.global_temperature, self.input.cooling_schedule, current_generation, Global=True)
            logger.info(f"Updated Temperature: {self.global_temperature}")

        for i in range(len(self.local_temperatures)):
            self.local_temperatures[i] = self.global_temperature

    def Temperature_update_methods(self, temperature, cooling_schedule, current_step, Global):
        """
        Method for distributing to the requested SA cooling schedule

        NOTE: If a future developer adds a new cooling schedule in SA that can be utilized in
        PSA, it has to be explicitly added to this function in order for PSA to utilize it. 

        updated by Jake Mikouchi. 04/25/2024
        """

        if Global:
            # updates the global temperatures
            if cooling_schedule == 'exponential_decrease':
                temperature = SA_Cooling_Schedules.exponential_decrease(temperature)
            if cooling_schedule == 'linear_update':
                temperature = SA_Cooling_Schedules.linear_update(self.input.initial_temperature, current_step, self.input.num_generations)
            if cooling_schedule == 'log_update':
                temperature = SA_Cooling_Schedules.logarithmic_update(self.input.initial_temperature, current_step)
            if cooling_schedule == 'lam':
                temperature = PSA_Cooling_Schedules.lam(self, current_step)  
        else: 
            # updates the local temperatures
            if cooling_schedule == 'exponential_decrease':
                temperature = SA_Cooling_Schedules.exponential_decrease(temperature)
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
        """

        scalingfactor = 100
        # finds the sum of the probabilities of each individual being selected in the Buffer
        sumprob = 0
        for i in range(self.input.buffer_size):
            sumprob += np.exp((-1 * self.buffer[i].fitness_value) / scalingfactor)

        # finds the unique probability of each individual being selected in the Buffer
        positional_prob = []
        for i in range(self.input.buffer_size):
            currprob = np.exp((-1 * self.buffer[i].fitness_value) / scalingfactor)
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

    def selection(self, proc):
        """
        Selects the current indivdiual in the SA optimization.
        This is a little weird and may need to be addressed. MIDAS is constructed to always maximize the fitness
        but SA works by minimizing the cost. For now maintianing GA likeness is the priority so cost/fitness is maximized.
        
        Created by Jake Mikouchi. 04/22/2025
        """

        with open("info.txt","a") as f:
            for soln in self.selected_solutions:
                f.write(str(soln.fitness_value)+"  ")
            f.write("\n")

        try:
            if self.selected_solutions[proc].chromosome == self.active_solutions[proc][0].chromosome:
                primary = self.active_solutions[proc][0]
                challenger = self.active_solutions[proc][1]
            if self.selected_solutions[proc].chromosome == self.active_solutions[proc][1].chromosome:
                primary = self.active_solutions[proc][1]
                challenger = self.active_solutions[proc][0]
        except: 
            primary = self.active_solutions[proc][0]
            challenger = self.active_solutions[proc][0]
            
        # challenger = self.active_solutions[proc][0]
        # if len(self.active_solutions[proc]) < 2:
        #     primary = self.active_solutions[proc][0]
        # else:
        #     primary = self.active_solutions[proc][1]

        selected = self.active_solutions[proc][0]

        if challenger.fitness_value >= primary.fitness_value:
            selected = challenger
        else: 
            acceptance_prob = np.exp(-1 * (primary.fitness_value - challenger.fitness_value) / self.local_temperatures[proc]) 
            chance = random.random()
            if chance < acceptance_prob:
                selected = challenger
                self.total_moves += 1
            else: 
                selected = primary

        if challenger.fitness_value >= self.best_in_gen[proc].fitness_value:
            self.best_in_gen[proc] = challenger

        self.selected_solutions[proc] = selected

        return selected


class PSA_Cooling_Schedules(object):
    """
    Class for Parallel Simulated Annealing cooling schedules.

    The cooling schedule sets the tolerance for accepting new solutions.
    Thisdictates the "randomness" of the optimization and balances the exploration vs exploitation 
    of the optimization. Generally, it is best for the cooling schedule to
    start at a high temperature and gradually decrease throughout the optimization.

    The PSA algorithm can access all cooling schedules available in the SA script. 
    This class is for PSA specific schedules which SA cannot utilize.

    If a future developer wants to add cooling schedules that can be universally applied to both 
    SA and PSA please add it to the simulated_annealing.py file.

    Updated by Jake Mikouchi 04/23/2025
    """

    def __init__(self, generation):
        self.generation = generation

    def lam(self, current_generation):
        """
        lam cooling schedule. adaptively updates the temperature acroding to statistics gathered
        during the optimization. This is unique to PSA and is not available for SA.
        Look to https://nstopenresearch.org/articles/2-5 for more information on how it works

        updated by Jake Mikouchi. 4/30/2025
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
        updated by Jake Mikouchi. 4/30/2025
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

  