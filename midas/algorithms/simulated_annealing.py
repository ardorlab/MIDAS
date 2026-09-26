## Import Block ##
import logging
from copy import deepcopy
import random
import numpy as np
from midas.utils import optimizer_tools as optools
from midas.utils.mutator import Mutator

class Simulated_Annealing():
    """
    Class for performing optimization using the simulated annealing.
    Simulated annealing works by iterating over a single solution and continuously perturbing it.
    Note that within MIDAS each iterations is refered to as a generation. This is not necessarily the correct nomenclature 
    for SA but it is named in this way for consistencey.

    Written by Brian Andersen. 1/9/2020
    Updated by Jake Mikouchi. 04/21/2025
    """

    def __init__(self, input):
        self.input = input   
        self.temperature = input.initial_temperature
        self.generation = 0 
        self.selected_solution = 0  


    def reproduction(self, pop_list, current_generation):
        """
        Generates a new individual by perturbing the origional solution. The perturbation 
        operates in a similar manner to mutations in genetic algorithm. 
        
        Updated by Jake Mikouchi. 04/21/2025
        """
        logger = logging.getLogger("MIDAS_logger")
        
    ## Container for holding new list of child chromosomes
        primary_individual = [SA_reproduction.selection(self, self.temperature, pop_list).chromosome]
        individual_pairs = deepcopy(primary_individual)
    ## preserve core parameters 
        core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]
    ## Perform perturbation
        individual_pairs.append(Mutator.mutator_methods(self.input, primary_individual[0]))
        self.temperature = SA_reproduction.Temperature_update_methods(self, self.temperature, self.input.cooling_schedule)
        logger.info(f"Updated Temperature: {self.temperature}")
        self.generation += 1

        return individual_pairs


class SA_reproduction():
    """
    Functions for performing reproduction of chromosomes using SA methodologies. 
    The name is kept as "reproduction" for consistency across algorithms
     
    Written by Jake Mikouchi. 04/21/25
    Updated by Jake Mikouchi. 09/26/26
    """

    def selection(self, temperature, pop_list):
        """
        Selects the current indivdiual in the SA optimization.
        
        Created by Jake Mikouchi. 04/22/2025
        """

        # optimizer.py does some weird shifting due to inactive solutions
        # so challenger is index 0 while primary is index 1

        try:
            if self.selected_solution.chromosome == pop_list[0].chromosome:
                primary = pop_list[0]
                challenger = pop_list[1]
            if self.selected_solution.chromosome == pop_list[1].chromosome:
                primary = pop_list[1]
                challenger = pop_list[0]

        except: 
            primary = pop_list[0]
            challenger = pop_list[0]

        selected = pop_list[0]

        if challenger.fitness_value >= primary.fitness_value:
            selected = challenger
        else: 
            acceptance_prob = np.exp(-1 * (primary.fitness_value - challenger.fitness_value) / temperature)
            chance = random.random()
            if chance < acceptance_prob:
                selected = challenger
            else: 
                selected = primary  

        self.selected_solution = selected
        return selected

    def Temperature_update_methods(self, temperature, cooling_schedule):
        """
        Method for distributing to the requested cooling schedule method.
        
        updated by Jake Mikouchi. ~spring 2025
        """
        if cooling_schedule == 'exponential_decrease':
            temperature = Cooling_Schedule.exponential_decrease(self.input.update_factor, temperature)
        if cooling_schedule == 'linear_update':
            temperature = Cooling_Schedule.linear_update(self.input.initial_temperature, self.generation, self.input.num_generations)
        if cooling_schedule == 'log_update':
            temperature = Cooling_Schedule.logarithmic_update(self.input.initial_temperature, self.generation)

        return temperature 

class Cooling_Schedule(object):
    """
    Class for Simulated Annealing cooling schedules.
    The cooling schedule sets the tolerance for accepting new solutions.
    The cooling schedule dictates the "randomness" of the optimization and balances the 
    exploration vs exploitation of the optimization. Generally, it is best for the cooling schedule to
    start at a high temperature and gradually decrease throughout the optimization.
    All cooling schedules shown here can be accessed by both SA and PSA.
    Updated by Jake Mikouchi 04/23/2025
    """

    def __init__(self, generation):
        self.generation = generation

    def exponential_decrease(update_factor, temperature):
        """
        T = T0*alpha
        Where 0.9 < alpha < 1.0 
        
        Updated by Jake Mikouchi 04/23/2025
        Updated by Jake Mikouchi 09/13/2025
        """
        if temperature <= 0.0001:
            temperature = 0.0001
        else:
            temperature = temperature * update_factor
        return temperature

    def linear_update( initial_temperature, current_generation, total_generations):
        """
        linearly updates the temperature
        
        created by Jake Mikouchi 04/23/2025
        """
        temperature = initial_temperature + ((0 - initial_temperature) / total_generations) * (current_generation + 1)

        return temperature

    def logarithmic_update(initial_temperature, current_generation):
        """
        Logarithmically updates the temperature
        Note that the user defined inital temperature is used as a contant rather than the actual starting point.
        
        created by Jake Mikouchi 04/23/2025
        """
        temperature = initial_temperature / np.log10(2 + current_generation)

        return temperature