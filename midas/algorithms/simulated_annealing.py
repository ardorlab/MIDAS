## Import Block ##
import logging
from copy import deepcopy
import random
import numpy as np
from midas.utils import optimizer_tools as optools

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
    ## Container for holding new list of child chromosomes
        primary_individual = [SA_reproduction.selection(self, self.temperature, pop_list).chromosome]
        individual_pairs = deepcopy(primary_individual)
    ## preserve core parameters 
        LWR_core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]
    ## Perform perturbation
        if self.input.perturbation_type == "perturb_by_gene":
            individual_pairs.append(SA_reproduction.perturb_by_gene(self.input, primary_individual[0]))
        else:
            raise ValueError("Requested perturbation type not recognized.")

        self.temperature = SA_reproduction.Temperature_update_methods(self, self.temperature, self.input.cooling_schedule)
        self.generation += 1

        return individual_pairs


class SA_reproduction():
    """
    Functions for performing reproduction of chromosomes using SA
    methodologies, currently only contains one perturbation method.
    The perturbation method works in the same way as mutations in GA.
    The name is kept as "reproduction" for consistency across algorithms
     
    Written by Jake Mikouchi. 04/21/25
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
                primary = pop_lis[0]
                challenger = pop_list[0][1]
            if self.selected_solution.chromosome == pop_list[1].chromosome:
                primary = pop_list[1]
                challenger = pop_list[0]

        except: 
            primary = pop_list[0]
            challenger = pop_list[0]

        # challenger = pop_list[0]
        # if len(pop_list) < 2:
        #     primary = pop_list[0]
        # else:
        #     primary = pop_list[1]

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
        logger = logging.getLogger("MIDAS_logger")

        if cooling_schedule == 'exponential_decrease':
            temperature = Cooling_Schedule.exponential_decrease(temperature)
        if cooling_schedule == 'linear_update':
            temperature = Cooling_Schedule.linear_update(self.input.initial_temperature, self.generation, self.input.num_generations)
        if cooling_schedule == 'log_update':
            temperature = Cooling_Schedule.logarithmic_update(self.input.initial_temperature, self.generation)
        logger.info(f"Updated Temperature: {temperature}")

        return temperature 

## Mutation types ##
    def perturb_by_gene(input_obj, chromosome):
        """
        Generates a new solution by randomly mutating a single gene.
        
        Created by Jake Mikouchi. 04/21/2025
        """
        ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")

        core_parameters = [input_obj.nrow, input_obj.ncol, input_obj.num_assemblies, input_obj.symmetry, input_obj.calculation_type]

        if input_obj.calculation_type in ["eq_cycle"]:
            zone_chromosome = [loc[0] for loc in chromosome]
            child_zone_chromosome = deepcopy(zone_chromosome)
            old_soln = zone_chromosome
            new_soln = child_zone_chromosome
            all_gene_options = input_obj.batches
            all_genes_list = list(input_obj.batches.keys())
        else:
            child_chromosome = deepcopy(chromosome)
            old_soln = chromosome
            new_soln = child_chromosome
            all_gene_options = input_obj.genome
            all_genes_list = list(input_obj.genome.keys())

        num_mutations = 1 #!TODO: this was hardcoded to 1 in old MIDAS. Should probably be parameterized.
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            new_soln = deepcopy(old_soln) #in the case of abortion, start from scratch.
            while new_soln == old_soln:
                for i in range(num_mutations):
                    loc_to_mutate = random.randint(0, len(new_soln)-1) #choose a random gene
                    old_gene = new_soln[loc_to_mutate]
                    gene_options = optools.Gene_Validity_check.contraceptive_check(input_obj, all_genes_list, all_gene_options,
                                                                                    core_parameters, old_soln, [], loc_to_mutate)
                    new_gene = random.choice(gene_options)
                    if new_gene != old_gene:
                        if all_gene_options[new_gene]['map'][loc_to_mutate] == 1:
                            new_soln[loc_to_mutate] = new_gene
            chromosome_is_valid = optools.Gene_Validity_check.abortive_check(input_obj,all_genes_list,all_gene_options,\
                                                                            core_parameters,new_soln)
            if not chromosome_is_valid:
                attempts += 1
                if attempts > 100000:
                    logger.error("Mutate-by-Gene has failed after 100,000 attempts; the Individual will be restored. Consider relaxing the constraints on the input space.")
                    return chromosome

        if input_obj.calculation_type in ["eq_cycle"]:
            #recreate child_chromosome
            child_chromosome = []
            for i in range(len(new_soln)):
                if new_soln[i] == chromosome[i][0]:
                    child_chromosome.append(chromosome[i])
                else:
                    child_chromosome.append((new_soln[i],None))
            child_chromosome = optools.Solution.EQ_reload_fuel(input_obj.genome,core_parameters,child_chromosome)

        else: 
            child_chromosome = new_soln

        return child_chromosome


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

    def exponential_decrease(temperature):
        """
        T = T0*alpha
        Where 0.9 < alpha < 1.0 
        
        Updated by Jake Mikouchi 04/23/2025
        """
        alpha = 0.95 #TODO make this avialble to edit in input file
        if temperature <= 0.0001:
            temperature = 0.0001
        else:
            temperature = temperature * alpha
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