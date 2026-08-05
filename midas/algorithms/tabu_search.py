import logging
from copy import deepcopy
import random
import numpy as np
from midas.utils import optimizer_tools as optools


class Tabu_Search():
    """
    Class for performing tabu seach optimizations.
    Tabu search is a stochastic optimization algorithm exploring neighboring solutions and prohibiting certain moves
    using the tabu list.
    Note that within MIDAS each iterations is refered to as a generation. This is not necessarily the correct nomenclature 
    for TS but it is named in this way for consistencey.

    Written by Jake Mikouchi. 07/03/2026
    """

    def __init__(self, input):
        self.input = input   
        self.best_fitness = -np.inf
        self.prev_active_soln = []
        self.active_soln = []
        self.active_fitness = -np.inf
        self.tabu_list = []

    def reproduction(self, pop_list, current_generation):
        """
        """
        logger = logging.getLogger("MIDAS_logger")

        # select or update active solution
        if current_generation.current == 1: 
            self.active_soln = TS_reproduction.initial_selection(self, pop_list)   
        else:
            self.prev_active_soln = self.active_soln
            self.active_soln = TS_reproduction.update_active(self, pop_list)

            # update tabu list
            TS_reproduction.update_tabu(self, self.active_soln, self.prev_active_soln)

        chromosome_list = []
        chromosome_list.extend(TS_reproduction.neighborhood_methods(self, self.active_soln))

        return chromosome_list
    

class TS_reproduction():
    """
    Functions for performing reproduction of chromosomes using TS
     
    Written by Jake Mikouchi. 07/03/26
    """

    def initial_selection(self, pop_list):
        """
        Parses through the initial population to select starting point for TS optimization.
        To avoid bias in the starting point, a roulette style selection is utilized.

        Written by Jake Mikouchi. 07/03/26
        """

        # if initial solution is provided then use that otherwise select with roulette
        if self.input.initial_population:
            winner = pop_list[0].chromosome
        else:
            probability_sum = 0
            selection_probability = {}
            selection_probability['low_bound'] = []
            selection_probability['up_bound']  = []
            
            shift_fitness_scale = min([soln.fitness_value for soln in pop_list]) #shift fitness values to start at zero (this corrects for negative fitness values)
            for solution in pop_list:
                selection_probability['low_bound'].append(probability_sum)
                probability_sum += (solution.fitness_value - shift_fitness_scale)
                selection_probability['up_bound'].append(probability_sum)

            value = random.random()
            value = value*probability_sum
            for j, solution in enumerate(pop_list):
                if selection_probability['low_bound'][j] <= value <= selection_probability['up_bound'][j]:
                    winner = solution.chromosome

        # update best fitness 
        self.best_fitness = solution.fitness_value
        self.active_fitness = solution.fitness_value

        return winner
    
    def update_active(self, pop_list):
        """
        Updates the active solution by replacing it with the best performing solution in the set of neighbors (pop_list).
        If the improved solutions contain a tabu move then it will only replace the active solution if if passes the aspirational test.

        Written by Jake Mikouchi. 07/03/26
        """
        winner = []
        fitness_compare = -np.inf
        for soln in pop_list:
            if soln.fitness_value > fitness_compare:
                allow_tabu = TS_reproduction.aspiration_methods(self, self.active_soln, soln)
                if allow_tabu:
                    winner = soln
                    fitness_compare = soln.fitness_value

        if winner:
            # update best fitness
            if winner.fitness_value > self.best_fitness:
                self.best_fitness =  winner.fitness_value
            self.active_fitness = winner.fitness_value
            winner = winner.chromosome
        else:
            # catch for cases where no new active is found
            winner = self.active_soln

        return winner
    
    def update_tabu(self, active_soln, prev_active_soln):
        """
        Determines the inverse of the move from prev_active_soln to active_soln and adds the move to the tabu list.
        Will also cut solutions from the tabu list as it grows beyond the user specified length.

        Written by Jake Mikouchi 07/03/26
        """        

        if active_soln != prev_active_soln:
            for idx in range(len(active_soln)): 
                new_tabu = {}
                if active_soln[idx] != prev_active_soln[idx]:
                    new_tabu['position'] = idx  
                    new_tabu['move'] = prev_active_soln[idx]
            
                    self.tabu_list.append(new_tabu)

        while len(self.tabu_list) > self.input.num_tabu:
            self.tabu_list.pop(0)

    def neighborhood_methods(self, active_soln):
        """
        Distributes to neighborhood construction methods
        
        Written by Jake Mikouchi. 07/03/2026
        """

        method = self.input.neighborhood_construct
        if method["method"] == 'flip':
            pop_list = TS_perturbations.flip(self.input, self.input.population_size, method['num_flips'], active_soln) 
        
        return pop_list
    
    def aspiration_methods(self, active_soln, challenger_soln):
        """
        Distributes to aspiration methods
        
        Written by Jake Mikouchi. 07/03/2026
        """

        # checks to see if tabu moves were conducted
        if TS_aspirations.check_tabu(self, active_soln, challenger_soln):
            method = self.input.aspiration
            if method["method"] == 'improved_best':
                allow_tabu = TS_aspirations.improved_best(self.best_fitness, challenger_soln) 
            if method["method"] == 'relative_increase':
                allow_tabu = TS_aspirations.relative_increase(method, self.active_fitness, challenger_soln) 

        # returns true if no tabu changes were made
        else:
            allow_tabu = True

        return allow_tabu


class TS_perturbations():
    """
    Functions for performing perturbations of chromosomes to construct neighborhoods 
    during TS optimization.

    Written by Jake Mikouchi. 07/03/2026
    """
    def flip(input_obj, population_size, num_flips, active_soln):
        """
        performs perturbations on the active chromosome by selecting a random position in the chromosome and flipping the corresponding gene.
        The 'flipping' is performed by selecting a random valid gene.

        Written by Jake Mikouchi. 08/05/2026
        """

      ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")
        core_parameters = [input_obj.nrow, input_obj.ncol, input_obj.num_assemblies, input_obj.symmetry, input_obj.calculation_type]
        all_gene_options = input_obj.genome
        all_genes_list = list(input_obj.genome.keys())

        population = [] # neighbors
        for i in range(population_size):
            neighbor = deepcopy(active_soln) # neighbor is a chromosome
            gene_locations = [j for j in range(len(neighbor))]
            random.shuffle(gene_locations)
            chromosome_is_valid = False
            attempts = 0
            while not chromosome_is_valid:
                while neighbor == active_soln:
                    for j in range(num_flips):
                        loc_to_mutate = random.choice(gene_locations)
                        old_gene = neighbor[loc_to_mutate]
                        gene_options = optools.Gene_Validity_check.contraceptive_check(input_obj, all_genes_list, all_gene_options,
                                                                                core_parameters, active_soln, [], loc_to_mutate)
                        try:
                            if gene_options == [0,1]:
                                new_gene = random.uniform(0,1)
                            else:
                                new_gene = random.choice(gene_options)
                        except:
                            break
                        if new_gene != old_gene:
                            if input_obj.calculation_type in ["single_cycle","eq_cycle", "lattice_physics"] and all_gene_options[new_gene]['map'][loc_to_mutate] == 1:
                                neighbor[loc_to_mutate] = new_gene
                            else:
                                neighbor[loc_to_mutate] = new_gene

                        gene_locations.pop(gene_locations.index(loc_to_mutate))

                chromosome_is_valid = optools.Gene_Validity_check.abortive_check(input_obj,all_genes_list,all_gene_options,\
                                                                                core_parameters,neighbor)
                if not chromosome_is_valid:
                    attempts += 1
                    if attempts > 100000:
                        logger.error("Chromosome perturbation has failed after 100,000 attempts; the Individual will be restored. Consider relaxing the constraints on the input space.")
                        population.append(neighbor)
                        break

            if input_obj.calculation_type in ["eq_cycle"]:
                #recreate child_chromosome
                child_chromosome = []
                for i in range(len(neighbor)):
                    if neighbor[i] == active_soln[i][0]:
                        child_chromosome.append(active_soln[i])
                    else:
                        child_chromosome.append((neighbor[i],None))
                child_chromosome = optools.Solution.EQ_reload_fuel(input_obj.genome,core_parameters,child_chromosome)

            else: 
                population.append(neighbor)

        return population
    

class TS_aspirations():
    """
    class for storing aspiration methods in Tabu search

    Written by Jake Mikouchi. 07/03/26
    """
    
    def check_tabu(self, active_soln, challenger_soln):
        """
        Checks to see if any tabu moves were made

        written by. Jake Mikouchi 07/03/26
        """
        tabu_move = False
        challenge_chrome = challenger_soln.chromosome

        for idx in range(len(active_soln)):
            for tabu in self.tabu_list:
                if idx == tabu['position']:
                    # determine if variable is continous range
                    for gene_type in self.input.genome.keys():
                        if self.input.calculation_type == "numeric_variable":
                            if 'continuous_range' in self.input.gene_options[gene_type].keys():
                                if challenge_chrome[idx] >= tabu['move'] - self.input.tabu_bands and challenge_chrome[idx] <= tabu['move'] + self.input.tabu_bands:
                                    tabu_move = True
                                    break
                        else:
                            if challenge_chrome[idx] == tabu['move']:
                                tabu_move = True
                                break

        return tabu_move

    def improved_best(best_fitness, challenger_soln):
        """
        allows tabu moves to be used only if the tabu move directly improves the best fitness value.

        Written by. Jake Mikouchi 07/03/26
        """
        allowed = False
        if challenger_soln.fitness_value >= best_fitness:
            allowed = True

        return allowed
    
    def relative_increase(method, active_fit_val, challenger_soln):
        """
        allows tabu moves to be used only if the fitness of the challenger is above a user defined threshhold.

        Written by. Jake Mikouchi 07/03/26
        """
        allowed = False
        if challenger_soln.fitness_value >= active_fit_val + method['increase']:
            allowed = True

        return allowed