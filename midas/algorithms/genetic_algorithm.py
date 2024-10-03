## Import Block ##
import logging
from copy import deepcopy
import random
from midas.utils import optimizer_tools as optools


## Classes ##
class Genetic_Algorithm():
    """
    Class for performing optimization using the genetic algorithm.

    Parameters:
        population: Class
            Class that contains the population size and stores the current
            solutions in the parent and child populations.
        generation: Class
            Keeps track of the current and total number of generations that
        selection: Class _evaluate
            different solutions in the optimization.
        file_settings: Dictionary
            The settings file read into the optimization. Carried through because
            some information needed to be carried along, but pickling all of the
            information didn't seem like a good way to carrty it thorugh the optimization.

    Written by Brian Andersen. 1/9/2020
    Updated by Nicholas Rollins. 10/03/2024
    """
    def __init__(self, input):
        self.input = input
        
        #!TODO: check this functionality and add it back in if necessary.
        '''if 'cleanup' in file_settings['optimization']:
            self.number_generations_post_cleanup = 2 #Arbitrarily chosen default.
            if file_settings['optimization']['cleanup']['perform']:
                self.perform_cleanup = True
            else: 
                self.perform_cleanup = False
        else:
            self.perform_cleanup = False
        
        if 'neural_network' in file_settings:
            from crudworks import CRUD_Predictor
            self.crud = CRUD_Predictor(file_settings)'''
    
    def reproduction(self,pop_list):
        """
        #!TODO: write docstring.
        
        Updated by Nicholas Rollins. 09/27/2024
        """
    ## Select individuals for mutation. The rest undergo crossover.
        mutation_list  = []
        crossover_list = []
        for indv in pop_list:
            if random.random() < self.mutation.rate:
                mutation_list.append(indv.genome)
            else:
                crossover_list.append(indv.genome)

        if len(self.crossover_list)%2 == 1:
            swap_solution = random.choice(mutation_list)
            crossover_list.append(swap_solution)
            mutation_list.remove(swap_solution)
    
    ## Perform Crossover
        crossover_mates_lists = GA_reproduction.crossover_assign_mates(crossover_list)
        child_chromosome_list = []
        for mate_one, mate_two in zip(crossover_mates_lists[0], crossover_mates_lists[1]):
            child_one, child_two = GA_reproduction.crossover(mate_one, mate_two, self.mutation.rate)
            child_chromosome_list.extend([child_one, child_two])
    
    ## Perform Mutation
        for chromosome in self.mutation_list:
            child = mutate_by_chromosome.mutate_by_chromosome(chromosome)
            child_chromosome_list.append(child)
        
        return child_chromosome_list


class GA_reproduction():
    """
    #!TODO: write docstring.
    
    Written by Nicholas Rollins. 09/27/2024
    """
    def crossover_assign_mates(self, crossover_list):
        """
        #!TODO: write docstring.
        
        Written by Nicholas Rollins. 09/27/2024
        """
        first_mate_list  = []
        second_mate_list = []
        while len(crossover_list) > 0:
            parent_one = random.choice(crossover_list)

            partner_similarities = -1
            for suitor in crossover_list:
                if suitor is parent_one:
                    pass
                else:
                    new_suitor_similarities = self.check_intersection(parent_one, suitor)
                    if new_suitor_similarities > partner_similarities:
                        parent_two = suitor[:]
                        partner_similarities = new_suitor_similarities

            crossover_list.remove(parent_one)
            if parent_two in crossover_list:
                crossover_list.remove(parent_two)

            first_mate_list.append(parent_one)
            second_mate_list.append(parent_two)
        return (first_mate_list, second_mate_list)
    
    def crossover(self, chromosome_one, chromosome_two, mutation_rate):
        """
        Function for performing crossover of the mated solutions by swapping
        differing genes between the chromosomes. Edited as of 1/20/2020 to account
        for fixed gene counts in problems.

        Parameters:
            chromosome_one: list
                The first genome that is to undergo crossover.
            chromosome_two: list
                The second geneome that is to undergo crossover.
            mutation_rate: float
                The percent of solutions that undergo mutation. Used here to
                determine the number of genes that are swapped.

        Written by Brian Andersen. 8/29/2019
        """
        difference_positions = self.return_different_positions(chromosome_one,
                                                               chromosome_two)
        child_one = []
        child_two = []
        position_count = 0
        for i, j in zip(chromosome_one, chromosome_two):
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
    
## Mutation types ##
    def mutate_by_chromosome(self, chromosome):
        """
        Generates new solution through mutation in the optimization.
        
        Updated by Nicholas Rollins. 09/27/2024
        """

        child_chromosome = deepcopy(chromosome)

        while child_chromosome == chromosome:
            for i in range(self.number):
                mutate = random.randint(0, len(child_chromosome)-1)
                old_gene = child_chromosome[mutate]

                key_list = list(self.genome_map.keys())
                new_gene = random.choice(key_list)
                if new_gene == old_gene:
                    pass
                else:
                    if self.genome_map[new_gene][mutate] == 1:
                        child_chromosome[mutate] = new_gene

        return child_chromosome