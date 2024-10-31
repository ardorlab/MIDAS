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
    
    def reproduction(self, pop_list): #!TODO: apply constraints.
        """
        Generates a new generation of individuals by performing crossover and mutation 
        operations on the parents generation. These operations are performed on an ordered 
        list of parents, as biased by the selection method.
        
        Updated by Nicholas Rollins. 09/27/2024
        """
    ## Perform selection process of parents
        # Selection is established by biasing reproduction using an ordered list.
        pop_list = self.selection_methods(pop_list, self.input.selection)

    ## Container for holding new list of child chromosomes
        child_chromosome_list = []
        # preserve best performing parent in the next generation
        child_chromosome_list.append(max(pop_list, key=lambda soln: soln.fitness_value)
    
    ## Select individuals for mutation. The rest undergo crossover.
        mutation_list  = []
        crossover_list = []
        for soln in pop_list:
            if random.random() < self.input.mutation_rate['initial_rate']: #!TODO: reapply variable mutation rates
                mutation_list.append(soln.chromosome)
            else:
                crossover_list.append(soln.chromosome)

        if len(crossover_list)%2 == 1: #make sure there is an even number of individuals in the crossover list.
            try:
                soln_to_move = random.choice(mutation_list)
                crossover_list.append(soln_to_move)
                mutation_list.remove(soln_to_move)
            except IndexError: # mutation_list could be empty.
                soln_to_move = random.choice(crossover_list)
                mutation_list.append(soln_to_move)
                crossover_list.remove(soln_to_move)
    
        LWR_core_parameters = [self.input.nrow, self.input.ncol, self.input.symmetry]
    ## Perform Crossover
        crossover_mates_lists = GA_reproduction.crossover_assign_mates(crossover_list)
        for mate_one, mate_two in zip(crossover_mates_lists[0], crossover_mates_lists[1]):
            child_one, child_two = GA_reproduction.crossover(mate_one, mate_two, 
                                                            self.input.mutation_rate['initial_rate'], 
                                                            LWR_core_parameters, self.input.genome, self.input.batches) # default for 'batches' is None. #!TODO: this should be an adaptive mutation rate.
            child_chromosome_list.extend([child_one, child_two])
    
    ## Perform Mutation
        for chromosome in mutation_list:
            if self.input.mutation_type == "mutate_by_gene":
                child = GA_reproduction.mutate_by_gene(self.input, chromosome)
            else:
                raise ValueError("Requested mutation type not recognized.")
            child_chromosome_list.append(child)
    
    ## Check new population size
        while len(child_chromosome_list) > len(pop_list):
            del child_chromosome_list[random.randint(1,len(child_chromosome_list))] # remove any random child to reduce population size (except first child)
        
        return child_chromosome_list
    
    def selection_methods(self, pop_list, method):
        """
        Method for distributing to the requested GA selection method.
        
        Written by Nicholas Rollins. 10/08/2024
        """
        if method == 'roulette':
            pop_list = GA_selection.roulette(pop_list, len(pop_list)) #assume unchanging population size 
        elif method == 'tournament':
            pop_list = GA_selection.tournament(pop_list, len(pop_list)) #assume unchanging population size 
        return pop_list


class GA_reproduction():
    """
    Functions for performing reproduction of chromosomes using GA
    methodologies, including crossover and mutation. Does not include
    the random generation of new individuals.
    
    Written by Nicholas Rollins. 09/27/2024
    """
    def crossover_assign_mates(crossover_list,SS_decode=False):
        """
        Before crossover between chromosomes can be performed, the
        individuals in the population need to be paired off.
        
        Written by Nicholas Rollins. 09/27/2024
        """
        first_mate_list  = []
        second_mate_list = []
        while len(crossover_list) > 0:
            parent_one = random.choice(crossover_list)

            partner_similarities = -1
            for suitor in crossover_list:
                if not suitor is parent_one:
                    if SS_decode:
                        parent_one_zones = [loc[0] for loc in parent_one]
                        suitor_zones = [loc[0] for loc in suitor]
                        new_suitor_similarities = GA_reproduction.check_intersection(parent_one_zones, suitor_zones)
                    else:
                        new_suitor_similarities = GA_reproduction.check_intersection(parent_one, suitor)
                    if new_suitor_similarities > partner_similarities:
                        parent_two = suitor[:]
                        partner_similarities = new_suitor_similarities

            crossover_list.remove(parent_one)
            if parent_two in crossover_list:
                crossover_list.remove(parent_two)

            first_mate_list.append(parent_one)
            second_mate_list.append(parent_two)
        
        return (first_mate_list, second_mate_list)
    
    def crossover(chromosome_one, chromosome_two, mutation_rate, LWR_core_parameters, genome, batches=None):
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
        difference_positions = GA_reproduction.return_different_positions(chromosome_one, chromosome_two)
        genes_list = list(genome.keys())
        if batches:
            genes_list = list(batches.keys())
        
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            child_one = []
            child_two = []
            position_count = 0
            for i, j in zip(chromosome_one, chromosome_two):
                if position_count in difference_positions:
                    if random.random() < mutation_rate:
                        if batches:
                            if i[0] != j[0]: #batches don't match.
                                #constrain input for child_one and child_two
                                child1_zone = [loc[0] for loc in child_one+chromosome_one[len(child_one):]]
                                child2_zone = [loc[0] for loc in child_two+chromosome_two[len(child_two):]]
                                child1_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, batches, 
                                                                                             LWR_core_parameters, child1_zone)
                                child2_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, batches, 
                                                                                             LWR_core_parameters, child2_zone)
                                if j[0] in child1_gene_opts and i[0] in child2_gene_opts: #swap batches
                                    child_one.append((j[0],None))
                                    child_two.append((i[0],None))
                                else: #don't swap batches
                                    child_one.append(i)
                                    child_two.append(j)
                            elif i[1] != j[1]: #FA's don't match.
                                #swapping FA's won't be constrained here; any potential conflicts will be resolved by EQ_reload_fuel().
                                child_one.append(j)
                                child_two.append(i)
                        else:
                            #constrain input for child_one and child_two
                            child1_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, genome, LWR_core_parameters, 
                                                                                         child_one+chromosome_one[len(child_one):])
                            child2_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, genome, LWR_core_parameters, 
                                                                                         child_two+chromosome_two[len(child_two):])
                            if j in child1_gene_opts and i in child2_gene_opts: #swap genes
                                child_one.append(j)
                                child_two.append(i)
                            else: #don't swap genes
                                child_one.append(i)
                                child_two.append(j)
                    else: #don't swap genes
                        child_one.append(i)
                        child_two.append(j)
                else: #don't swap genes
                    child_one.append(i)
                    child_two.append(j)
                position_count += 1
                
            attempts += 1
            if optools.Constrain_Input.check_constraints(genes_list,batches,LWR_core_parameters,\
                                                         [loc[0] for loc in child_one]) and \
               optools.Constrain_Input.check_constraints(genes_list,batches,LWR_core_parameters,\
                                                         [loc[0] for loc in child_two]):
                chromosome_is_valid = True
            if attempts > 1000:
                raise ValueError("Crossover has failed after 1,000 attempts. Consider relaxing the constraints on the input space.")
                
        if batches: #reload fuel in 'None' locations.
            child_one = optools.Constrain_Input.EQ_reload_fuel(genome,LWR_core_parameters,child_one)
            child_two = optools.Constrain_Input.EQ_reload_fuel(genome,LWR_core_parameters,child_two)

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
            if not i == j:
                diff_position_list.append(position_count)
            position_count += 1
        return diff_position_list
    
## Mutation types ##
    def mutate_by_gene(input_obj, chromosome):
        """
        Generates a new solution by randomly mutating a single gene.
        
        Updated by Nicholas Rollins. 09/27/2024
        """
        ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")
        
        LWR_core_parameters = [input_obj.nrow, input_obj.ncol, input_obj.symmetry]
        
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
            while new_soln == old_soln:
                for i in range(num_mutations):
                    loc_to_mutate = random.randint(0, len(new_soln)-1) #choose a random gene
                    old_gene = new_soln[loc_to_mutate]
                    gene_options = optools.Constrain_Input.calc_gene_options(all_genes_list, all_gene_options, LWR_core_parameters, old_soln) #constraint input
                    new_gene = random.choice(gene_options)
                    if new_gene != old_gene:
                        if all_gene_options[new_gene]['map'][loc_to_mutate] == 1:
                            new_soln[loc_to_mutate] = new_gene
            chromosome_is_valid = optools.Constrain_Input.check_constraints(all_genes_list,all_gene_options,\
                                                                            LWR_core_parameters,new_soln)
            attempts += 1
            if attempts > 100000:
                logger.error(f"Mutate-by-Gene has failed after 100,000 attempts; Parent will be restored. Consider relaxing the constraints on the input space.")
                new_soln = deepcopy(old_soln)

        if input_obj.calculation_type in ["eq_cycle"]:
            #recreate child_chromosome
            child_chromosome = []
            for i in range(len(new_soln)):
                if new_soln[i] == chromosome[i][0]:
                    child_chromosome.append(chromosome[i])
                else:
                    child_chromosome.append((new_soln[i],None))
            child_chromosome = optools.Constrain_Input.EQ_reload_fuel(input_obj.genome,LWR_core_parameters,child_chromosome)

        return child_chromosome


class GA_selection():
    """
    #!TODO: write docstring.
    
    Written by Nicholas Rollins. 10/08/2024
    """
    def roulette(pop_list, desired_pop_size):
        """
        Performs a roulette selection of the solutions. Can be used for determining
        solution front or parents for the next generation or whatever.
        
        Updated by Nicholas Rollins. 10/08/2024
        """
        unused_solutions = deepcopy(pop_list)
        winners = []
        for i in range(desired_pop_size):
            probability_sum = 0
            selection_probability = {}
            selection_probability['low_bound'] = []
            selection_probability['up_bound']  = []
            
            # the probability of selecting an individual is weighted by its fitness
            for solution in unused_solutions:
                selection_probability['low_bound'].append(probability_sum)
                probability_sum += solution.fitness
                selection_probability['up_bound'].append(probability_sum)

            value = random.random()
            value = value*probability_sum
            for j, solution in enumerate(unused_solutions):
                if selection_probability['low_bound'][j] <= value <= selection_probability['up_bound'][j]:
                    winners.append(solution)
                    unused_solutions.remove(solution)

            # if we run out of parents, continue from a fresh list
            if not unused_solutions:
                unused_solutions = deepcopy(pop_list)

        return winners
    
    def tournament(pop_list, desired_pop_size):
        """
        Performs a tournament selection of the solutions. Each tournament selects up to half of the 
        list as winners. If additional solutions are needed, a second tournament is held with the losers
        from the first tournament.
        
        Parameters:
            pop_list: list
                The solutions that will compete in the tournament
            desired_pop_size: int
                The number of solutions chosen at the end of the tournament.

        Written by Brian Andersen. 1/9/2020
        Updated by Nicholas Rollins. 10/08/2024
        """
        unused_solutions = pop_list[:]
        used_solutions = []
        winners = []
        for i in range(desired_pop_size):
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
                    unused_solutions = deepcopy(pop_list)
                    used_solutions = []

            unused_solutions.remove(one)
            unused_solutions.remove(two)

        return winners