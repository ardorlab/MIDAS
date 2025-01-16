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
    
    def reproduction(self, pop_list, current_generation): #!TODO: apply constraints.
        """
        Generates a new generation of individuals by performing crossover and mutation 
        operations on the parents generation. These operations are performed on an ordered 
        list of parents, as biased by the selection method.
        
        Updated by Nicholas Rollins. 09/27/2024
        updatde by Jake Mikouchi. 01/06/2025
        """
    ## Container for holding new list of child chromosomes
        child_chromosome_list = []
   
    ## create list of solutions for crossover.
        crossover_list = []
        for soln in pop_list:
            crossover_list.append(soln.chromosome)

    # determine elites from previous generation if requested in input file
        Elites = GA_reproduction.Determine_Elites(self, pop_list)

        if len(crossover_list)%2 == 1: #make sure there is an even number of individuals in the crossover list.
            soln_to_move = random.choice(crossover_list)
            crossover_list.append(soln_to_move)
    
    ## preserve core parameters 
        LWR_core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]

    ## TODO add selection method for quelling

    ## Perform Crossover
        Num_children = int(self.input.population_size) - len(Elites)
        while len(child_chromosome_list) < Num_children:
            parents = self.selection_methods(pop_list, self.input.selection)
            mate_one = parents[0]
            mate_two = parents[1]
            child_one, child_two = self.crossover_methods(mate_one, mate_two, 
                                                            self.input.crossover, 
                                                            LWR_core_parameters, self.input.genome, self.input.batches) 
            if len(child_chromosome_list) == Num_children - 1: 
                child_chromosome_list.append(random.choice([child_one,child_two]))
            else:
                child_chromosome_list.extend([child_one, child_two])

    ## Perform Mutation
        curr_mutation_rate = GA_reproduction.linear_update(self.input.mutation_rate['initial_rate'], self.input.mutation_rate['final_rate'], 
                                                           current_generation, self.input.num_generations)
        mutation_list  = [] 
        for soln in child_chromosome_list:
            if random.random() < curr_mutation_rate:
                mutation_list.append(soln)
                child_chromosome_list.remove(soln)

        for chromosome in mutation_list:
            if self.input.mutation_type == "mutate_by_gene":
                child = GA_reproduction.mutate_by_gene(self.input, chromosome)
            else:
                raise ValueError("Requested mutation type not recognized.")
            child_chromosome_list.append(child)
    
    ## add elites into population
        child_chromosome_list.extend(Elites)

        return child_chromosome_list
    
    def selection_methods(self, pop_list, method):
        """
        Method for distributing to the requested GA selection method.
        
        Written by Nicholas Rollins. 10/08/2024
        updated by Jake Mikouchi. 12/24/2024
        """

        ## TODO modify selection methods for use in queling
        if method["method"] == 'roulette':
            pop_list = GA_selection.roulette(pop_list, 2) #assume no quelling (only selecting parents) 
        elif method["method"] == 'tournament':
            pop_list = GA_selection.tournament(pop_list, 2) #assume no quelling (only selecting parents) 
        elif method["method"] == 'ktournament':
            pop_list = GA_selection.ktournament(pop_list, 2, method) #assume no quelling (only selecting parents) 
        elif method["method"] == 'truncation':
            pop_list = GA_selection.truncation(pop_list, 2) #assume no quelling (only selecting parents) 
        elif method["method"] == 'sus':
            pop_list = GA_selection.sus(pop_list, 2) #assume no quelling (only selecting parents) 
        elif method["method"] == 'random':
            pop_list = GA_selection.random(pop_list, 2) #assume no quelling (only selecting parents) 
        return pop_list
    
    def crossover_methods(self, mate_one, mate_two, crossover, LWR_core_parameters, genome, batches):
        """
        Method for distributing to the requested GA crossover method.
        
        Written by Jake Mikouchi 1/5/2025
        """
        if crossover['method'] == 'uniform':
            child_one, child_two = GA_reproduction.uniform_crossover(mate_one, mate_two, 
                                                            crossover['crossover_rate'], 
                                                            LWR_core_parameters, genome, batches)
        elif crossover['method'] == 'random_element':
            child_one, child_two = GA_reproduction.random_element_crossover(mate_one, mate_two, 
                                                            crossover['num_swaps'], 
                                                            LWR_core_parameters, genome, batches)
        elif crossover['method'] == 'one_point':
            child_one, child_two = GA_reproduction.one_point_crossover(mate_one, mate_two, 
                                                            LWR_core_parameters, genome, batches)
        elif crossover['method'] == 'two_point':
            child_one, child_two = GA_reproduction.two_point_crossover(mate_one, mate_two, 
                                                            LWR_core_parameters, genome, batches)

        return child_one, child_two
    
class GA_reproduction():
    """
    Functions for performing reproduction of chromosomes using GA
    methodologies, including crossover and mutation. Does not include
    the random generation of new individuals.
    
    Written by Nicholas Rollins. 09/27/2024
    """
    
    def uniform_crossover(chromosome_one, chromosome_two, crossover_rate, LWR_core_parameters, genome, batches=None):
        """
        Function for performing crossover of the mated solutions by swapping
        differing genes between the chromosomes. Edited as of 1/20/2020 to account
        for fixed gene counts in problems.

        Parameters:
            chromosome_one: list
                The first genome that is to undergo crossover.
            chromosome_two: list
                The second geneome that is to undergo crossover.
            crossover_rate: float
                Used here to determine the number of genes that are swapped.

        Written by Brian Andersen. 8/29/2019
        Updated by Jake Mikouchi 12/24/2024
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
                    if random.random() < crossover_rate:
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
    

    def random_element_crossover(chromosome_one, chromosome_two, num_swaps, LWR_core_parameters, genome, batches=None):
        """
        genes are randomly selected within each chromosome to be swapped. 
        This is intended to have greater randomness than uniform crossover and be better for problems with fewer number of genes (assemblies)
        This method is not exactly "crossover", and is more a combination of crossover and mutation.
        This is because the position of swapped genes is not preserved in the chromosomes.
        Despite this, this method works quite well for loading patterns.

        Written by Jake Mikouci. 1/05/25
        """

        genes_list = list(genome.keys())
        if batches:
            genes_list = list(batches.keys())

        # creates list of potential positions to be swapped
        chromosome_elements = [i for i in range(len(chromosome_one))]
        # list to store genes which have already been swapped
        c1_swapped_elements = []
        c2_swapped_elements = []

        child_one = deepcopy(chromosome_one)
        child_two = deepcopy(chromosome_two)
        for i in range(num_swaps):

            c1_gene = float("Nan")
            c2_gene = float("Nan")
            child1_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, genome, LWR_core_parameters, 
                                                            child_one+chromosome_one[len(child_one):])
            child2_gene_opts = optools.Constrain_Input.calc_gene_options(genes_list, genome, LWR_core_parameters, 
                                                            child_two+chromosome_two[len(child_two):])

            attempts = 0
            chromosome_is_valid = False
            # selects genes to be swapped and ensures that children are valid solutions
            while chromosome_is_valid == False:
                c1_gene_position = random.choice(chromosome_elements)
                while c1_gene_position in c1_swapped_elements:
                    c1_gene_position = random.choice(chromosome_elements)
                c1_gene = child_one[c1_gene_position]

                c2_gene_position = random.choice(chromosome_elements)
                while c2_gene_position in c2_swapped_elements:
                    c2_gene_position = random.choice(chromosome_elements)
                c2_gene = child_two[c2_gene_position]

                attempts += 1
                if optools.Constrain_Input.check_constraints(genes_list,batches,LWR_core_parameters,\
                                                         [loc[0] for loc in child_one]) and \
                   optools.Constrain_Input.check_constraints(genes_list,batches,LWR_core_parameters,\
                                                         [loc[0] for loc in child_two]):
                    chromosome_is_valid = True
                if attempts > 1000:
                    raise ValueError("Crossover has failed after 1,000 attempts. Consider relaxing the constraints on the input space.")
                

            # swaps genes
            child_one[c1_gene_position] = c2_gene
            child_two[c2_gene_position] = c1_gene

            # stores positions that have been previously swapped
            c1_swapped_elements.append(c1_gene_position)
            c2_swapped_elements.append(c2_gene_position)
            
        return child_one, child_two

    def one_point_crossover(chromosome_one, chromosome_two, LWR_core_parameters, genome, batches=None):
        """
        the entire gene sequence of chromosome_one and chromosome_two is split and grafted at a single random point
        This method preserves gene positions within the chromosome, however this is not ideal for loading patterns.
        This is due to the spatial dependence of assemblies in loading patterns. 

        Written by Jake Mikouchi. 1/05/25
        """
        genes_list = list(genome.keys())
        if batches:
            genes_list = list(batches.keys())

        # creates list of potential positions to be swapped
        chromosome_elements = [i for i in range(len(chromosome_one))]
        # selects position to be swapped
        crossover_position = random.choice(chromosome_elements)

        child_one_seq_a = chromosome_one[:crossover_position]
        child_one_seq_b = chromosome_one[crossover_position:]
        child_two_seq_a = chromosome_two[:crossover_position]
        child_two_seq_b = chromosome_two[crossover_position:]

        # creates children with swapped genes
        child_one = child_one_seq_a + child_two_seq_b
        child_two = child_two_seq_a + child_one_seq_b

        return child_one, child_two

    def two_point_crossover(chromosome_one, chromosome_two, LWR_core_parameters, genome, batches=None):
        """
        the entire gene sequence of chromosome_one and chromosome_two is split and grafted at a single random point on each sequence

        Written by Jake Mikouchi. 1/05/25
        """
        genes_list = list(genome.keys())
        if batches:
            genes_list = list(batches.keys())

        # creates list of potential positions to be swapped
        chromosome_elements = [i for i in range(len(chromosome_one))]
        # selects both position to be swapped
        crossover_position_1 = 0
        crossover_position_2 = 0
        while crossover_position_1 == crossover_position_2:
            crossover_position_1 = random.choice(chromosome_elements)
            crossover_position_2 = random.choice(chromosome_elements)

        crossover_positions = sorted([crossover_position_1, crossover_position_2])

        child_one_seq_a = chromosome_one[:crossover_positions[0]]
        child_one_seq_b = chromosome_one[crossover_positions[0]:crossover_positions[1]]
        child_one_seq_c = chromosome_one[crossover_positions[1]:]

        child_two_seq_a = chromosome_two[:crossover_positions[0]]
        child_two_seq_b = chromosome_two[crossover_positions[0]:crossover_positions[1]]
        child_two_seq_c = chromosome_two[crossover_positions[1]:]

        # creates children with swapped genes
        child_one = child_one_seq_a + child_two_seq_b + child_one_seq_c
        child_two = child_two_seq_a + child_one_seq_b + child_two_seq_c

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
            new_soln = deepcopy(old_soln) #in the case of abortion, start from scratch.
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
                logger.error("Mutate-by-Gene has failed after 100,000 attempts; Parent will be restored. Consider relaxing the constraints on the input space.")
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

    def linear_update(initial_rate, final_rate, current_generation, num_generations):
        """
        linearly updates the mutation rate at a given generation
        
        written by Jake Mikouchi. 12/21/2024
        """
        curr_rate = initial_rate + ((final_rate - initial_rate) / num_generations) * (current_generation + 1)

        return curr_rate
    
    def Determine_Elites(self, pop_list):
        """
        Determines the best solutions (elites) of the previous generation and stores them
        The elites are used to ensure that "good genes" are maintained in the population

        written by Jake Mikouchi 12/31/24 
        """
        if self.input.elites < 1.0:
            Num_Elites = round(self.input.elites * self.input.population_size)
        elif self.input.elites > 1.0:
            Num_Elites = int(self.input.elites)

        Elites = []

        if Num_Elites > 0:
            pop_list.sort(key=lambda x: x.fitness_value, reverse=True)
            for i in pop_list:
                if i.chromosome not in Elites:
                    Elites.append(i.chromosome)
                if len(Elites) >= Num_Elites:
                    break

        return Elites


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
                    winners.append(solution.chromosome)
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
            if one.fitness_value > two.fitness_value:
                winners.append(one.chromosome)
                used_solutions.append(two)
                unused_solutions.remove(one)
            else:
                winners.append(two.chromosome)
                used_solutions.append(one)
                unused_solutions.remove(two)

            if len(unused_solutions) < 2:
                if len(used_solutions) > 2:
                    unused_solutions = used_solutions
                    used_solutions = []
                else:
                    unused_solutions = deepcopy(pop_list)
                    used_solutions = []

        return winners
    
    def ktournament(pop_list, desired_pop_size, method):
        """
        k-way tournament selection
        places k random solutions in a tournament. The solution with the highest fitness is added to the winners.
        The tournament size can be anything as long as its < len(pop_list)
        Jake Mikouchi 12/27/24
        """
        kway = method['k']
        unused_solutions = deepcopy(pop_list)
        winners = []

        for i in range(desired_pop_size):
            ksols = random.choices(unused_solutions, k=kway)
            ksolsfit = [i.fitness_value for i in ksols]

            winners.append(ksols[ksolsfit.index(max(ksolsfit))].chromosome)
            unused_solutions.remove(ksols[ksolsfit.index(max(ksolsfit))])

            if len(unused_solutions) < 2:
                if len(winners) > 2:
                    unused_solutions = used_solutions
                    used_solutions = []
                else:
                    unused_solutions = deepcopy(pop_list)
                    used_solutions = []

        return winners
    
    
    def truncation(pop_list, desired_pop_size):
        """
        Truncation Selection
        orders the population list by fitness in descending order. The top performers are then chosen as the winner.
        This method is rarely used due to the lack of randomness, and because of this it tends to perform worse than other selection methods.
        Truncation selection is sometimes refered to as "elitism selection"
        Jake Mikouchi 5/10/24
        """
        unused_solutions = deepcopy(pop_list)
        winners = []
        unused_solutions_fitness = [i.fitness_value for i in unused_solutions]

        for i in range(0, len(unused_solutions_fitness)):
            for j in range(i+1, len(unused_solutions_fitness)):
                if unused_solutions_fitness[i] <= unused_solutions_fitness[j]:
                    unused_solutions_fitness[i], unused_solutions_fitness[j] = unused_solutions_fitness[j],unused_solutions_fitness[i]
                    unused_solutions[i], unused_solutions[j] = unused_solutions[j],unused_solutions[i]

        for i in unused_solutions:
            if i.chromosome not in winners:
                winners.append(i.chromosome)
            if len(winners) >= desired_pop_size:
                break

        return winners
    
    def sus(pop_list, desired_pop_size):
        """
        stochastic universal sampling (SUS)
        SUS uses a single random value to sample all of the solutions by choosing them at evenly spaced intervals.
        Jake Mikouchi 5/10/24
        """
        unused_solutions = pop_list[:]
        winners = []

        # orders all solutions in by fitness
        unused_solutions_fitness = [i.fitness_value for i in unused_solutions]
        for i in range(0, len(unused_solutions_fitness)):
            for j in range(i+1, len(unused_solutions_fitness)):
                if unused_solutions_fitness[i] <= unused_solutions_fitness[j]:
                    unused_solutions_fitness[i], unused_solutions_fitness[j] = unused_solutions_fitness[j],unused_solutions_fitness[i]
                    unused_solutions[i], unused_solutions[j] = unused_solutions[j],unused_solutions[i]


        while len(winners) < desired_pop_size:
            # calculates average fitness of all solutions
            average_unused_solutions_fitness = sum(unused_solutions_fitness) / len(unused_solutions_fitness)
            # calculates alpha and delta parameters
            alpha = random.uniform(0,1)
            delta = average_unused_solutions_fitness * alpha
            Sum = unused_solutions_fitness[0]

            for i in range(len(unused_solutions)):
                if delta < Sum:
                    winners.append(unused_solutions[i].chromosome)
                    unused_solutions.remove(unused_solutions[i])
                    unused_solutions_fitness.remove(unused_solutions_fitness[i])
                    delta = delta + average_unused_solutions_fitness
                    break

                else: 
                    Sum += unused_solutions_fitness[i]

        return winners

    def random(pop_list, desired_pop_size):
        """
        Randomly selects the next set of parents
        Written by Jake Mikouchi 12/27/24
        """

        unused_solutions = deepcopy(pop_list)
        winners = []
        for i in range(desired_pop_size):
            random_index = random.randint(0, len(unused_solutions)-1)
            winners.append(unused_solutions[random_index].chromosome)
            unused_solutions.pop(random_index)

        return winners