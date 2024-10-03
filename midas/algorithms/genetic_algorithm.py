## Import Block ##
import os
from pathlib import Path
from shutil import rmtree
import logging
from copy import deepcopy
from multiprocessing import Pool
from itertools import repeat
import random
import csv
from midas.utils import optimizer_tools as optools
from midas.applications import parcs332 as parcs


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
    """
    def __init__(self,
                 population,
                 generation,
                 fitness,
                 #!selection,
                 input):
        
        self.population = population
        self.generation = generation
        self.fitness = fitness
        #!self.selection = selection
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
    
    def generate_solution(self,name,chromosome=None):
        """
        Function for creating new solutions/indiviuals/particles in the optimization space.
        
        Written by Nicholas Rollins. 09/27/2024
        """
        #generate blank solution object
        soln = optools.Solution(f"{name}")
        
        #copy objectives and constraints from input
        soln.parameters = deepcopy(self.input.objectives)
        
        if chromosome:
            soln.chromosome = chromosome
        else: #generate random chromosome
            soln.chromosome = soln.generate_initial(self.input.genome)

        return soln
    
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

    def main(self):
        """
        #!TODO: write docstring.
        
        Written by Nicholas Rollins. 09/27/2024
        """
    ## Create results directory and/or clear old
        cwd = Path(os.getcwd())
        results_dir = cwd.joinpath(self.input.results_dir_name)
        if os.path.exists(results_dir):
            rmtree(results_dir, ignore_errors=True)
            os.mkdir(results_dir)
        else:
            os.mkdir(results_dir)
        
    ## Initialize beginning population
        logger.info("Generating initial population of %s individuals...", self.input.population_size)
        self.population.current = []
        for i in range(self.population.size):
            self.population.current.append(self.generate_solution(f'Gen_0_Indv_{i}'))
        
        pool = Pool(processes=self.input.num_procs) #initialize parallel execution
        
    ## Evaluate fitness
        logger.info("Calculating fitness for generation %s...", self.population.current)
        ## Execute and parse objective/constraint values
        self.population.current = pool.starmap(parcs.evaluate, zip(self.population.current, repeat(self.input)))
        ## Calculate fitness from objective/constriant values
        for soln in self.population.current:
            soln.fitness_value = self.fitness.calculate(soln.parameters)
        logger.info("Done!")
    
    ## Archive initial results
        archive_header = "Generation,Individual,Fitness Value"
        for param in self.input.parameters.keys():
            archive_header += ',' + str(param)
        archive_header += ",Chromosome"
        ## write output file
        with open("optimizer_results.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(archive_header)
        
        best_soln_index = [s['value'] for s in self.population.current].index(max([s['value'] for s in self.population.current]))
        for i in range(len(self.population.current)):
            soln = self.population.current[i]
            soln_result_string = str(self.generation.current) + ',' + str(i) + ',' + str(soln.fitness_value)
            for param in soln.parameters.keys():
                soln_result_string += ',' + str(soln.parameter[param]['value'])
            for gene in soln.chromosome:
                soln_result_string += ',' + str(gene)
                ## write to output file
                with open("optimizer_results.csv", 'a') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter=',')
                    csvwriter.writerow(soln_result_string)
            if i == best_soln_index:
                best_soln_string = ",".join(soln_result_string.split(',')[2:]) #remove generation and individual index.
        
        
        logger.info("Generation %s Best Individual:", self.population.current)
        logger.info(best_soln_string+'\n') #!should there be a header to label entries?
        
    ## Iterate over generations #!TODO: I believe selection is missing from this process?
        for self.generation.current in range(self.generation.total):
        ## Create new generation
            logger.info("Creating population of %s individuals for generation %s...", self.input.population_size, self.population.current+1)
            new_chromosome_list = self.reproduction(self.population.current)
            self.population.current = []
            for i in range(len(new_chromosome_list)):
                self.population.current.append(self.generate_solution(f'Gen_%s_Indv_{i}', self.population.current+1))
        
        ## Evaluate fitness
            logger.info("Calculating fitness for generation %s...", self.population.current)
            ## Execute and parse objective/constraint values
            self.population.current = pool.starmap(parcs.evaluate, zip(self.population.current, repeat(self.input)))
            ## Calculate fitness from objective/constriant values
            for soln in self.population.current:
                soln.fitness_value = self.fitness.calculate(soln.parameters)
            logger.info("Done!")
        
        ## Archive results
            best_soln_index = [s['value'] for s in self.population.current].index(max([s['value'] for s in self.population.current]))
            for i in range(len(self.population.current)):
                soln = self.population.current[i]
                soln_result_string = str(self.generation.current) + ',' + str(i) + ',' + str(soln.fitness_value)
                for param in soln.parameters.keys():
                    soln_result_string += ',' + str(soln.parameter[param]['value'])
                for gene in soln.chromosome:
                    soln_result_string += ',' + str(gene)
                    ## write to output file
                    with open("optimizer_results.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_string)
                if i == best_soln_index:
                    best_soln_string = ",".join(soln_result_string.split(',')[2:]) #remove generation and individual index.
            
            logger.info("Generation %s Best Individual:", self.population.current)
            logger.info(best_soln_string+'\n') #!should there be a header to label entries?
    
        return


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