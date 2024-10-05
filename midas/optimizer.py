## Import Block ##
import os
from pathlib import Path
from shutil import rmtree
from copy import deepcopy
from multiprocessing import Pool
from itertools import repeat
import csv

from midas.algorithms import genetic_algorithm as GA
from midas.utils import optimizer_tools as optools
from midas.codes import parcs332 as parcs
#!from midas.applications import parcs_332


## Classes ##
class Optimizer():
    """
    #!TODO: write docstring.
    
    Written by Nicholas Rollins. 09/23/2024
    """
    def __init__(self, inp_lines):
        self.input = inp_lines
    
    def build_optimizer(self):
        """
        Assembles all the pieces of the optimization and create the
        algorithm object.
        
        Written by Brian Andersen. 01/07/2019
        SA added by Johnny Klemes. 03/22/2020
        RL added by G. K. Delipe.  03/24/2023
        Updated by Nicholas Rollins. 09/11/2024
        """
        methodology = self.input.methodology
        num_gene_combos = self.calculate_number_gene_combinations(self.input.genome)
        
        self.population = optools.Population(self.input.population_size, num_gene_combos)
        self.generation = optools.Generation(self.input.num_generations, num_gene_combos)
        self.fitness    = optools.Fitness()
        if self.input.calculation_type == "parcs":
            self.eval_func  = parcs.evaluate #assign, don't execute.
        
        if methodology == 'genetic_algorithm':
            self.algorithm = GA.Genetic_Algorithm(self.input)
        #!TODO: Add the other algorithms back in.

        return
    
    def calculate_number_gene_combinations(self, genome_map):
        """
        Calculates number of possible gene combinations for inputs.

        parameters:
            genome_map: From the input file data dictionary would be the two
            keys: [decision_variables][parameters]. 
        
        Written by Brian Andersen. 1/7/2019
        Updated by Nicholas Rollins. 09/24/2024
        """
        
        number_changes = 0
        for chromosome in genome_map:
            if chromosome == 'symmetry_list':
                pass
            else:
                number_changes += len(genome_map[chromosome]['map'])  

        return number_changes
    
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
        self.population.current = pool.starmap(self.eval_func, zip(self.population.current, repeat(self.input)))
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
            new_chromosome_list = self.algorithm.reproduction(self.population.current)
            self.population.current = []
            for i in range(len(new_chromosome_list)):
                self.population.current.append(self.generate_solution(f'Gen_%s_Indv_{i}', self.population.current+1))
        
        ## Evaluate fitness
            logger.info("Calculating fitness for generation %s...", self.population.current)
            ## Execute and parse objective/constraint values
            self.population.current = pool.starmap(self.eval_func, zip(self.population.current, repeat(self.input)))
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