## Import Block ##
from midas.algorithms import genetic_algorithm as GA
from midas.utils import optimizer_tools as optools
#!from midas.applications import parcs_332


## Classes ##
class Optimizer_Factory():
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
        
        if methodology == 'genetic_algorithm':
            population_    = optools.Population(self.input.population_size, num_gene_combos)
            generation_    = optools.Generation(self.input.num_generations, num_gene_combos)
            fitness_       = optools.Fitness()
            self.optimization = GA.Genetic_Algorithm(population   = population_,
                                                     generation   = generation_,
                                                     fitness      = fitness_,
                                                     input       = self.input)
        #!TODO: Add the other algorithms back in.

        return self.optimization
    
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