## Import Block ##
import math
import random
from copy import deepcopy
"""
These are generic optimizer classes that are shared by all algorithms. #!TODO: can this be solved with the super().__init__ method?
"""

## Classes ##
class Population():
    """
    Container to hold the population of optimizer solutions.
    
    Updated by Nicholas Rollins. 09/24/2024
    """
    def __init__(self, pop_size, number_genes):
        if pop_size == 'calculate_from_genes':
            self.size = self.calculate_size(number_genes)
        else:
            self.size = pop_size
        
        self.current = []
        self.archive = {'solutions':[], 'fitnesses':[], 'parameters':[]}

    def calculate_size(self, number_genes):
        """
        Calculates the population size based on total number genes.
        """
        return int(10*math.sqrt(number_genes))


class Generation(): #!TODO: an object for holding two integers is silly; fold this into Population?
    """
    Container for tracking generation number of the optimizer.
    
    Updated by Nicholas Rollins. 09/24/2024
    """
    def __init__(self, num_gens, number_changes):
        if num_gens == 'calculate_from_genes':
            self.total = self.calculate_total_generations(number_changes)
        else:
            self.total = num_gens
        
        self.current = 0

    def calculate_total_generations(self, number_changes):
        """
        Calculates the total number of generations in an optimization
        based on the number of changes.
        """
        return int(5*math.sqrt(number_changes))


class Solution():
    """
    This is the generic solution class to represent solutions in the optimization.

    Parameters: None

    Written by Brian Andersen. 1/7/2019
    Updated by Nicholas Rollins. 09/11/2024
    """
    def __init__(self, name=None):
        self.name = name
        self.parameters = {}
        self.chromosome = []
        self.fitness_value = float("NaN")
    
    def generate_initial(self,calc_type,genome):
        """
        Generates the initial solutions to the optimization problem by
          randomly generating a new chromosome.

        Parameters: 
            genome: Dictionary
                The genome portion of the dictionary settings file. 

        Written by Nicholas Rollins. 10/11/2024
        """
        if calc_type in ['single_cycle']:
            return self.LP_chromosome(genome)
        elif calc_type == 'eq_cycle':
            return self.EQ_chromosome(genome)
        else:
            return None
    
    def LP_chromosome(self,genome):
        genes_list = list(genome.keys())
        chromosome_length = []
        for gene in genes_list:
            chromosome_length.append(len(genome[gene]['map']))
        
        chromosome = []
        for i in range(max(chromosome_length)):
                gene_options = Constrain_Input.calc_gene_options(genes_list, genome, chromosome)
                gene = random.choice(gene_options)
                if genome[gene]['map'][i]: #check that the selected gene option is viable at this location.
                    chromosome.append(gene)
        
        return chromosome
    
    def EQ_chromosome(self,genome):
        raise ValueError("DEBUG STOP")#!


class Constrain_Input():
    def calc_gene_options(genes_list, genome, chromosome):
        """
        Constrain the available options for the chromosome based on
        the existing inventory.
        
        Written by Nicholas Rollins. 10/10/2024
        """
        valid_genes_list = []
        for gene in genes_list:
            if genome[gene]['constraint']:
                ctype = genome[gene]['constraint']['type']
                cvalue = genome[gene]['constraint']['value']
                if ctype == 'max_quantity':
                    if chromosome.count(gene) < cvalue: #only include option if less than the max quantity have been already used.
                        valid_genes_list.append(gene)
                elif ctype == 'less_than_variable':
                    if chromosome.count(gene) < chromosome.count(cvalue): #only include option if fewer than the target option have been already used.
                        valid_genes_list.append(gene)
            else:
                valid_genes_list.append(gene)
        
        return valid_genes_list


class Fitness(object):
    """
    The generic fitness function. Requires user specified weights for every 
    design objective. Takes the Form F =   W1(minimize objectives) 
                                         - W2(maximize objectives)
                                         + W3(meet targets)
                                         + W4(satisfy objectives)
    Written by Brian Andersen. 1/18/2020
    Updated by Nicholas Rollins. 09/27/2024
    """
    def __init__(self):
        """
        Fitness object does need to be initialized with no attributes, it serves
        as a container for the Fitness.calculate() function.
        """
        pass

    def calculate(self,parameters):
        """
        Calculates the generic fitness function, based on the listed fitness
        function above. Returns the solution list with evaluated fitnesses.

        Parameters:
            parameters: Dict
                Dictionary of objectives/constraints values to be included 
                in the fitness function.

        Written by Nicholas Rollins. 09/27/2024
        """
        fitness = 0.0
        
        for param in parameters:
            pgoal = parameters[param]['goal']
            pweight = parameters[param]['weight']
            pvalue = parameters[param]['value']
            if pgoal == 'maximize':
                fitness += pvalue
            elif pgoal == 'minimize':
                fitness -= pvalue
            elif pgoal == 'meet_target':
                ptarget = parameters[param]['target']
                fitness -= abs(ptarget - pvalue)
            elif pgoal == 'greater_than_target':
                ptarget = parameters[param]['target']
                fitness += pvalue - ptarget
            elif pgoal == 'less_than_target':
                ptarget = parameters[param]['target']
                fitness += ptarget - pvalue
        
        return fitness