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
        self.archive = []

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
    
    def generate_initial(self,genome):
        """
        Generates the initial solutions to the optimization problem by
          randomly generating a new chromosome.

        Parameters: 
            genome: Dictionary
                The genome portion of the dictionary settings file. 

        Written by Brian Andersen. 1/9/2020
        Updated by Nicholas Rollins. 09/27/2024
        """
        chromosome_list = list(genome.keys())
        chromosome_length = []
        for gene in chromosome_list:
            chromosome_length.append(len(genome[gene]['map']))
        
        chromosome = []
        for i in range(max(chromosome_length)):
                gene = random.choice(chromosome_list)
                if genome[gene]['map'][i]:
                    chromosome.append(gene)
        
        return chromosome


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

    def calculate(self,solution):
        """
        Calculates the generic fitness function, based on the listed fitness
        function above. Returns the solution list with evaluated fitnesses.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Nicholas Rollins. 09/27/2024
        """
        solution.fitness = 0.0
        
        for param in solution.parameters:
            pgoal = solution.parameters[param]['goal']
            pweight = solution.parameters[param]['weight']
            pvalue = solution.parameters[param]['value']
            if pgoal == 'maximize':
                solution.fitness += pvalue
            elif pgoal == 'minimize':
                solution.fitness -= pvalue
            elif pgoal == 'meet_target':
                solution.fitness -= abs(ptarget - pvalue)
            elif pgoal == 'greater_than_target':
                solution.fitness += pvalue - ptarget
            elif pgoal == 'less_than_target':
                solution.fitness += ptarget - pvalue
        
        return solution.fitness