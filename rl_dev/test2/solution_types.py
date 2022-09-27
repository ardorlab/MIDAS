import random
import numpy

"""
File for containing the base solution type supported by the modular optimization 
framework. A seperate file is used because the original implementation caused 
lots of problems that I didn't understand. My best guess at the moment is a weird 
recursive loop that python didn't really want to handle. So independent file
here. 

Written by Brian Andersen. 1/21/2021
"""
class Solution(object):
    """
    This is the generic solution class to represent solutions in the optimization. 
    This specific solution isn't really meant to be implemented. It serves as a 
    parent for the other specific solution classes.

    Parameters: None

    Written by Brian Andersen. 1/7/2019
    """
    def __init__(self):
        self.name = None
        self.parameters = {}
        self.genome = []
        self.fitness = None
        self.neural_network = {}
        self.fixed_genome = False

    def add_additional_information(self,settings):
        pass

    def generate_initial(self,chromosome_map):
        """Generate initial genome for assembly optimization."""
        chromosome_length = None
        chromosome_list = list(chromosome_map.keys())

        for chromosome in chromosome_map:
            if chromosome_length is None:
                chromosome_length = len(chromosome_map[chromosome]['map'])
            elif len(chromosome_map[chromosome]['map']) == chromosome_length:
                pass
            else:
                raise ValueError("Chromosome Maps are of unequal length")
                
        for i in range(chromosome_length):
            no_gene_found = True
            while no_gene_found:
                gene = random.choice(chromosome_list)
                if chromosome_map[gene]['map'][i]:
                    self.genome.append(gene)
                    no_gene_found = False

    def generate_initial_uniform(self,chromosome_map):
        chromosome_length = None
        chromosome_list = list(chromosome_map.keys())

        for chromosome in chromosome_map:
            if chromosome_length is None:
                chromosome_length = len(chromosome_map[chromosome]['map'])
            elif len(chromosome_map[chromosome]['map']) == chromosome_length:
                pass
            else:
                raise ValueError("Chromosome Maps are of unequal length")
                
        main_gene = random.choice(chromosome_list)
        for i in range(chromosome_length):
            if chromosome_map[main_gene]['map'][i] == 1:
                self.genome.append(gene)
            else:
                no_gene_found = True
                while no_gene_found:
                    gene = random.choice(chromosome_list)
                    if chromosome_map[gene]['map'][i] == 1:
                        self.genome.append(gene)
                        no_gene_found = False

    def evaluate(self):
        """
        Evaluates the solution to the optimization problem.
        """
        pass

    def raise_value_error(self):
        """
        Raises an error if there is no value for the given parameter.

        Parameters: 
            solution: class
                Instance of the solution class in the optimization.

        Written by Brian Andersen. 1/7/2019
        """
        casmo_keys = ['max_kinf','peak_kinf','eoc_kinf','peak_pin_power',
                      'enrichment','input_enrichment']
        param_list = []
        for param in self.parameters:
            if 'value' in self.parameters[param]:
                pass
            else:
                param_list.append(param)

        help_ = "No value has been written to the following parameters: "
        for par in param_list:
            help_ += " {}  ".format(par)
        help_ += ".\n"
        help_ += "The most likely error is an incorrect key.\n"
        
        help_ += "Supported parameter keys for casmo lattices are:\n"
        for par in casmo_keys:
            help_ += "{}\n".format(par)

        raise ValueError(help_)

def evaluate_function(solution):
    """
    Executes the method solution.evaluate(). This intermediate step is 
    necessary because Pool.map() is used to evaluate solutions in the 
    genetic algorithm. It needs to call a function in order to work, at
    least to the best of my knowledge. So this function is used.

    Parameters: 
        solution: class
            Instance of the problem type being optimized.

    Written by Brian Andersen. 1/7/2019
    """
    solution.evaluate()

    return solution

def test_evaluate_function(solution):
    """
    Executes the method solution.evaluate(). This intermediate step is 
    necessary because Pool.map() is used to evaluate solutions in the 
    genetic algorithm. It needs to call a function in order to work, at
    least to the best of my knowledge. So this function is used.

    Parameters: 
        solution: class
            Instance of the problem type being optimized.

    Written by Brian Andersen. 1/7/2019
    """
    #solution.evaluate()
    for param in solution.parameters:
        solution.parameters[param]['value'] = 10.*random.random()

    return solution

def convert_list_to_matrix(list_):
    """
    Converts lists to a NxN numpy matrix. Returns NxN numpy matrix.
    
    Developed specifically for number box testing but may have other uses as well.
    
    Parameters
        list_: list
        The genome that is to be turned into a numpy matrix.

    Written by Brian Andersen. 1/7/2019
    """
    chromosome_length = len(list_)
    number_rows=1
    number_columns = int(chromosome_length/number_rows)
    while number_rows < number_columns:
        number_rows +=1
        if chromosome_length % number_rows == 0:
            number_columns = int(chromosome_length/number_rows)

    matrix_ = numpy.zeros((number_rows,number_columns),dtype=int)

    position = 0
    for i in range(number_rows):
        for j in range(number_columns):
            matrix_[i][j] = list_[position]
            position += 1

    return matrix_

def return_triangular_string(list_):
    """
    Returns the genome as a single triangular string
    Used for octant symmetry in PWR fuel lattices and diagonal symmetry
    in BWR fuel lattices.

    Parameters
        list_: list
        The genome that is to be turned into a triangle.

    Written by Brian Andersen. 1/7/2019
    """
    triangular_string = ''
    column = 0
    row    = 1
    for item in list_:
        if type(item) is str:
            entry = item
        elif type(item) == list:
            print_line = "Function casmo_simulate_functions.return_triangular_string does not accept lists as an indice in the triangle"
            raise TypeError(print_line)
        elif type(item) == set:
            print_line = "Function casmo_simulate_functions.return_triangular_string does not accept sets as an indice in the triangle"
            raise TypeError(print_line)
        elif type(item) == dict:
            print_line = "Function casmo_simulate_functions.return_triangular_string does not accept dictionaries as an indice in the triangle"
            raise TypeError(print_line)
        else:             
            entry = str(item)
        triangular_string += entry.ljust(3)
        column +=1
        if row == column:
            triangular_string += "\n"
            column = 0
            row +=1

    return triangular_string

def compare_substrings(string1,string2):
    """
    Compares two strings. Returns true if string1 is a subset of string2.

    Parameters:
        string1: str
            The subset being searched for in string2.
        string2: str
            The string in which the subset is being searched for.

    Written by Brian Andersen. 1/7/2019 
    """
    string1elems = string1.strip().split()
    string2elems = string2.strip().split()
    desired_length = len(string1elems)
    match_count = 0
    for elem in string1elems:
        if elem in string2elems:
            match_count += 1
    if desired_length == match_count:
        return True
    else:
        return False

class Unique_Solution_Analyzer(object):
    """
    Ensures that all solutions are unique, and or haven't been previously analyzed.
    """
    def __init__(self,mutator):
        self.previous_solution_list = []
        self.mutator = mutator
        self.list_size_limit = 300  #arbitrarily set and can be changed. Implemented because in massive optimizations 
                                    #solutions at the beginning of the optimization are likely to be vastly different then
                                    #current solutions. So there is no need to keep them in the list.

    def analyze(self,new_solutions):
        """
        Analyzes all solutions in the provided solution list. 
        If any solutions are found to be in the previous solution list the solution is altered to 
        be one that has not been previously analyzed. 
        """
        for solution in new_solutions:
            while solution.genome in self.previous_solution_list:
                print("Uniqueness function Activated")
                for i in range(30):
                   solution.genome = self.mutator.reproduce(solution.genome)
            self.previous_solution_list.append(solution.genome)

        while len(self.previous_solution_list) > self.list_size_limit:
            self.previous_solution_list.pop()  #Simply removes the top solution list. The top solution should
                                               #be the oldest solution. 

        return new_solutions
