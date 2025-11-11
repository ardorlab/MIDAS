## Import Block ##
import math
import random
import logging
import numpy as np
from copy import deepcopy
import csv
import matplotlib.pyplot as plt
import os
from midas.utils.problem_preparation import LWR_Core_Shapes
"""
These are generic optimizer classes that are shared by all algorithms. #!TODO: can this be solved with the super().__init__ method?
"""


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

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
        
        self.current = 0 # dynamically tracks the iteration number of the optimizer
        self.initial = 0 # stores the starting generation number in the case of a restarted calculation

    def calculate_total_generations(self, number_changes):
        """
        Calculates the total number of generations in an optimization
        based on the number of changes.
        """
        return int(5*math.sqrt(number_changes))


class Solution():
    """
    This is the generic solution class to represent solutions in the optimization.
    "generate_initial" is the filtering function allowing for multiple generation approaches
    to be taken, based on the calculation type.

    Parameters: None

    Written by Brian Andersen. 1/7/2019
    Updated by Nicholas Rollins. 09/11/2024
    """
    def __init__(self,input_data, name=None):
        self.name = name
        self.parameters = {}
        self.chromosome = []
        self.fitness_value = float("NaN")
        self.input = input_data
    
    def generate_initial(self,calc_type,LWR_core_parameters,genome,batches=None):
        """
        Generates the initial solutions to the optimization problem by
          randomly generating a new chromosome.

        Parameters: 
            genome: Dictionary
                The genome portion of the dictionary settings file. 

        Written by Nicholas Rollins. 10/11/2024
        """
        if calc_type in ['single_cycle']:
            return self.LP_chromosome(genome, LWR_core_parameters)
        elif calc_type == 'eq_cycle':
            return self.EQ_chromosome(genome, batches, LWR_core_parameters)
        elif calc_type == 'lattice_physics':
            return self.lat_chromosome(genome,LWR_core_parameters)
        else:
            raise ValueError("Calculation Type not recognized; potential solution not generated.")
    
    def LP_chromosome(self,genome,LWR_core_parameters):
        genes_list = list(genome.keys())
        chromosome_length = []
        for gene in genes_list:
            chromosome_length.append(len(genome[gene]['map']))
        
        chromosome = []
        for i in range(max(chromosome_length)):
                gene_options = Gene_Validity_check.contraceptive_check(self.input, genes_list, genome, LWR_core_parameters,
                                                                                        [], chromosome, i)
                invalid = True
                antihang = 0
                while invalid:
                    antihang +=1
                    if antihang > 1000:
                        raise ValueError("Random solution generation failed after 1000 attempts. Check the variables maps and constraints.")
                    gene = random.choice(gene_options)
                    if genome[gene]['map'][i]: #check that the selected gene option is viable at this location.
                        chromosome.append(gene)
                        invalid = False
                    else:
                        gene_options.remove(gene)
        
        return chromosome
    
    def EQ_chromosome(self,genome,batches,LWR_core_parameters):
        """
        Encodes the shuffling scheme as a chromosome for Equilibrium cycle 
        calculations. Each entry is a tuple containing the batch name loaded in 
        that location, and either the name of the fuel assembly type (if fresh fuel)
        or the index in the chromosome of the shuffling source for the fuel assembly
        (if reloaded fuel).
        
        chromosome: list, e.g. [('batch_0','FAtype_1'),('batch_1',0)]
        
        Written by Nicholas Rollins. 10/14/2024
        """
        batches_list = list(batches.keys())
        chromosome_length = []
        for batch in batches_list:
            chromosome_length.append(len(batches[batch]['map']))
        
    ## Randomly generate zoning map.
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            zone_chromosome = [None]*max(chromosome_length) #zone map only, not encoded.
            chromosome_randindex = list(range(max(chromosome_length)))
            random.shuffle(chromosome_randindex)
            for i in chromosome_randindex:
                batch_options = Gene_Validity_check.contraceptive_check(self.input, batches_list, batches, LWR_core_parameters,
                                                                                        [], zone_chromosome, i)
                valid = False
                while not valid:
                    try:
                        batch = random.choice(batch_options)
                    except IndexError:
                        raise IndexError("Random chromosome generation has no valid solution. Is the input space over-constrained?")
                    if batches[batch]['map'][i]: #check that the selected gene option is viable at this location.
                        zone_chromosome[i] = batch
                        valid = True
                    else:
                        batch_options.remove(batch)
            chromosome_is_valid = Gene_Validity_check.abortive_check(self.input,batches_list, batches, LWR_core_parameters, zone_chromosome)
            attempts += 1
            if attempts > 10000:
                raise ValueError("Random chromosome generation has failed after 10,000 attempts. Is the input space over-constrained?")
        
    ## Assemble list of gene options for each batch.
        gene_options_dict = {}
        gene_options_dict[0] = list(genome.keys()) #list of fresh fuel options for batch 0
        for i in range(len(zone_chromosome)):
            batch_num = int(str(zone_chromosome[i]).replace(' ','_').split('_')[-1])
            # add loc as viable option for reloading for a subsequent batch.
            if batch_num+1 not in gene_options_dict:
                gene_options_dict[batch_num+1] = [i]
            else:
                gene_options_dict[batch_num+1].append(i)
        
    ## Randomly load fuel into shuffling scheme.
        chromosome = []
        for i in range(len(zone_chromosome)):
            chromosome.append((zone_chromosome[i],None))
        chromosome =Solution.EQ_reload_fuel(genome,LWR_core_parameters,chromosome)
        
        return chromosome
    
    def lat_chromosome(self,genome,core_parameters):
        
        """
        Generates an initial solution for the rod configuration for a lattice physics calculation.
        
        Written by Nicholas Rollins. 03/10/2025
        """
        symmetry = core_parameters[3]
        
        genes_list = list(genome.keys())
        chromosome_length = []
        for gene in genes_list:
            chromosome_length.append(len(genome[gene]['map']))
        
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            attempts += 1
            if attempts > 10000:
                raise ValueError("Random solution generation has failed after 10,000 attempts. Consider checking the variables maps and constraints.")

            chromosome = [None]*max(chromosome_length)
            #randomize the order by which locations are sampled to avoid bias
            randindexlist = [x for x in range(max(chromosome_length))]
            random.shuffle(randindexlist)
            for i in randindexlist:
                    gene_options = Gene_Validity_check.contraceptive_check(self.input, genes_list, genome, symmetry,
                                                                                        [], chromosome, i)
                    invalid = True
                    while invalid:
                        try:
                            gene = random.choice(gene_options)
                        except IndexError:
                            raise IndexError("Random solution generation failed after 1000 attempts. Check the variables maps and constraints.")
                        if genome[gene]['map'][i] == 1: #check that the selected gene option is viable at this location.
                            chromosome[i] = gene
                            invalid = False
                        else:
                            gene_options.remove(gene)
                            
            if Gene_Validity_check.abortive_check(self.input,genes_list,genome,core_parameters,chromosome):
                chromosome_is_valid = True
        
        return chromosome

    def SS_decoder(chromosome):
        """
        Extracts the encoded loading pattern from a chromosome 
        that represents a shuffling scheme.
        #!TODO: I believe this currently works under EQ cycle assumptions.
        
        Written by Nicholas Rollins. 10/14/2024
        """
        decoded_LP =  [None]*len(chromosome)
        for i in range(len(chromosome)):
            gene = chromosome[i][1]
            feed_fuel = False
            antihang = 0
            while not feed_fuel:
                if isinstance(gene, int): #FA is an index, implying reloaded fuel.
                    gene = chromosome[gene][1]
                    antihang += 1
                else:
                    decoded_LP[i] = gene
                    feed_fuel = True
                if antihang > 100:
                    raise ValueError("SS_decoder has discovered a circular dependency in a malformed shuffling scheme. This is highly irregular.")
        return decoded_LP

    def EQ_reload_fuel(genome, LWR_core_parameters, chromosome):
        """
        For an equilibrium cycle solution with batches selected for each location,
        missing FA types (i.e. None) are randomly filled from the available options.
        #!TODO: this could be expanded to general Shuffling Scheme reloading if the 
                gene_options_dict is decided differently.
        
        Written by Nicholas Rollins. 10/15/2024
        """
        ## Extract zones map
        zone_chromosome = [loc[0] for loc in chromosome]
        
        ## fetch the duplication multiplicity of each location when expanded to the full core.
        num_rows = LWR_core_parameters[0]
        num_cols = LWR_core_parameters[1]
        num_FA   = LWR_core_parameters[2]
        symmetry = LWR_core_parameters[3]
        multdict = LWR_Core_Shapes.get_symmetry_multiplicity(num_rows, num_cols, num_FA, symmetry)
        
        ## Assemble list of gene options for each batch.
        gene_options_dict = {}
        gene_options_dict[0] = list(genome.keys()) #list of fresh fuel options for batch 0
        for i in range(len(zone_chromosome)):
            batch_num = int(str(zone_chromosome[i]).replace(' ','_').split('_')[-1])
            # add loc as viable option for reloading for a subsequent batch.
            if batch_num+1 not in gene_options_dict:
                gene_options_dict[batch_num+1] = [i]
            else:
                gene_options_dict[batch_num+1].append(i)

        ## Remove already used gene options and resolve conflicts
        for i in range(len(zone_chromosome)):
            batch_num = int(str(zone_chromosome[i]).replace(' ','_').split('_')[-1])
            if chromosome[i][1] and batch_num != 0:
                try:
                    gene_options_dict[batch_num].remove(chromosome[i][1])
                except ValueError: #previous selection at this location made invalid by mutation.
                    chromosome[i] = (chromosome[i][0],None)
        
        ## Randomly load fuel into empty locations in shuffling scheme.
        for i in range(len(zone_chromosome)):
            batch_num = int(zone_chromosome[i].replace(' ','_').split('_')[-1])
            
            if not chromosome[i][1]: #location is missing a FA
            ## choose valid loading option before continuing.
                if not gene_options_dict[batch_num]:
                    raise ValueError(f"Failed to reload fuel; no source locations available for unassigned location of batch {batch_num}.\n{chromosome}") #!TODO: remove chromosome printout?
                valid = False
                attempt = 0
                while not valid:
                    attempt += 1
                    gene = random.choice(gene_options_dict[batch_num])
                    if batch_num == 0:
                        #!if genome[gene]['map'][i]: #check that the selected gene option is viable at this location. this requires decoding.
                        valid = True
                    #there must be enough symmetrical locs in the source to fill the symmetric locs in the target.
                    elif multdict[i] <= multdict[gene]:
                        valid = True
                    if attempt > 1000:
                        raise ValueError(f"Failed to reload fuel in shuffling scheme after 1,000 attempts for unassigned location of batch {batch_num}.\n{chromosome}") #!TODO: remove chromosome printout?
                    
                chromosome[i] = (chromosome[i][0],gene)
                if batch_num != 0:
                    gene_options_dict[batch_num].remove(gene)

        return chromosome


class Gene_Validity_check():  
    """
    Class used to constrain the parameters that make up a potential solution.
    The type of constraint is determined by the type of optimization being done. All solutions 
    should be passed to the gene_option_check function and that function will then distribute the 
    solution to the correct contraint checker. all solutions should be passed through the checks even
    if the constraint is trivial (i.e. None, or "unconstrained").
    
    Written by Jake Mikouchi. 09/04/2025
    """
    def contraceptive_check(input_obj, genes_list, genome, parameters, child, chromosome, indx):
        """
        Method for distributing to the correct contraceptive gene checking function
        Every type of contraceptive gene check should be accessed through this distributor.
        
        Written by Jake Mikouchi. 09/04/2025
        """
        if input_obj.calculation_type == 'single_cycle':
            valid_genes_list = Gene_Validity_check.calc_LWR_gene_options(genes_list, genome, parameters, child+chromosome[len(child):], indx)
        elif input_obj.calculation_type == 'eq_cycle':
            if None in chromosome:
                child_zone = chromosome
            else: 
                child_zone = [loc[0] for loc in child+chromosome[len(child):]]
            valid_genes_list = Gene_Validity_check.calc_LWR_gene_options(genes_list, genome, parameters, child_zone, indx)            
        elif input_obj.calculation_type == 'lattice_physics':
            valid_genes_list = Gene_Validity_check.calc_lat_gene_options(genes_list, genome, parameters, child+chromosome[len(child):], indx)
        else: 
            logger.warning('Unconstrained optimization')  
        
        return valid_genes_list

    def calc_LWR_gene_options(genes_list, genome, LWR_core_parameters, chromosome, index):
        """
        Constrain the available options for the chromosome based on
        the existing inventory.
        
        Written by Nicholas Rollins. 10/10/2024
        """
        ## fetch the duplication multiplicity of each location when expanded to the full core.
        num_rows = LWR_core_parameters[0]
        num_cols = LWR_core_parameters[1]
        num_FA   = LWR_core_parameters[2]
        symmetry = LWR_core_parameters[3]
        multdict = LWR_Core_Shapes.get_symmetry_multiplicity(num_rows, num_cols, num_FA, symmetry)
        gene_counts = LWR_Core_Shapes.count_in_LP(multdict,chromosome)

        ## if chromosome represents a shuffling scheme, not a loading pattern, the LP needs to be extracted before this step.
        valid_genes_list = []
        for gene in genes_list:
            if genome[gene]['constraint']:
                ctype = genome[gene]['constraint']['type']
                cvalue = genome[gene]['constraint']['value']
                if gene not in gene_counts:
                    gene_counts[gene] = 0
                if cvalue not in gene_counts:
                    gene_counts[cvalue] = 0
                if ctype == 'max_quantity':
                    #only include option if less than the max quantity have been already used.
                    if gene_counts[gene] < cvalue and (cvalue - gene_counts[gene]) > 1:
                        valid_genes_list.append(gene)
                elif ctype == 'less_than_variable':
                    #only include option if fewer than the target option have been already used.
                    if gene_counts[gene] < gene_counts[cvalue] and (gene_counts[cvalue] - gene_counts[gene]) > 1:
                        valid_genes_list.append(gene)
            else:
                valid_genes_list.append(gene)
        
        ## make sure that each gene option is valid for the gene location
        temp_gene_list = deepcopy(valid_genes_list)
        for gene in temp_gene_list:
            if not genome[gene]['map'][index] == 1:
                valid_genes_list.remove(gene)
        
        return valid_genes_list

    def calc_lat_gene_options(genes_list, genome, symmetry, chromosome, index):
        if isinstance(symmetry,list):
            gene_counts = LWR_Core_Shapes.count_in_lattice(symmetry[3],chromosome)  
        else: 
            gene_counts = LWR_Core_Shapes.count_in_lattice(symmetry,chromosome) 

        valid_genes_list = []
        for gene in genes_list:
            if genome[gene]['constraint']:
                ctype = genome[gene]['constraint']['type']
                cvalue = genome[gene]['constraint']['value']
                if gene not in gene_counts:
                    gene_counts[gene] = 0
                if cvalue not in gene_counts:
                    gene_counts[cvalue] = 0
                if ctype == 'max_quantity':
                    #only include option if less than the max quantity have been already used.
                    if gene_counts[gene] < cvalue and (cvalue - gene_counts[gene]) > 1:
                        valid_genes_list.append(gene)
                elif ctype == 'less_than_variable':
                    #only include option if fewer than the target option have been already used.
                    if gene_counts[gene] < gene_counts[cvalue] and (gene_counts[cvalue] - gene_counts[gene]) > 1:
                        valid_genes_list.append(gene)
            else:
                valid_genes_list.append(gene)
        
        ## make sure that each gene option is valid for the gene location
        temp_gene_list = deepcopy(valid_genes_list)
        for gene in temp_gene_list:
            if not genome[gene]['map'][index] == 1:
                valid_genes_list.remove(gene)
        
        return valid_genes_list

    def abortive_check(input_obj, genes_list, genome, parameters, child):
        """
        Method for distributing to the correct abortive gene checking function
        Every type of abortive gene check should be accessed through this distributor.
        
        Written by Jake Mikouchi. 09/08/2025
        """

        if input_obj.calculation_type == 'single_cycle':
            valid_chromosome = Gene_Validity_check.check_constraints(genes_list, genome, parameters, child)
        elif input_obj.calculation_type == 'eq_cycle':
            # input arguments are different based on if the solution is generated in the initial population or not
            is_tuple = False 
            for item in child:
                if isinstance(item, tuple):
                    is_tuple = True
            if is_tuple:
                child_zone = [loc[0] for loc in child]
                valid_chromosome = Gene_Validity_check.check_constraints(genes_list, genome, parameters, child_zone)
            else: 
                valid_chromosome = Gene_Validity_check.check_constraints(genes_list, genome, parameters, child)

        elif input_obj.calculation_type == 'lattice_physics':
            valid_chromosome = Gene_Validity_check.check_constraints(genes_list, genome, parameters, child)  
        else: 
            logger.warning('Unconstrained optimization')  
        
        return valid_chromosome

    def check_constraints(genes_list, genome, core_parameters, solution):
        """
        Check solution parameters against user-specified constraints on the input space.
        Returns True if the solution is valid and False if a constraint is violated.
        
        Written by Nicholas Rollins. 10/15/2024
        """
        if not genome: #! this implies that there are no constraints, but also no valid choices?
            return True
        
        ## fetch the duplication multiplicity of each location when expanded to the full core.
        num_rows = core_parameters[0]
        num_cols = core_parameters[1]
        num_FA   = core_parameters[2]
        symmetry = core_parameters[3]
        calc_type = core_parameters[4]
        if calc_type in ['single_cycle', 'eq_cycle']:
            multdict = LWR_Core_Shapes.get_symmetry_multiplicity(num_rows, num_cols, num_FA, symmetry)
        
        ## make sure that quantities of each gene type appearing in the solution are allowed.
        if calc_type != "lattice_physics":
            gene_counts = LWR_Core_Shapes.count_in_LP(multdict,solution)
        else:
            gene_counts = LWR_Core_Shapes.count_in_lattice(symmetry,solution)
        for gene in genes_list:
            if genome[gene]['constraint']:
                ctype = genome[gene]['constraint']['type']
                cvalue = genome[gene]['constraint']['value']
                if gene not in gene_counts:
                    gene_counts[gene] = 0
                if ctype == 'max_quantity': #quantity of gene type must be less than the max allowed quantity.
                    if gene_counts[gene] > cvalue:
                        return False
                elif ctype == 'less_than_variable': #quantity of gene type must be less than the quantity of the target variable.
                    if cvalue not in gene_counts:
                        gene_counts[cvalue] = 0
                    if gene_counts[gene] > gene_counts[cvalue]:
                        return False
        
        ## make sure that each gene option is valid for the gene
        for i in range(len(solution)):
            if not genome[solution[i]]['map'][i] == 1: #gene value not allowed at this location in chromosome.
                return False
            
        return True #if you haven't exited with "False" by this point, all constraints were passed.

  
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
            if param in ["assembly_burnup"]:
                pass
            else:
                pgoal = parameters[param]['goal']
                pweight = parameters[param]['weight']
                pvalue = parameters[param]['value']
                if not pvalue:
                    logger.error("No value was provided by MIDAS for the objective parameter '%s'. This is highly irregular.",param)
                
                if pgoal == 'maximize':
                    fitness += pvalue*pweight
                elif pgoal == 'minimize':
                    fitness -= pvalue*pweight
                elif pgoal == 'meet_target':
                    ptarget = parameters[param]['target']
                    fitness -= abs(ptarget - pvalue)*pweight
                elif pgoal == 'greater_than_target':
                    ptarget = parameters[param]['target']
                    penalty = (ptarget - pvalue)*pweight if ptarget - pvalue > 0.0 else 0.0
                    fitness -= penalty
                elif pgoal == 'less_than_target':
                    ptarget = parameters[param]['target']
                    penalty = (pvalue - ptarget)*pweight if pvalue - ptarget > 0.0 else 0.0
                    fitness -= penalty
        return fitness


class Solution_Reporting():
    """
    Class used to find and report best solution, optimization statistics, etc.

    Written by Cole Howard. 2/3/2025
    """
    def __init__(self):
        self.info = {}  # Dictionary to store all row information
        self.statistics = {} #Dictionary to store problem statistics
        self.generations = []

    def best_solution_information(self, filename):
        """
        Function used to parse output CSV file, store all data in a dictionary, and report best solution

        Written by Cole Howard. 2/3/2025
        """
        with open(filename, 'r') as file:
            reader = csv.reader(file)
            header = next(reader)  # Read the header row
            max_value = float('-inf')
            max_row = None

            for row in reader:
                try:
                    # Get the value in the fitness column
                    third_value = float(row[2])
                    
                    # Store the entire row info in self.info grouped by generation
                    generation = f"Gen_{row[0]}"  # Format generation key
                    individual = f"Ind_{row[1]}"  # Format individual key
                    row_data = {}
                    chromosome_start_index = header.index('Chromosome')
                    
                    for i, key_name in enumerate(header[:chromosome_start_index]):
                        if key_name not in ['Generation', 'Individual']:
                            try:
                                row_data[key_name] = float(row[i]) if '.' in row[i] else int(row[i])
                            except ValueError:
                                row_data[key_name] = row[i]
                    
                    # Store all chromosome values from the start index to the end of the row
                    row_data['Chromosome'] = row[chromosome_start_index:]
                    
                    if generation not in self.info:
                        self.info[generation] = {}
                    self.info[generation][individual] = row_data  # Store row under generation -> individual

                    # Check if it's the highest value so far
                    if third_value > max_value:
                        max_value = third_value
                        max_row = row
                except (ValueError, IndexError):
                    continue

        if max_row:
            result = {}
            chromosome_start_index = header.index('Chromosome')

            for i, key in enumerate(header[:chromosome_start_index]):
                try:
                    result[key] = float(max_row[i]) if '.' in max_row[i] else int(max_row[i])
                except ValueError:
                    result[key] = max_row[i]

            result['Chromosome'] = max_row[chromosome_start_index:]
            return result
        else:
            return None
    
    def optimization_statistics(self):
        """
        Function to calculate average fitness, max fitness, and standard deviation of fitness per generation
        and store within a dictionary.

        Written by Cole Howard. 2/3/2025
        """
        for gen_key, individuals in self.info.items():
            fitness_values = []
            
            for individual_data in individuals.values():
                if 'Fitness Value' in individual_data:
                    fitness_values.append(individual_data['Fitness Value'])
            
            if fitness_values:
                self.statistics[gen_key] = {
                    'Average_Fitness': np.mean(fitness_values),
                    'Max_Fitness': np.max(fitness_values),
                    'Std_Fitness': np.std(fitness_values)
                }
        #Get generations as an array to plot with
        self.generations = sorted([int(x.split('_')[1]) for x in self.statistics.keys()])
        return self.statistics
    
    def plot_optimization_statistics(self):
        #Store fitness statistics as arrays to plot with
        avg_fitness = [self.statistics[f'Gen_{gen}']["Average_Fitness"] for gen in self.generations]
        max_fitness = [self.statistics[f'Gen_{gen}']["Max_Fitness"] for gen in self.generations]
        std_fitness = [self.statistics[f'Gen_{gen}']["Std_Fitness"] for gen in self.generations]

        #Plot average fitness per generation
        plt.figure()
        plt.plot(self.generations, avg_fitness)
        plt.xlabel("Generation")
        plt.ylabel("Average Fitness Value")
        plt.title("Average Fitness per Generation")
        plt.grid()
        #Delete the plot if it already exists
        if os.path.exists("average_fitness_plot.png"):
            os.remove("average_fitness_plot.png")
        plt.savefig("average_fitness_plot.png")
        plt.close()

        #Plot max fitness per generation
        plt.figure()
        plt.plot(self.generations, max_fitness)
        plt.xlabel("Generation")
        plt.ylabel("Maximum Fitness Value")
        plt.title("Maximum Fitness per Generation")
        plt.grid()
        #Delete the plot if it already exists
        if os.path.exists("max_fitness_plot.png"):
            os.remove("max_fitness_plot.png")
        plt.savefig("max_fitness_plot.png")
        plt.close()

        #Plot standard deviation of fitness per generation
        plt.figure()
        plt.plot(self.generations, std_fitness)
        plt.xlabel("Generation")
        plt.ylabel("Fitness Value Standard Deviation")
        plt.title("Standard Deviation of the Fitness per Generation")
        plt.grid()
        #Delete the plot if it already exists
        if os.path.exists("std_fitness_plot.png"):
            os.remove("std_fitness_plot.png")
        plt.savefig("std_fitness_plot.png")
        plt.close()

    def plot_optimization_convergence(self):
        """Plot convergence by tracking the best fitness up to each generation."""
        best_fitness_values = []
        best_so_far = float('-inf')
        
        for gen in self.generations:
            if 'Max_Fitness' in self.statistics[f'Gen_{gen}']:
                best_so_far = max(best_so_far, self.statistics[f'Gen_{gen}']['Max_Fitness'])
                best_fitness_values.append(best_so_far)
        
        plt.figure()
        plt.plot(self.generations, best_fitness_values)
        plt.xlabel("Generation")
        plt.ylabel("Best Fitness")
        plt.title("Optimization Convergence")
        plt.grid()
        #Delete the plot if it already exists
        if os.path.exists("convergence_plot.png"):
            os.remove("convergence_plot.png")
        plt.savefig("convergence_plot.png")  # Save the plot
        plt.show()