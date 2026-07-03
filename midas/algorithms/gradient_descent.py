import logging
from copy import deepcopy
import random
import numpy as np
from midas.utils import optimizer_tools as optools


class Gradient_Descent():
    """
    Class for performing optimization using the gradient descent module.
    Gradient Descent is a determinisitic optimization algorithm always converging to the same solution given the same inputs.
    Note that within MIDAS each iterations is refered to as a generation. This is not necessarily the correct nomenclature 
    for GD but it is named in this way for consistencey.

    Written by Jake Mikouchi. 07/01/2026
    """

    def __init__(self, input):
        self.input = input   
        self.learning_rate = input.learning_rate
        self.epsilon = input.epsilon
        self.active_sol = []
          
    def reproduction(self, pop_list, current_generation):
        """
        Computes the graient through numerical differentiation by perturbing each gene. 
        The computed gradient is then used to choose then next solution in the optimization process.
        
        Written by Jake Mikouchi. 07/01/2026
        """
        logger = logging.getLogger("MIDAS_logger")

        # select or update active solution
        if current_generation.current == 1: 
            self.active_sol = GD_reproduction.initial_selection(self, pop_list)    
        else:
            self.active_sol = GD_reproduction.gradient_update(self, pop_list)
        
        # create perturbations for fitness calculations
        chromosome_list = [self.active_sol]
        chromosome_list.extend(GD_reproduction.numerical_perturbations(self, self.active_sol))

        return chromosome_list
    
class GD_reproduction():
    """
    Functions for performing reproduction of chromosomes using GD
     
    Written by Jake Mikouchi. 07/01/26
    """

    def initial_selection(self, pop_list):
        """
        Parses through the initial population to select starting point for GD optimization.
        To avoid bias in the starting point, a roulette style selection is utilized.

        Written by Jake Mikouchi. 07/01/26
        """

        # if initial solution is provided then use that otherwise select with roulette
        if self.input.initial_population:
            winner = pop_list[0].chromosome
        else:
            probability_sum = 0
            selection_probability = {}
            selection_probability['low_bound'] = []
            selection_probability['up_bound']  = []
            
            shift_fitness_scale = min([soln.fitness_value for soln in pop_list]) #shift fitness values to start at zero (this corrects for negative fitness values)
            for solution in pop_list:
                selection_probability['low_bound'].append(probability_sum)
                probability_sum += (solution.fitness_value - shift_fitness_scale)
                selection_probability['up_bound'].append(probability_sum)

            value = random.random()
            value = value*probability_sum
            for j, solution in enumerate(pop_list):
                if selection_probability['low_bound'][j] <= value <= selection_probability['up_bound'][j]:
                    winner = solution.chromosome

        return winner


    def numerical_perturbations(self, chromosome):
        """
        generates 2 solutions by perturbing active solution.
        This is required so that the gradient can later be calculated.
        This is the only place where exection pathways for GD and SGD differ.

        Written by Jake Mikouchi. 07/01/26
        Updated by Jake Mikouchi 07/03/26
        """
        epsilon = self.epsilon
        core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry, self.input.calculation_type]
        all_gene_options = self.input.genome
        all_genes_list = list(self.input.genome.keys())

        pertrubs = []

        # standard gradient descent
        if not self.input.sgd:
            for i in range(len(chromosome)):
                perturbed_plus = deepcopy(chromosome)
                perturbed_minus = deepcopy(chromosome)

                gene_options = optools.Gene_Validity_check.contraceptive_check(self.input, all_genes_list, all_gene_options,
                                                                                        core_parameters, chromosome, [], i)
                # continous variable 
                if gene_options == [0,1]:
                    perturbed_plus[i] = perturbed_plus[i] + epsilon
                    perturbed_minus[i] = perturbed_minus[i] - epsilon
                    # extra check to ensure no solution goes outside of bounds
                    if perturbed_plus[i] > 1.0:
                        perturbed_plus[i] = 1.0 
                    elif perturbed_minus[i] < 0.0:
                        perturbed_minus[i] = 0.0 
                    pertrubs.append(perturbed_plus)
                    pertrubs.append(perturbed_minus)
                # discrete variable
                else: 
                    original_value = perturbed_plus[i]
                    lower_vals = gene_options[0:gene_options.index(original_value)]
                    higher_vals = gene_options[gene_options.index(original_value)+1:]
                    if not lower_vals:
                        perturbed_minus[i] = original_value
                    else:
                        perturbed_minus[i] = lower_vals[-1]
                    if not higher_vals:
                        perturbed_plus[i] = original_value
                    else:
                        perturbed_plus[i] = higher_vals[0]
                    pertrubs.append(perturbed_plus)
                    pertrubs.append(perturbed_minus)
        # stochastic gradient descent 
        else:
            gene_locations = [i for i in range(self.input.c_length)]
            random.shuffle(gene_locations)
            while len(pertrubs) < self.input.population_size-1:
                idx = random.choice(gene_locations)
                perturbed_plus = deepcopy(chromosome)
                perturbed_minus = deepcopy(chromosome)

                gene_options = optools.Gene_Validity_check.contraceptive_check(self.input, all_genes_list, all_gene_options,
                                                                                        core_parameters, chromosome, [], idx)
                # continous variable 
                if gene_options == [0,1]:
                    perturbed_plus[idx] = perturbed_plus[idx] + epsilon
                    perturbed_minus[idx] = perturbed_minus[idx] - epsilon
                    # extra check to ensure no solution goes outside of bounds
                    if perturbed_plus[idx] > 1.0:
                        perturbed_plus[idx] = 1.0 
                    elif perturbed_minus[idx] < 0.0:
                        perturbed_minus[idx] = 0.0 
                    pertrubs.append(perturbed_plus)
                    pertrubs.append(perturbed_minus)
                # discrete variable
                else: 
                    original_value = perturbed_plus[idx]
                    lower_vals = gene_options[0:gene_options.index(original_value)]
                    higher_vals = gene_options[gene_options.index(original_value)+1:]
                    if not lower_vals:
                        perturbed_minus[idx] = original_value
                    else:
                        perturbed_minus[idx] = lower_vals[-1]
                    if not higher_vals:
                        perturbed_plus[idx] = original_value
                    else:
                        perturbed_plus[idx] = higher_vals[0]
                    pertrubs.append(perturbed_plus)
                    pertrubs.append(perturbed_minus)
                
                gene_locations.pop(gene_locations.index(idx))

        return pertrubs

    def gradient_update(self, pop_list):
        """
        Uses fitness values from previously generated solutions to update the active soluition and explore the design space.

        Written by Jake Mikouchi 07/01/26
        """
        all_gene_options = self.input.genome
        all_genes_list = list(self.input.genome.keys())
        core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry, self.input.calculation_type]

        epsilon = self.epsilon
        lr = self.learning_rate

        active = []        
        # identify active solutions as list may have been reorganized by optimizer
        for soln in pop_list:
            if soln.chromosome == self.active_sol:
                active = soln.chromosome
                active_soln = soln

        grad = list(np.zeros_like(active))

        for idx in range(len(active)):
            perturbed_plus = None
            perturbed_minus = None
            for soln in pop_list:
                if soln.chromosome != self.active_sol:
                    if soln.chromosome[idx] > active[idx]:
                        perturbed_plus = soln
                    elif soln.chromosome[idx] < active[idx]:
                        perturbed_minus = soln
                
            if not perturbed_plus:
                perturbed_plus = active_soln
            elif not perturbed_minus:
                perturbed_minus = active_soln

            # calculate gradient for each gene
            grad[idx] = (perturbed_plus.fitness_value - perturbed_minus.fitness_value) / (2 * epsilon)

        # calculate new active solution using gradient
        new_active = list(np.zeros_like(active))
        for i in range(len(active)):
            gene_options = optools.Gene_Validity_check.contraceptive_check(self.input, all_genes_list, all_gene_options,
                                                                                    core_parameters, new_active, [], i)
            # continous variable
            if gene_options == [0,1]:
                new_active[i] = active[i] + lr * grad[i]
            # discrete variable
            else:
                original_value = active[i]
                lower_vals = gene_options[0:gene_options.index(original_value)]
                higher_vals = gene_options[gene_options.index(original_value)+1:]
                perturbation = lr * grad[i]
                temp_val = active[i] + perturbation
                if perturbation >= 0.0:
                    if not higher_vals:
                        new_active[i] = original_value
                    else:
                        if np.abs(higher_vals[0] - temp_val) < np.abs(temp_val - active[i]):
                            new_active[i] = higher_vals[0]
                        else:
                            new_active[i] = original_value
                elif perturbation < 0.0:
                    if not lower_vals:
                        new_active[i] = original_value
                    else:
                        if np.abs(temp_val - lower_vals[-1]) < np.abs(active[i] - temp_val):
                            new_active[i] = lower_vals[-1]
                        else:
                            new_active[i] = original_value                
        
        for i in range(len(active)):
            if new_active[i] > 1.0:
                new_active[i] = 1.0 
            elif new_active[i] < 0.0:
                new_active[i] = 0.0   

        return new_active