from copy import deepcopy
from scipy.stats import spearmanr

class Termination_Criteria():
    """
    This class contains all termination criteria (TC) available in MIDAS.
    TC can be used with any optimation algorithm but was implemented specifically for use with GA.
    TC works by terminating an optimization run after certain cririteria are met.
    This can be used to save time and computational resources when an optimzation is failing to converge or has converged already.

    Written by Jake Mikouchi 1/2/25
    """
    
    def __init__(self):
        self.current_best_fitness = -10000
        self.previous_best_fitness = -10000
        self.consecutive_generations = 0
        self.current_burnup = [0, 0, 0, 0, 0, 0, 0, 0]
        self.previous_burnup = [0, 0, 0, 0, 0, 0, 0, 0]
        self.terminate = False
        pass

    def TC_methods(self, pop_list, TC):
        """
        Method for distributing to the requested Termination Criteria method.
        
        Written by Jake Mikouchi. 1/02/25
        """
        terminate = False

        if TC['method'] == 'consecutive':
            terminate = self.consecutive(pop_list, TC['termination_generations'])
        elif TC['method'] == 'spearman':
            terminate = self.spearman(pop_list, TC['termination_generations'])

        return terminate


    def consecutive(self, pop_list, termination_generations):
        """
        This Termination Criterion is meant to be very simple
        If the best fitness does not increase over a set number of generations, the optimization will automatically stop
        
        Jake Mikouchi 2/28/24
        """ 

        self.previous_best_fitness = deepcopy(self.current_best_fitness)

        for solution in pop_list:
            if solution.fitness_value > self.current_best_fitness:
                self.current_best_fitness = solution.fitness_value
        
        if self.current_best_fitness <= self.previous_best_fitness:
            self.consecutive_generations += 1
        else:
            self.consecutive_generations = 0
        
        if self.consecutive_generations < termination_generations:
            self.terminate = False
        else:
            self.terminate = True

        return self.terminate
    

    def spearman(self, pop_list, termination_generations):
        """
        This Termination Criterion utilizes the Spearman rank coefficient
        If the spearman rank does not change over a set number of generations, the optimization will automatically stop
        the spearman rank coefficient is a statistical method to measure the similarity between two datasets
        this uses the assembly burnup of the whole population to calculate the spearman rank and compares the burnup profile of a population to the previous population
        
        uses self.current_best_fitness to store the spearman rank coefficient in order to avoid unnecessary variables

        Jake Mikouchi 3/6/24
        # """

        self.previous_best_fitness = deepcopy(self.current_best_fitness)
        self.previous_burnup = deepcopy(self.current_burnup)
        self.current_burnup = [0, 0, 0, 0, 0, 0, 0, 0]

        # reads burnup and accounts for burunp for each assembly in the FULL core
        for solution in pop_list:
                LP_Burnup = solution.parameters["assembly_burnup"]["value"]
                for i in range(len(LP_Burnup)):
                    if i == 0:
                        self.current_burnup[i] += LP_Burnup[i]
                    elif i == 4 or i == 7:
                        self.current_burnup[i] += LP_Burnup[i]*8
                    else:
                        self.current_burnup[i] += LP_Burnup[i]*4
        # calculates the spearman rank coefficient
        self.current_best_fitness = spearmanr(self.previous_burnup,self.current_burnup).statistic
        
        if self.current_best_fitness == self.previous_best_fitness:
            self.consecutive_generations += 1
        else:
            self.consecutive_generations = 0
        
        if self.consecutive_generations < termination_generations:
            self.terminate = False
        else:
            self.terminate = True

        return self.terminate

    