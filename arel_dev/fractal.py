import os
import sys
import copy
import numpy
import random


"""
So, I think this class is less about fractalization now, cause I just can't figure out how to make it work. I really can't. I mean
the algorithm will remain unchanged and identical. It's just that I can't quite figure out how to justify it being a fractal. Because
if I could produce a picture that would be cool, but how do you make an N-dimensional picture? What you can easily make however is a 
lava flow. Then the number of different directions in which the lava can flow is based on the number of decision variables in a given point. 

Then, if the initial implementation of this algorithm doesn't work, a simulated annealing style cooling schedule could be implemented to ensure 
the algorithm isn't focusing on a local optimal solution. So if I get this idea off the ground I will re-make this file and call it something like lava or
whatever.
"""

class Fractalizer(object):
    """
    Generates new solutions for the fractal optimization method.
    """
    def __init__(self,decision_map):
        decision_length = None
        for key in decision_map:
            if not decision_length:
                decision_length = len(decision_map[key]['map'])
            elif len(decision_map[key]['map']) == decision_length:
                pass
            else:
                raise ValueError("Decision Maps are of unequal length")
            
        self.total_decision_variables = decision_length #The total number of decision variables.
        self.current_decision_variable = random.randint(0,self.total_decision_variables-1)      #Which decistion variable is to be edited in Fractilizer
        self._moving_forward = random.choice([True,False])
        
        self.variable_map = {}
        for i in range(self.total_decision_variables):
            self.variable_map[i] = []
            for key in decision_map:
                if decision_map[key]['map'][i]:
                    self.variable_map[i].append(key)

    def generate_new_solutions(self,solution_map):
        """
        Generates the new solutions based off the fractal optimization 
        methodology. Note that the actual number of new solutions returned
        is variable because it is based on the number of choices possible for that 
        decision variable.

        Parameters:
            solution_map: The layout of decision variables for the current best solution.
        """
        new_solution_list = []
        current_decision = solution_map[self.current_decision_variable]
        available_decisions = self.variable_map[self.current_decision_variable]
        for decision in available_decisions:
            if decision == current_decision:
                pass
            else:
                new_solution = copy.deepcopy(solution_map)
                new_solution[self.current_decision_variable] = decision
                new_solution_list.append(new_solution) 

        self._advance_decision_variable()

        return new_solution_list

    def _advance_decision_variable(self):
        """
        Sets what the next decision variable to alter will be. 
        """   
        if self._moving_forward:
            self.current_decision_variable += 1
            if self.current_decision_variable == self.total_decision_variables:
                self._moving_forward = False
        else:
            self.current_decision_variable -= 1
            if not self.current_decision_variable:
                self._moving_forward = True

class Fractal_Optimization(object):
    """
    Testing out this idea I had.
    """
    def __init__(self,solution,
                      number_solutions,
                      settings,
                      objective,
                      number_processors):

        self.solution = solution
        self.number_solutions  = number_solutions
        self.number_processors = number_processors
        self.file_settings      = settings
        self.objective_function = objective
        self.unique_solutions = 0
        self.converged_solutions = {}
    
    def main_in_parallel(self):
        raise NotImplementedError

    def main_in_serial(self):
        """
        Runs the optimization code in serial
        """
        current_best_list = {}
        fractal_list = {}
        for i in range(self.number_solutions):
            current_best = self.solution()
            current_best.name = f"Initial_{i}"
            current_best.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
            current_best.generate_initial()
            current_best.add_additional_information(self.file_settings)
            current_best.evaluate()
            current_best_list[i] = current_best

            fractal = Fractalizer(self.file_settings['genome']['chromosomes'])
            fractal_list[i] = fractal
        self.unique_solutions = i

        solution_unconverged = True
        while solution_unconverged:
            new_solution_dictionary = {}
            for i in range(self.number_solutions):
                new_solution_dictionary[i] = []
                new_decision_maps = fractal_list[i].generate_new_solutions(current_best_list[i].genome)
                for map_ in new_decision_maps:
                    new_solution = self.solution()
                    new_solution.name = f"solution_{self.unique_solutions}"
                    new_solution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                    new_solution.genome = map_
                    new_solution.add_additional_information(self.file_settings)
                    new_solution.evaluate()
                    self.unique_solutions += 1
                    new_solution_dictionary[i].append(new_solution)

            fractal_list = self._determine_position_in_space(fractal_list,new_solution_dictionary)
            solution_unconverged = self.determine_if_converged(fractal_list)

        

                #Need to write a new solution class instance block               
                #solution evaluator, and solution convergence function.

    def _determine_position_in_space(self,fractal_list,new_solutions):
        """
        Determines which solutions are going to become the current best solutions. Note, methodology has changed from
        fractals to lava. So it determines what the best way is for lava to flow down the mountain, with the solution 
        fitness representing position on the mountain. 
        """
        evaluated_solutions = {}
        for i in range(self.number_solutions):
            evaluated_solutions = self.objective_function(new_solutions[i])
            for solution in evaluated_solutions:
                if solution.fitness < fractal_list[i].fitness:
                    fractal_list[i] = copy.deepcopy(solution)
                

        return fractal_list

    def _determine_if_converged(self,fractal_list,fractallers):
        """
        Determines if the optimization should continue to progress or not.
        """
        if not self.converged_solutions:
            self.converged_solutions = copy.deepcopy(fractal_list)
        
        differences_found = False
        for one,two in zip(self.converged_solutions,fractal_list):
            if one.genome == two.genome:
                pass
            else:
                differences_found = True
                break

        














