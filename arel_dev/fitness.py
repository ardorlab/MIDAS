import os
import sys
import copy
import numpy

"""
This file is for storing all the possible ways to calculate the fitness or 
objective functions. In general, fitness or objective functions are designed
so that the smaller the fitness or objective value, the better the solution.

For methods that seek to increase the fitness of the optimized solution, a 
specialized class iswritten, which basically takes the inverse of the 
calculated solution fitness so that the optimization works as intended.

NOTE: All fitness classes should have the calculate method, which accepts 
itself and solution_list as the arguments. solution_list should be all 
solutions that need to have the fitness calculated. All solutions are 
calculated in this manner because some fitness methods, such as ranked fitness,
are based upon the corresponding solutions as well. 
So yeah
"""

class Fitness(object):
    """
    The generic fitness function. Requires user specified weights for every 
    design objective. Takes the Form F =   W1(minimize objectives) 
                                         - W2(maximize objectives)
                                         + W3(meet targets)
                                         + W4(satisfy objectives)
    WRitten by Brian Andersen. 1/18/2020
    """
    def __init__(self):
        pass

    def calculate(self,solution_list):
        """
        Calculates the generic fitness function, based on the listed fitness
        function above. Returns the solution list with evaluated fitnesses.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        for solution in solution_list:
            solution.fitness = 0              #Fitness is probably set to none right now.
            for param in solution.parameters:
                if solution.parameters[param]['goal'] == 'meet_target':
                    temp = copy.deepcopy(solution.parameters[param]['target'])
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                    temp = abs(temp)
                elif solution.parameters[param]['goal'] == 'less_than_target':
                    temp = copy.deepcopy(solution.parameters[param]['value'])
                    temp -= solution.parameters[param]['target']
                    if temp <= 0.0:
                        temp = 0.0
                elif solution.parameters[param]['goal'] == 'greater_than_target':
                    temp = copy.deepcopy(solution.parameters[param]['target'])
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                    if temp <= 0.:
                        temp = 0.0
                elif solution.parameters[param]['goal'] == 'minimize':
                    temp = copy.deepcopy(solution.parameters[param]['value'])
                elif solution.parameters[param]['goal'] == 'maximize':
                    temp = 0
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                temp *= copy.deepcopy(solution.parameters[param]['weight'])

                solution.fitness += temp

        return solution_list

class Genetic_Algorithm_Weighted(object):
    """
    Does the exact opposite of the regular fitness function. Since Simulated Annealing wants to minimize it
    and genetic algorithm is supposed to maximize it. 
    """
    def __init__(self):
        pass

    def calculate(self,solution_list):
        """
        Calculates the generic fitness function, based on the listed fitness
        function above. Returns the solution list with evaluated fitnesses.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 5/24/2020
        """
        for solution in solution_list:
            solution.fitness = 0              #Fitness is probably set to none right now.
            for param in solution.parameters:
                if solution.parameters[param]['goal'] == 'meet_target':
                    temp = copy.deepcopy(solution.parameters[param]['target'])
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                    temp = abs(temp)
                elif solution.parameters[param]['goal'] == 'less_than_target':
                    temp = copy.deepcopy(solution.parameters[param]['value'])
                    temp -= solution.parameters[param]['target']
                    if temp <= 0.0:
                        temp = 0.0
                elif solution.parameters[param]['goal'] == 'greater_than_target':
                    temp = copy.deepcopy(solution.parameters[param]['target'])
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                    if temp <= 0.:
                        temp = 0.0
                elif solution.parameters[param]['goal'] == 'minimize':
                    temp = copy.deepcopy(solution.parameters[param]['value'])
                elif solution.parameters[param]['goal'] == 'maximize':
                    temp = 0
                    temp -= copy.deepcopy(solution.parameters[param]['value'])
                temp *= copy.deepcopy(solution.parameters[param]['weight'])

                solution.fitness -= temp

        return solution_list


class Ranked_Fitness(Fitness):
    """
    Calculates Fitness based on rank of solutions. 

    Written by Brian Andersen. 1/18/2020
    """
    def __init__(self):
        Fitness.__init__(self)

    def calculate(self,solution_list):
        """
        Calculates the ranked fitness for the solution list.
        Returns the solution list with evaluated fitnesses.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        rank_dictionary = {}
        for solution in solution_list:
            for param in solution.parameters:
                if param in rank_dictionary:
                    pass
                else:
                    rank_dictionary[param] = {}
                    rank_dictionary[param]['goal'] = solution.parameters[param]['goal']
                    rank_dictionary[param]['list'] = []
                rank_dictionary[param]['list'].append(solution.parameters[param]['value'])

        for param in rank_dictionary:
            if rank_dictionary[param]['goal'] == 'minimize':
                rank_dictionary[param]['list'].sort()
            elif rank_dictionary[param]['goal'] == 'maximize':
                rank_dictionary[param]['list'].sort(reverse=True)
            elif rank_dictionary[param]['goal'] == 'meet_target':
                raise NotImplementedError
            elif rank_dictionary[param]['goal'] == 'less_than_target':
                rank_dictionary[param]['list'].sort()
            elif rank_dictionary[param]['goal'] == 'greater_than_target':
                rank_dictionary[param]['list'].sort(reverse=True)

        for solution in solution_list:
            temp_count = 1
            for param in solution.parameters:
                foo = solution.parameters[param]['value']
                for i,goo in enumerate(rank_dictionary[param]['list']):
                    if goo == foo:
                        if solution.parameters[param]['goal'] == 'less_than_target':
                            if solution.parameters[param]['target'] < foo:
                                temp_count += float(i)*solution.parameters[param]['weight']
                        elif solution.parameters[param]['goal'] == 'greater_than_target':
                            if solution.parameters[param]['target'] > foo:
                                temp_count += float(i)*solution.parameters[param]['weight']
                        else:
                            temp_count += float(i)*solution.parameters[param]['weight']
                        break

            solution.fitness = temp_count

        return solution_list

class Binned_Fitness(Fitness):
    """
    Calculates Fitness based on binned solution space values.

    Parameters
    ------------
    parameter_list : list
        List of the various parameters that are to be optimized.
    objective_list : list
        Whether each of the objectives should be minimized,maximized,
        or meet a target value. Follows the order of the parameter_list.
    bin_limits : list
        Each entry in the list should be a sublist of two values
        where the first entry in the sublist is the lower bin 
        limit. The second entry in the sublist is the upper bin
        limit. Follows the order of the parameter list. 

    Written by Brian Andersen. 1/18/2020
    """
    def __init__(self):
        Fitness.__init__(self)

    def calculate(self,solution_list):
        """
        Calculates the binned fitness for the solution list.
        Returns the solution list with evaluated fitnesses.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        for solution in solution_list:
            solution.fitness = 0              #Fitness is probably set to none right now.
            for param in solution.parameters:
                solution.fitness += solution.parameters[param]['bin']

        return solution_list

class ACDPF(Fitness):
    """
    Class for performing optimization using adaptively constrained
    discontinuous penalty function.

    Written by Brian Andersen. 1/18/2020
    """
    def __init__(self):
        Fitness.__init__(self)
        self.normalizing_count = 0 #Updated when the normilizing count is updated.
        self.limiting_values = {}
        self.normalizing_term = {}

    def calculate(self,solution_list):
        """
        Calculate solution fitness based on adaptive discontinuous penalty
        function.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        self.update_limiting_value(solution_list)
        self.update_normalizing_term(solution_list)

        track_file = open('optimization_track_file.txt','a')
        for param in solution_list[0].parameters:
            track_file.write("Normalizing term for parameter: {}  {} \n".format(param,self.normalizing_term[param]['term']))
            track_file.write("Limiting value for parameter: {} {} \n".format(param,self.limiting_values[param]))
        track_file.close()

        fitness_value = 0
        for solution in solution_list:
            fitness_value = 0
            for param in solution.parameters:

                goal = solution.parameters[param]['goal']
                constant = solution.parameters[param]['constant']
                value = solution.parameters[param]['value']
                norm = self.normalizing_term[param]['term']
                if(goal == "minimize"):
                    temp = self.limiting_values[param] - value
                elif(goal == "maximize"):
                    temp =  value - self.limiting_values[param]
                elif(goal == "meet_target"):
                    raise NotImplementedError
    
                temp = constant*(temp*temp)
                temp /= norm
                fitness_value += temp + 1.000
            
            solution.fitness = fitness_value
            
        return solution_list

    def update_normalizing_term(self,solution_list):
        """
        Calculates the normalizing term for the given solutions in
        the solution list for the given parameter.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """

        for solution in solution_list:
            for parameter in solution.parameters:
                goal = solution.parameters[parameter]['goal']
                value = solution.parameters[parameter]['value']
                if parameter in self.normalizing_term:
                    pass
                else:
                    self.normalizing_term[parameter] = {}
                    self.normalizing_term[parameter]["list"] = []
                    self.normalizing_term[parameter]['term'] = None
                if(goal == "minimize"):
                    temp = self.limiting_values[parameter] - value
                    temp *= self.limiting_values[parameter] - value
                elif(goal == "maximize"):
                    temp =  value - self.limiting_values[parameter]
                    temp *=  value - self.limiting_values[parameter]
                elif(goal == "meet_target"):
                    raise NotImplementedError
                self.normalizing_term[parameter]['list'].append(temp)
        
        for parameter in self.normalizing_term:
            self.normalizing_term[parameter]['term'] = numpy.mean(self.normalizing_term[parameter]['list'])

    def update_limiting_value(self,solution_list):
        """
        Determines the limiting value of the parameter for ACDPF 
        fitness function.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        min_max_dictionary = {}
        for param in solution_list[0].parameters:
            min_max_dictionary[param] = {}
            min_max_dictionary[param]['min'] = solution_list[0].parameters[param]['value']
            min_max_dictionary[param]['max'] = solution_list[0].parameters[param]['value']

        for solution in solution_list:
            for parameter in solution.parameters:
                value = solution.parameters[parameter]['value']
                if value < min_max_dictionary[parameter]['min']:
                    min_max_dictionary[parameter]['min'] = value
                if value  > min_max_dictionary[param]['max']:
                    min_max_dictionary[param]['max'] = value
                goal = solution.parameters[parameter]['goal']
                if parameter in self.limiting_values:
                    pass
                else:
                    self.limiting_values[parameter] = value
                if(goal == "minimize"):                          #Because this is for genetic algorithms the goal is not to minimize the    
                    if value < self.limiting_values[parameter]:  #end function over time, but to increase it as much as possible.
                        self.limiting_values[parameter] = value  #
                elif(goal == "maximize"):                        #
                    if value > self.limiting_values[parameter]:  # 
                        self.limiting_values[parameter] = value  #
                elif(goal == "meet_target"):
                    raise NotImplementedError

        track_file = open('optimization_track_file.txt','a')
        for param in min_max_dictionary:
            track_file.write("Minimum value for parameter: {} {} \n".format(param,min_max_dictionary[param]['min']))
            track_file.write("Maximum value for parameter: {} {} \n".format(param,min_max_dictionary[param]['max']))
        track_file.close()

class Fitness_Maxizer(Fitness):
    """
    Class used when the fitness of the optimization is to be maximized, rather
    than minimized.

    Parameters:
        method: class
            The standard fitness function that the optimization is using.

    Written by Brian Andersen. 1/18/2020
    """
    def __init__(self,method):
        Fitness.__init__(self)
        self.method = method

    def calculate(self,solution_list):
        """
        Calculates the fitness using the specified method, and then calculates
        the inverse to maximize the fitness of the population.

        Parameters:
            solution_list: List
                List of solutions whose fitness is to be calculated.

        Written by Brian Andersen. 1/18/2020
        """
        solution_list = self.method.calculate(solution_list)

        for solution in solution_list:
            solution.fitness = 10000./solution.fitness

        return solution_list
        
class Quantum_Fitness(Fitness):
    """
    Dr. Palmtag's fitness function doesn't fit with any of my other
    fitness functions so I am just adding it 
    """
    def __init__(self):
        Fitness.__init__(self)

    def calculate(self,solution_list):
        """
        Calculates the fitness for the quantum optimization project.
        """
        for solution in solution_list:
            exp = solution.parameters['exposure']['value']
            Fq = solution.parameters['PinPowerPeaking']['value']

            solution.fitness = 0.12*exp - Fq
            if Fq > 2.2:
                solution.fitness -= 2*(Fq-2.2)

        return solution_list