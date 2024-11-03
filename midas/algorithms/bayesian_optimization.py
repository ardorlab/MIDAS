import importlib.metadata as importlib_metadata
from skopt import Optimizer as skOptimizer
from skopt.space import Real, Integer, Categorical

class Bayesian_Optimization():
   #!TODO: Finish docstrings (put in parameters)
    """
    Class for performing optimization using Bayesian Optimization, utilizing scikit-optimize python library

    Written by Cole Howard, 10/21/2024
    """
    def __init__(self, input):
        self.input = input
        self.optimizer = skOptimizer(
            dimensions = self.parse_dimensions(),
            n_initial_points = 0,
            acq_func = self.input.acquisition_function
        )

    def parse_dimensions(self):
        """
        Function to parse genome data and create dimension space for problem. This is needed because the Bayesian Optimization
        needs to know what possible values each variable can take. This function formats it so that scikit-optimize can understand it

        Written by Cole Howard, 10/21/2024
        """
        dimensions = []

        for param, attributes in self.input.genome.items():
            # Check if 'map' exists for the parameter (specifying physical locations)
            if 'map' in attributes:
                variable_map = attributes['map']
                valid_positions = [pos for pos, val in enumerate(variable_map) if val == 1]

                # For each valid position, treat the parameter as a categorical variable (param name)
                for pos in valid_positions:
                    # Create a categorical dimension using the parameter name
                    dim = Categorical([param])
                    dimensions.append(dim)
            else:
                # If no map exists, infer the type from the provided values
                variable_values = attributes.get('values', []) #!TODO: Using 'values' string as a placeholder until I figure out what to use here

                if all(isinstance(v, int) for v in variable_values):
                    dim = Integer(min(variable_values), max(variable_values))
                elif all(isinstance(v, float) for v in variable_values):
                    dim = Real(min(variable_values), max(variable_values))
                elif all(isinstance(v, str) for v in variable_values):
                    dim = Categorical(variable_values)
                else:
                    raise ValueError(f"Mixed or unknown types in variable '{param}' values: {variable_values}")

                dimensions.append(dim)

        return dimensions
    
    def reproduction(self, generation):
        """
        This function uses the ask and tell functions to update the prior distribution and generate another generation of points to sample from in conjunction
        with optimizer.py

        Written by Cole Howard, 10/21/2024
        """
        #Retrieve generation of individuals and respective fitness values from input self.population.current
        pop_list = [soln.chromosome for soln in generation]
        fitness_list = [soln.fitness_value for soln in generation]

        #Update prior distribution
        self.optimizer.tell(pop_list,fitness_list)

        #Create next set of points from which to sample
        next_generation=self.optimizer.ask(n_points=self.input.batch_size)
        return next_generation