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
        dimension_space = []

        # Assuming all maps have the same layout
        map_shape = None
        location_options = {}  # Tracks which assemblies can go in each location

        # Build location options based on assembly maps
        for assembly, data in self.input.genome.items():
            assembly_map = data['map']

            # Initialize map shape if not yet set
            if map_shape is None:
                map_shape = len(assembly_map)

            # Iterate over the map and store permissible assemblies for each position
            for idx, loc in enumerate(assembly_map):
                if loc == 1:
                    if idx not in location_options:
                        location_options[idx] = []
                    location_options[idx].append(assembly)

        # Create categorical dimensions for each map position
        for idx in range(map_shape):
            if idx in location_options:
                dimension_space.append(Categorical(location_options[idx], name=f"position_{idx}"))
            else:
                # Positions without valid assemblies can be left empty or filled as needed
                dimension_space.append(Categorical([None], name=f"position_{idx}"))

        return dimension_space
    
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