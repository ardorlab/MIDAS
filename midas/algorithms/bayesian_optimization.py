import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from scipy.optimize import minimize
import logging

logger = logging.getLogger("MIDAS_logger")
class Bayesian_Optimization:
    """
    Class for performing Bayesian Optimization using a Gaussian Process surrogate
    model and binary encoding for categorical variables.

    Parameters:
        population: Class
            Class that contains the population size and stores the current solutions
        generation: Class
            Class that contains the current generation and their fitness values
    
    Written by Cole Howard. 1/1/2025
    """

    def __init__(self, input, noise=1e-6, random_state=None):
        self.input = input
        self.dimensions = self.parse_dimensions()  # List of length D, each an array of categories
        self.population = []       # Stores categorical points
        self.fitness_values = []   # Stores fitness values
        self.random_state = np.random.RandomState(random_state) #Initialize random state used throughout problem

        # Precompute the total number of binary features
        self.binary_dims = [self.calculate_binary_length(len(dim_cats)) for dim_cats in self.dimensions]
        self.total_features = sum(self.binary_dims)  # Total number of binary features

        # Scalers for input (X) and fitness values (y)
        self.scaler_x = StandardScaler()
        self.scaler_y = StandardScaler()

        #Noise for the GP model
        self.noise = noise

        #Define kernel function being used by GP as Matern
        self.kernel = Matern(
            length_scale=1.0,
            length_scale_bounds=(1e-4, 1e4),
            nu=self.input.kernel_smoothness
        )

        #Define GP being used
        self.gp = GaussianProcessRegressor(
            kernel=self.kernel,
            alpha=self.noise,
            n_restarts_optimizer=5,
            random_state=self.random_state
        )

        self.iterations = 0
        self.hyperparams = []
        self.hyperparam_convergence = False

    def parse_dimensions(self):
        """
        Parse genome data to create a dimension space for the optimization problem.
        Each element in dimension_space is a list of possible categories for that dimension.

        Written by Cole Howard. 10/29/2024
        """
        dimension_space = []
        map_shape = None
        location_options = {}

        for assembly, data in self.input.genome.items():
            assembly_map = data['map']
            if map_shape is None:
                map_shape = len(assembly_map)
            for idx, loc in enumerate(assembly_map):
                if loc == 1:
                    if idx not in location_options:
                        location_options[idx] = []
                    location_options[idx].append(assembly)

        for idx in range(map_shape):
            if idx in location_options:
                dimension_space.append(location_options[idx])
            else:
                dimension_space.append([None])

        return dimension_space

    @staticmethod
    def calculate_binary_length(num_categories):
        """
        Calculate the number of bits required to represent the number of categories in one dimension in binary.

        Written by Cole Howard. 1/1/2025
        """
        return int(np.ceil(np.log2(num_categories)))

    def categorical_to_binary_vector(self, categorical_individual):
        """
        Convert one categorical individual into a binary vector so the gaussian process can fit it

        Written by Cole Howard. 1/1/2025
        """
        binary_vector = []
        for i, cat in enumerate(categorical_individual):
            dim_cats = self.dimensions[i]
            if cat not in dim_cats:
                idx = 0  # Default to the first category if invalid
            else:
                idx = dim_cats.index(cat)

            num_bits = self.binary_dims[i]
            binary_rep = [int(x) for x in f"{idx:0{num_bits}b}"]
            binary_vector.extend(binary_rep)

        return np.array(binary_vector, dtype=float)

    def binary_vector_to_categorical(self, binary_vector):
        """
        Convert a 1D binary vector created by the GP into a categorical solution.

        Written by Cole Howard. 1/1/2025
        """
        categorical_sol = []
        start = 0
        for i in range(len(self.dimensions)):
            num_bits = self.binary_dims[i]
            end = start + num_bits
            binary_block = binary_vector[start:end]
            idx = int("".join(str(int(b)) for b in binary_block), 2)  # Convert binary to decimal

            # Constrain idx to the valid range of categories
            idx = min(idx, len(self.dimensions[i]) - 1)

            categorical_sol.append(self.dimensions[i][idx])
            start = end
        return categorical_sol

    def binary_encode_population(self, population):
        """
        Encode an entire population of categorical solutions into a 2D array
        of binary vectors.

        Written by Cole Howard. 1/1/2025
        """
        encoded_list = []
        #Recursively converting each individual within the population into a binary vector and appending that to the list
        for indiv in population:
            binary_vec = self.categorical_to_binary_vector(indiv)
            encoded_list.append(binary_vec)
        return np.array(encoded_list)

    def fit_gaussian_process(self):
        """
        Fit a Gaussian Process model to the current population and fitness values,
        using the binary-encoded version of each solution

        Written by Cole Howard. 1/1/2025
        """
        #Encode population in binary and scale in preparation to fit data
        encoded_population = self.binary_encode_population(self.population)
        scaled_population = self.scaler_x.fit_transform(encoded_population)

        #Scale fitness values and reshape the array as needed to fit
        scaled_fitness = self.scaler_y.fit_transform(
            np.array(self.fitness_values).reshape(-1, 1)
        ).ravel()
            
        #Fit GP model with most recent population/fitness data
        self.gp.fit(scaled_population, scaled_fitness)
        hyperparameters = self.gp.kernel_.get_params()
        lengthscale = hyperparameters['length_scale']
        self.hyperparams.append(lengthscale)
        if self.iterations>=10 and self.hyperparam_convergence == False:
            if abs((self.hyperparams[self.iterations-1]-self.hyperparams[self.iterations-9])/self.hyperparams[self.iterations-9]) <= 0.01:
                self.hyperparam_convergence = True
                logger.info("Hyperparameters have met convergence criteria; Hyperparameter fitting is now turned off")
                self.kernel = Matern(length_scale = self.hyperparams[self.iterations-1])
                self.gp = GaussianProcessRegressor(
                kernel=self.kernel,
                alpha=self.noise,
                optimizer=None,
                random_state=self.random_state
        )

    def tell(self, pop_list, fitness_list):
        """
        Update the surrogate model with the latest generation and its fitness values

        Written by Cole Howard. 1/1/2025
        """
        #Add the latest population and fitness values to their objects
        self.population.extend(pop_list)
        self.fitness_values.extend(fitness_list)
        #Fit GP model with the updated data set
        self.fit_gaussian_process()

    def expected_improvement(self, candidates, kappa):
        """
        Calculate EI for a set of candidate points in binary-encoded form.

        Written by Cole Howard. 1/1/2025
        """
        #Scale the candidates
        scaled_candidates = self.scaler_x.transform(candidates)
        #Predict mean and standard deviation for use in EI formula
        mean, std = self.gp.predict(scaled_candidates, return_std=True)

        #Shift fitness to avoid a negative 'best' value
        shifted_fitness = np.array(self.fitness_values) - np.min(self.fitness_values)
        best = np.min(shifted_fitness)

        #Replace standard deviations of zero with a very small number to avoid a division by zero error
        std = np.where(std == 0.0, 1e-9, std)
        #Calculate expected improvement
        z = (best - mean - kappa) / std
        ei = (best - mean - kappa) * norm.cdf(z) + std * norm.pdf(z)
        return ei

    def probability_of_improvement(self, candidates, kappa):
        """
        Calculate PI for a set of candidates in binary-encoded form.

        Written by Cole Howard. 1/1/2025
        """
        #Scale the candidates
        scaled_candidates = self.scaler_x.transform(candidates)
        #Predict mean and standard deviation for use in PI formula
        mean, std = self.gp.predict(scaled_candidates, return_std=True)

        #Shift fitness values to avoid negative 'best' value
        shifted_fitness = np.array(self.fitness_values) - np.min(self.fitness_values)
        best = np.min(shifted_fitness)

        #Replace zeros in standard deviation array with very small numbers to avoid division by zero error
        std = np.where(std == 0.0, 1e-9, std)
        #Calculate PI formula
        z = (best - mean - kappa) / std
        pi = norm.cdf(z)
        return pi

    def upper_confidence_bound(self, candidates, kappa):
        """
        Upper confidence bound acquisition function. Calculates the predicted minimum at upper bounds 

        Written by Cole Howard. 1/1/2025
        """
        #Scale the candidates to fit the data
        scaled_candidates = self.scaler_x.transform(candidates)
        #Predictions of means and standard deviations
        mean, std = self.gp.predict(scaled_candidates, return_std=True)
        #UCB formula
        return mean + kappa * std

    def lower_confidence_bound(self, candidates, kappa):
        """
        Lower confidence bound acquisition function. Calculates the predicted minimum at lower bounds

        Written by Cole Howard. 1/1/2025
        """
        #Scale the candidates to fit the data
        scaled_candidates = self.scaler_x.transform(candidates)
        #Get predicted means and standard deviations
        mean, std = self.gp.predict(scaled_candidates, return_std=True)
        #LCB formula
        return mean - kappa * std

    def real_vector_to_binary(self, x):
        """
        Convert a vector of real numbers between 0-1 into a binary vector

        Written by Cole Howard. 1/1/2025
        """
        binary_vector = []
        start = 0
        for num_bits in self.binary_dims:
            end = start + num_bits
            binary_block = x[start:end]  # Slice the real vector
            binary_block = (binary_block >= 0.5).astype(float)  # Threshold at 0.5
            binary_vector.extend(binary_block)
            start = end
        return np.array(binary_vector)

    def get_new_points(self, kappa, n_restarts):
        """
        Use L-BFGS-B in a continuous space to find a candidate
        that maximizes/minimizes an acquisition function,
        then convert it into a categorical vector

        Written by Cole Howard. 1/1/2025
        """
        #Determine which acquisition function to use based on the problem input
        def acquisition_function(x):
            binary_vec = self.real_vector_to_binary(x)
            if self.input.acquisition_function == 'EI':
                #Best point from EI has the highest value but we are using minimize, so we invert its value
                return -self.expected_improvement(binary_vec.reshape(1, -1), kappa)
            elif self.input.acquisition_function == 'PI':
                #Best point from PI has the highest value but we are using minimize, so we invert its value
                return -self.probability_of_improvement(binary_vec.reshape(1, -1), kappa)
            elif self.input.acquisition_function == 'UCB':
                return self.upper_confidence_bound(binary_vec.reshape(1, -1), kappa)
            elif self.input.acquisition_function == 'LCB':
                return self.lower_confidence_bound(binary_vec.reshape(1, -1), kappa)
        #Scaled problem bounds
        bounds = [(0, 1)] * self.total_features

        best_x = None
        #Initialize the best acq func value as positive inf so that lower values will replace it until the lowest value
        best_fun = float('inf')
        #Minimize the acquisition function with random restarts
        for _ in range(n_restarts):
            x0 = self.random_state.uniform(low=0.0, high=1.0, size=(self.total_features,))
            res = minimize(acquisition_function, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 1000000})

            if res.fun < best_fun:
                best_fun = res.fun
                best_x = res.x

        #Convert candidate point into a binary vector, then convert the binary vector into a categorical point we can use
        best_binary = self.real_vector_to_binary(best_x)
        candidate = self.binary_vector_to_categorical(best_binary)
        return candidate

    def reproduction(self, generation, gen):
        """
        Takes in the current population and their respective fitness values, and uses Bayesian Optimization to 
        suggest the next points.

        Written by Cole Howard. 1/1/2025

        Parameters:
            pop_list: object
                The current population.current list containing the list of individuals and their fitness values
        """
        self.iterations = gen
        pop_list = [soln.chromosome for soln in generation]
        #Invert incoming fitness because BO is a minimization tool and MIDAS prefers a maximum fitness
        fitness_list = [-1 * soln.fitness_value for soln in generation]
        if gen <= 100:
            self.tell(pop_list, fitness_list) #Fit the new data to the surrogate model

        candidates = []
        #Create the batch of candidate points for the next generation
        while len(candidates) < self.input.population_size:
            candidate = self.get_new_points(
                kappa=self.input.exploration_exploitation_factor,
                n_restarts=5
            )
            candidates.append(candidate)
        return candidates