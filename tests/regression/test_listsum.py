"""
Optimization problem: list summation 

Given a list size N and possible values for its elements in the range [0,9]
the objective is to find the list that maximizes the sum of its elements, which is [9,9,9,...,9]

There are 10^N possible solutions. 


From simple studies it seems we report the minimum, median, and maximum fitness values 
for different M number of random trials of a problem size of N=10:
   M = 1: 22/45/65
   M = 10: 47/60/81
   M = 100: 59/67/79
   M = 1000: 68/73/79
   M = 10000: 72/77/81
   M = 100000: 78/81/86
"""

import numpy as np
import logging

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):
    """
    This function will compute the sum of the list elements write them into a dictionary in the solution.parameters object

    Written by Gregory Delipei. 09/10/2025
    """
    #Create separate dictionary with parameters
    new_dict = {}
    new_dict["list_sum"] = list_summation(np.array(solution.chromosome, dtype=int))
   
    #Only give the optimizer the parameters which were included in the input file
    for key in new_dict:
        if key in input.objectives:
            solution.parameters[key]["value"] = new_dict[key]
        else:
            logger.info(f'Parameter {key} is available in NuScale database but not currently used')
    return solution

def list_summation(list_solution):
    return np.sum(list_solution)


def random_best_solution(M):
    current_best_solution = None
    current_best_fitness = -np.inf
    for _ in range(M):
        rnd_solution = np.random.randint(0, 10, N)
        rnd_fitness = list_summation(rnd_solution)
        if rnd_fitness > current_best_fitness:
            current_best_solution = rnd_solution
            current_best_fitness = rnd_fitness
    return current_best_solution, current_best_fitness

if __name__ == "__main__":
    N = 10
    best_solution = np.ones(N, dtype=int) * 9
    best_fitness = list_summation(best_solution)
    print(f"Best solution: {best_solution}")
    print(f"Best fitness: {best_fitness}")


    # Evaluate average performance for random search
    M = 100000
    Mout = 100
    rnd_solutions_list = []
    rnd_fitnesses_list = []
    for i in range(Mout):
        rnd_solution , rnd_fitness = random_best_solution(M)
        rnd_solutions_list.append(rnd_solution)
        rnd_fitnesses_list.append(rnd_fitness)

    # Convert lists to numpy arrays for easier indexing if not already
    rnd_fitnesses_array = np.array(rnd_fitnesses_list)
    rnd_solutions_array = np.array(rnd_solutions_list)

    # Get the index of the median fitness
    median_idx = np.argsort(rnd_fitnesses_array)[len(rnd_fitnesses_array) // 2]
    median_fitness = rnd_fitnesses_array[median_idx]
    median_solution = rnd_solutions_array[median_idx]

    # Get the index of the maximum fitness
    max_idx = np.argmax(rnd_fitnesses_array)
    max_fitness = rnd_fitnesses_array[max_idx]
    max_solution = rnd_solutions_array[max_idx]

    # Get the index of the minimum fitness
    min_idx = np.argmin(rnd_fitnesses_array)
    min_fitness = rnd_fitnesses_array[min_idx]
    min_solution = rnd_solutions_array[min_idx]

    
    print(f"Min random solution: {min_solution}")
    print(f"Min random fitness: {min_fitness}")


    print(f"Median random solution: {median_solution}")
    print(f"Median random fitness: {median_fitness}")

    print(f"Max random solution: {max_solution}")
    print(f"Max random fitness: {max_fitness}")

    



