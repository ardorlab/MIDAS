from copy import deepcopy
from zipfile import Path

from midas.codes import surrogatemodel_BK
from midas.utils import optimizer_tools
from midas.utils import problem_preparation
from midas.utils.optimizer_tools import LP_map_loader

import os
import logging
import time as TT

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, mc_input, xsdict):
    """
    interface for evaluating a multicycle solution by running each cycle through the optimizer and collecting the results.

    C.King 2026-11
    """

    # initialize parameter result collector
    cycle_param_results = {}
    for obj_key in mc_input.objectives.keys():
        cycle_param_results[obj_key] = {}

    print("MULTICYCLE INPUT:")
    print("Chromosome:",solution.chromosome)
    print("Output directory:",mc_input.results_dir_name)
    # print("Num cycles:", mc_input.num_cycles)
    # print("Optimizer:",mc_input.methodology)
    # print("Pop Size:",mc_input.population_size)
    # print("Num Gens:",mc_input.num_generations)
    # print("Objectives:",mc_input.objectives)
    # print("Perturbs:",mc_input.perturbation_type)   
    print("Assembly Options:",mc_input.fa_options)
    print("Assembly Params/Genome:",mc_input.genome)

    print("Attempting to import the Optimizer!!")
    from midas.optimizer import Optimizer # import in function to prevent circular input
    for cycle_num in range(mc_input.num_cycles):
        # include a more interpretable printout of each 
        logger.info(f"Evaluating cycle {cycle_num+1} of {mc_input.num_cycles+1}")
        cycle_input = deepcopy(mc_input) # copy all settings from the input object to avoid modifying the original multicycle problem
        
        # Change settings to be for single-cycle. Not all settings need to be changed.
        # the current multicycle problem is not affected as this a copy

        # cycle_input.calculation_type = "single_cycle" # switch to single-cycle for evaluating the single-cycle problem
        # we don't have to switch to single_cycle since 

        # Some of the inputs for multicycle can be copied over to single-cycle. Others need to be adjusted
        # Since the inputs have been copied, we do not have to worry about breaking the original multicycle problem

        # GENERAL BLOCK
        # optimizer: simulated_annealing
        cycle_input.methodology = "simulated_annealing"
        # set initial_temperature and other PSA defaults
        # code_type: surrogate (maybe?)
        # results_dreictory_name: smoething not output_files
        cycle_input.results_dir_name = "output_files" + "_cycle_" + str(cycle_num + 1) # might need to adjust this for each cycle that is made
        # BIG ISSUE: might blow up file size
        
        # OPTIMIZER BLOCK
        cycle_input.population_size = 1 # pop size is 1 for simulated annealing
        cycle_input.num_generations = 100 # set a reasonable number of generations for simulated annealing
        cycle_input.objectives = {'max_boron': {'goal': 'less_than_target', 'target': 1300, 'weight': 1.0},
                                  'PinPowerPeaking': {'goal': 'less_than_target', 'target': 5000, 'weight': 2.1},
                                  'FDeltaH': {'goal': 'less_than_target', 'target': 1.48, 'weight': 6000},
                                  'cycle_length': {'goal': 'maximize', 'weight': 1.0}                                  
        }

        # ALGORITHM BLOCK
        cycle_input.perturbation_type = {"method": "perturb_by_swap", "num_permutations": 4}

        cycle_input.initial_temperature = 100
        if cycle_input.num_procs > 1:
            cycle_input.cooling_schedule = 'lam' 
        if cycle_input.num_procs <= 1:
            cycle_input.cooling_schedule = 'exponential_decrease'
            if cycle_input.cooling_schedule == "lam":
                raise ValueError(f"Cooling schedule '{cycle_input.cooling_schedule}' only available for parallel simulated annealing.")
        cycle_input.secondary_cooling_schedule = 'exponential_decrease'

        cycle_input.genome = solution.single_cycle_from_multicycle(cycle_num,cycle_input.fa_options) # single_cycle_from_mutlicycle modifies asssembly_options already
        # add burn assemblies will have to be in the process above!
        # parcs_data can remain the same
        
        print(f"CYCLE {cycle_num+1} INPUT:")
        print("Output directory:",cycle_input.results_dir_name)
        # print("Num cycles:", cycle_input.num_cycles) # this may be an issue for anything that uses the number of cycles as a check
        # print("Optimizer:",cycle_input.methodology)
        # print("Pop Size:",cycle_input.population_size)
        # print("Num Gens:",cycle_input.num_generations)
        # print("Objectives:",cycle_input.objectives)
        # print("Perturbs:",cycle_input.perturbation_type)
        print("Assembly Options:",cycle_input.fa_options)


        print("Assembly Params/Genome Before:",cycle_input.genome)
        print("Running LP map loader")
        LP_map_loader(cycle_input)

        print("Assembly Params/Genome After:",cycle_input.genome)

        


        # cycle_input = problem_preparation.Prepare_Problem_Values.prepare_calc(cycle_input) # prepare calculation values
        print(f"Optimizer for Cycle {cycle_num+1} is activated! Beware!")
        optimizer2 = Optimizer(cycle_input) # Run the optimizer with the cycle_input parameters
        optimizer2.build_optimizer()
        logger.info(f"Completed secondary Optimizer assembly for cycle {cycle_num+1}. Running secondary optimizer...\n")
        ts = TT.time()
        optimizer2.main() # reinitializing is probably not the most efficient way of doing this...
        # nevertheless, this can be resolved later
        logger.info(f'Execution completed after {TT.time()-ts:.5f} seconds')

        # extract info somehow

    print("MULTICYCLE INPUT MODIFICATION CHECK")
    print("Chromosome:",solution.chromosome)
    print("Output directory:",mc_input.results_dir_name)
    # print("Num cycles:",mc_input.num_cycles)
    # print("Optimizer:",mc_input.methodology)
    # print("Pop Size:",mc_input.population_size)
    # print("Num Gens:",mc_input.num_generations)
    # print("Objectives:",mc_input.objectives)
    # print("Perturbs:",mc_input.perturbation_type)
    print("Assembly Options:",mc_input.fa_options)
    print("Assembly Params/Genome:",mc_input.genome)

    # Sum up results of each cycle
    for item in cycle_param_results:
        cycle_param_results[item] = 5
        solution.parameters[item]['value'] = cycle_param_results[item]

    # process cycle_param_results to get final numbers for solution fitness function params

    return solution # return solution object and return to standard MIDAS execution
