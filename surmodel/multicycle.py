from copy import deepcopy
from zipfile import Path

from skopt import Optimizer
from midas.codes import surrogatemodel_BK
from midas.utils import optimizer_tools
from midas.utils import problem_preparation
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

    cycle_param_results = {}

    print("MULTICYCLE INPUT:",mc_input)

    for cycle_num in range(mc_input.num_cycles):
        # include a more interpretable printout of each 
        logger.info(f"Evaluating cycle {cycle_num+1} of {mc_input.num_cycles+1}")
        cycle_input = deepcopy(mc_input) # copy all settings from the input object to avoid modifying the original multicycle problem
        
        # Change settings to be for single-cycle. Not all settings need to be changed.
        # the current multicycle problem is not affected as this a copy

        # cycle_input.calculation_type = "single_cycle" # switch to single-cycle for evaluating the single-cycle problem
        # we don't have to switch to single_cycle since 

        cycle_input = problem_preparation.prep_inp.prepare_calc(cycle_input) # prepare calculation values
        # Some of the inputs for multicycle can be copied over to single-cycle. Others need to be adjusted
        # Since the inputs have been copied, we do not have to worry about breaking the original multicycle problem

        # GENERAL BLOCK
        # optimizer: simulated_annealing
        cycle_input.optimizer = "simulated_annealing"
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

        cycle_input.assembly_parameters = solution.single_cycle_from_multicycle(solution,cycle_num,solution.chromosome) # single_cycle_from_mutlicycle modifies asssembly_options already
        # add burn assemblies will have to be in the process above!
        # parcs_data can remain the same
        
        print(f"CYCLE {cycle_num+1} INPUT:",cycle_input)
        # optimizer = Optimizer(cycle_input) # Run the optimizer with the cycle_input parameters
        # optimizer.build_optimizer()
        # logger.info(f"Completed secondary Optimizer assembly for cycle {cycle_num+1}. Running secondary optimizer...\n")
        # ts = TT.time()
        # optimizer.main() # reinitializing is probably not the most efficient way of doing this...
        # nevertheless, this can be resolved later
        # logger.info(f'Execution completed after {TT.time()-ts:.5f} seconds')

        # extract info somehow

    print("MULTICYCLE INPUT MODIFICATION CHECK:",mc_input)
    # process cycle_param_results to get final numbers for solution fitness function params
    cycle_param_results["cycle_length"]['value'] = 5 # TODO: this

    return cycle_param_results # solution parameters