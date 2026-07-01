import numpy as np
import logging
import os

from midas.utils import optimizer_tools as optools

logger = logging.getLogger("MIDAS_logger")


def evaluate(solution, input_obj):
    """
    Interface for the Levy Function.
    The levy function has an assymetric shape making optimization slightly more challenging than styblinksi-tang.
    evaluate function manages all operations requrired to interface with MIDAS

    Created by Jake Mikouchi. 07/01/2026
    """
    
    #Create results directory if it doesn't exist
    results_path = os.path.join(input_obj.results_dir_name, solution.name)
    os.makedirs(results_path, exist_ok=True)

    # get objectives for optimization
    solution.parameters = get_results(solution.parameters, solution.chromosome, input_obj, job_failed=False)

    return solution

def get_results(parameters, chromosome, input_obj, job_failed=False):
    """
    function used to retrieve all potential output parameters that may be used in an optimization

    Parameters: 
        1 - Levy Function value

    Written by Jake Mikouchi. 07/01/2026
    """

    realized_chromosome = optools.Solution.chromosome_realization(input_obj, chromosome)

    #dictionary with all parameters
    results_dict = {}

    # populate dictionary based on if the calcualtion was successful
    # there is likely no scenario where the calculation will not be successful for this code but still good practice
    if not job_failed:
        results_dict["stybtang"] = {'value':Stybtang(realized_chromosome, input_obj), 'output_index':1}
    else:
        results_dict["stybtang"] = {'value':1000000, 'output_index':1}

    # populate parameters based on output index 
    for key in results_dict:
        output_index = results_dict[key]['output_index']
        for parameter_key in input_obj.objectives:
            if 'output_index' in input_obj.objectives[parameter_key].keys():
                if output_index == input_obj.objectives[parameter_key]['output_index']:
                    parameters[parameter_key]["value"] = results_dict[key]['value']

    return parameters


def Stybtang(realized_chromosome, input_obj):
    """
    Directly calculates value of Levy function for dimension d
    where d is determined by the chromosome length

    Created by Jake Mikouchi. 07/01/2026
    """
    dimension = input_obj.c_length

    omega = []
    for i in range(dimension):
        omega.append(1 + (realized_chromosome[i]-1)/4 )

    sol_val = 0 

    sol_val += np.sin(np.pi * omega[0])**2
    sol_val += (omega[-1] - 1)**2 * (1 + np.sin(2 * np.pi * omega[-1])**2)

    for w in omega[1:-1]:
        sol_val += (w-1)**2 * (1 + 10 * np.sin(np.pi * w + 1)**2)

    return sol_val    