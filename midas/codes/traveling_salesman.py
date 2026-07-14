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
        1 - salesman_distance: total distance traveled
        2 - missed_desitnation: counts number of locations not visited
        
    Written by Jake Mikouchi. 07/01/2026
    """

    #dictionary with all parameters
    results_dict = {}

    # populate dictionary based on if the calcualtion was successful
    # there is likely no scenario where the calculation will not be successful for this code but still good practice
    if not job_failed:
        results_dict["distance"] = {'value': salesman_distance(chromosome, input_obj), 'output_index':1}
        results_dict["missed"] = {'value': missed_destination(chromosome, input_obj), 'output_index':2}
    else:
        results_dict["distance"] = {'value':10000, 'output_index':1}
        results_dict["missed"] = {'value': 10000, 'output_index':2}

    # populate parameters based on output index 
    for key in results_dict:
        output_index = results_dict[key]['output_index']
        for parameter_key in input_obj.objectives:
            if 'output_index' in input_obj.objectives[parameter_key].keys():
                if output_index == input_obj.objectives[parameter_key]['output_index']:
                    parameters[parameter_key]["value"] = results_dict[key]['value']

    return parameters


def salesman_distance(chromosome, input_obj):
    """
    calculates total distance traveled by the salesman

    Created by Jake Mikouchi. 07/06/2026
    """
    dimension = input_obj.c_length

    total_distance = 0
    for idx in range(dimension):
        if idx > 0:
            start_loc = input_obj.gene_options[chromosome[idx-1]]['coordinate']
            end_loc = input_obj.gene_options[chromosome[idx]]['coordinate']
            step_distance = np.sqrt( (end_loc[0] - start_loc[0])**2 + (end_loc[1] - start_loc[1])**2 )
            total_distance += step_distance
    return total_distance

def missed_destination(chromosome, input_obj):
    """
    checks to see if a destination (gene) was missed (if not all genes are used in the chromosome)

    Created by Jake Mikouchi. 07/06/2026
    """

    missed = 0
    for key in input_obj.genome.keys():
        if chromosome.count(key) < 1:
            missed += 1
    return missed
