import h5py
import os
import numpy as np
import logging

logger = logging.getLogger("MIDAS_logger")

def evaluate(soln, input): #!TODO: Put parameters in docstring, implement burnup and cycle cost into MIDAS
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    updated by Jake Mikouchi. 1/3/25
    """
    #Each objective is stored as one index in a single array within the hdf5 file, so I am getting each specific value
    objectives, BU = read_hdf5(soln.chromosome, input) #TODO!: Add cost back in

    results_dict = {}
    for res in ["cycle_length", "pinpowerpeaking", "fdeltah", "max_boron", "assembly_burnup"]:
        results_dict[res] = {}
        results_dict[res]['value'] = []

    results_dict["cycle_length"]["value"] = objectives[0]
    results_dict["fdeltah"]["value"] = objectives[1]
    results_dict["pinpowerpeaking"]["value"] = objectives[2]
    results_dict["max_boron"]["value"] = objectives[3]
    results_dict["assembly_burnup"]["value"] = BU 
    #solution.parameters['cycle_cost'] = cost

    for param in soln.parameters.keys():
        if param in results_dict:
            soln.parameters[param]['value'] = results_dict[param]["value"]
        else:
            logger.warning(f"Parameter '{param}' not supported in NuScale look up table results parsing.")

    return soln

def read_hdf5(soln, input): #TODO!: Comment code better
    """
    This function will search the hdf5 files containing the NuScale SMR database of loading patterns for the user-defined loading pattern,
    then return the parameters for that solution

    Parameters:
        individual: array
            array containing loading pattern for NuScale SMR

    Written by Cole Howard. 10/29/2024
    """
    individual = [input.fa_options['fuel'][sol]['type'] for sol in soln if sol in input.fa_options['fuel']]
    for ind in individual:
        if int(ind) not in [2, 3, 4, 5, 6, 7]:
            raise ValueError('"Type" under assembly options under "fuel" must be between 2-7. NuScale database only has these fuel types')
    assembly_name = ''.join(map(str,individual))
    file_number = f'{individual[0]}{individual[1]}' #The number in the hdf5 file is the same as the first and second number of the loading pattern
    filepath = f"/cm/shared/databases/SMR_IPWR_DATABASE/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath): #Check here to make sure the file being looked up actually exists
        raise FileNotFoundError(f'The file Solutions_{file_number}.hdf5 does not exist. Make sure all assemblies in the LP are between 2-7.')
    
    with h5py.File(filepath,'r') as hdf5_file:
        if assembly_name in hdf5_file:

            assembly = hdf5_file[assembly_name]

            objectives = assembly["Objectives"][:] #Store all objectives for that LP
            #cost = assembly["cost"] #Store cost of LP #TODO!: Add cost back in once new dataase is generated
            BU = assembly["BU"][:] #Store burnup for each assembly in the LP
            hdf5_file.close()
            return objectives, BU
        else:
            raise KeyError(f'The assembly "{assembly_name}" was not found in the file Solutions_{file_number}.hdf5. Please check that the assembly map was typed correctly.')
