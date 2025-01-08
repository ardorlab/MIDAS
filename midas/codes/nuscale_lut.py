import h5py
import os
import numpy as np
import logging

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    """
    #Each objective is stored as one index in a single array within the hdf5 file, so I am getting each specific value
    objectives, BU, cost = read_hdf5(solution.chromosome, input)
    #Create separate dictionary with parameters
    new_dict = {}
    new_dict["cycle_length"] = objectives[0]
    new_dict["fdeltah"] = objectives[1]
    new_dict["pinpowerpeaking"] = objectives[2]
    new_dict["max_boron"] = objectives[3]
    new_dict["cycle_cost"] = cost

    # Adding in burnup parameters
    new_dict["max_burnup"] = np.max(BU)
    new_dict["min_burnup"] = np.min(BU)
    new_dict["average_burnup"] = np.mean(BU)

    #Only give the optimizer the parameters which were included in the input file
    for key in new_dict:
        if key in input.objectives:
            solution.parameters[key]["value"] = new_dict[key]
        else:
            logger.info(f'Parameter {key} is available in NuScale database but not currently used')

    return solution

def read_hdf5(soln, input): #TODO!: Comment code better
    """
    This function will search the hdf5 files containing the NuScale SMR database of loading patterns for the user-defined loading pattern,
    then return the parameters for that solution

    Parameters:
        individual: array
            array containing loading pattern for NuScale SMR (8 assemblies)

    Written by Cole Howard. 10/29/2024
    """
    individual = [input.fa_options['fuel'][sol]['type'] for sol in soln if sol in input.fa_options['fuel']]
    assembly_name = ''.join(map(str,individual))
    file_number = f'{individual[0]}{individual[1]}' #The number in the hdf5 file is the same as the first and second number of the loading pattern
    filepath = f"/cm/shared/databases/SMR_IPWR_DATABASE/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath): #Check here to make sure the file being looked up actually exists
        raise FileNotFoundError(f'The file Solutions_{file_number}.hdf5 does not exist. Make sure all assemblies in the LP are between 2-7.')
    #Open hdf5 file and index to the correct value, output the parameters
    with h5py.File(filepath,'r') as hdf5_file:

        assembly = hdf5_file[assembly_name]
        objectives = assembly["Objectives"][:] #Store all objectives for that LP
        #cost = assembly["cost"] #Store cost of LP #TODO!: Add cost back in once new dataase is generated
        BU = assembly["BU"][:] #Store burnup for each assembly in the LP
        cost = assembly["Cost"][()]
        hdf5_file.close()
        return objectives, BU, cost