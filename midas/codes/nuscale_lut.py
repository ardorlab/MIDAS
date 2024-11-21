import h5py
import os
import numpy as np
import logging

def evaluate(soln, input): #!TODO: Put parameters in docstring, implement burnup and cycle cost into MIDAS
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    """
    #Each objective is stored as one index in a single array within the hdf5 file, so I am getting each specific value
    objectives, BU = read_hdf5(soln.chromosome, input) #TODO!: Add cost back in
    soln.parameters["cycle_length"]["value"] = objectives[0]
    soln.parameters["fdeltah"]["value"] = objectives[1]
    soln.parameters["pinpowerpeaking"]["value"] = objectives[2]
    soln.parameters["max_boron"]["value"] = objectives[3]
    #solution.parameters['cycle_cost'] = cost

    # Adding in burnup parameters in case it is used later on
    #soln.parameters["max_burnup"]["value"] = max(BU) #TODO!: Add burnup parameters to solution object
    #soln.parameters["min_burnup"]["value"] = min(BU)
    #soln.parameters["average_burnup"]["value"] = np.mean(BU)

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
