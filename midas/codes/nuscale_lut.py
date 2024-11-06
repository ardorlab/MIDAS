import h5py
import os
import numpy as np
import logging

def evaluate(solutions, input): #!TODO: Put parameters in docstring
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    """
    for soln in solutions:
        #Each objective is stored as one index in a single array withon the hdf5 file, so I am getting each specific value
        objectives, cost, BU = read_hdf5(solutions)
        soln.parameters['cycle_length'] = objectives[0]
        soln.parameters['fdeltah'] = objectives[1]
        soln.parameters['pinpowerpeaking'] = objectives[2]
        soln.parameters['max_boron'] = objectives[3]
        soln.parameters['cycle_cost'] = cost

        # Adding in burnup parameters in case it is used later on
        soln.parameters['max_burnup'] = max(BU)
        soln.parameters['min_burnup'] = min(BU)
        soln.parameters['average_burnup'] = np.mean(BU)

    return solutions

def read_hdf5(individual):
    """
    This function will search the hdf5 files containing the NuScale SMR database of loading patterns for the user-defined loading pattern,
    then return the parameters for that solution

    Parameters:
        individual: array
            array containing loading pattern for NuScale SMR

    Written by Cole Howard. 10/29/2024
    """
    assembly_name = ''.join(map(str,individual))
    file_number = f'{individual[0]}{individual[1]}' #The number in the hdf5 file is the same as the second number in the LP array for now, may change later
    filepath = f"/cm/shared/databases/SMR_IPWR_DATABASE/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath): #Check here to make sure the file being looked up actually exists
        raise FileNotFoundError(f'The file Solutions_{file_number}.hdf5 does not exist. Make sure all assemblies in the LP are between 2-7.')
    
    with h5py.File(filepath,'r') as hdf5_file:
        if assembly_name in hdf5_file:

            assembly = hdf5_file[assembly_name]

            objectives = assembly["Objectives"][:] #Store all objectives for that LP
            cost = assembly["Cost"] #Store cost of LP
            BU = assembly["BU"][:] #Store burnup for each assembly in the LP
            return objectives, cost, BU
        else:
            raise KeyError(f'The assembly "{assembly_name}" was not found in the file Solutions_{file_number}.hdf5. Please check that the assembly map was typed correctly.')
