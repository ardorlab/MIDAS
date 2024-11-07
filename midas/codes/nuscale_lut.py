import h5py
import os
import numpy as np
import logging

def evaluate(solution, input): #!TODO: Put parameters in docstring
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    """
    solutions = solution.chromosome
    #Each objective is stored as one index in a single array within the hdf5 file, so I am getting each specific value
    for soln in solutions:    
        objectives, BU = read_hdf5(soln) #TODO!: Add cost back in
        soln.parameters['cycle_length'] = objectives[0]
        soln.parameters['fdeltah'] = objectives[1]
        soln.parameters['pinpowerpeaking'] = objectives[2]
        soln.parameters['max_boron'] = objectives[3]
        #solution.parameters['cycle_cost'] = cost

        # Adding in burnup parameters in case it is used later on
        soln.parameters['max_burnup'] = max(BU)
        soln.parameters['min_burnup'] = min(BU)
        soln.parameters['average_burnup'] = np.mean(BU)

    return solutions

def read_hdf5(soln):
    """
    This function will search the hdf5 files containing the NuScale SMR database of loading patterns for the user-defined loading pattern,
    then return the parameters for that solution

    Parameters:
        individual: array
            array containing loading pattern for NuScale SMR

    Written by Cole Howard. 10/29/2024
    """
    individual = []
    for sol in soln:
        if sol == 'NSFA23':
            individual.append(2)
        elif sol == 'NSFA23GAD':
            individual.append(3)
        elif sol == 'NSFA39':
            individual.append(4)
        elif sol == 'NSFA39GAD':
            individual.append(5)
        elif sol == 'NSFA455':
            individual.append(6)
        elif sol == 'NSFA455GAD':
            individual.append(7)
    raise ValueError(soln)
    assembly_name = ''.join(map(str,individual))
    file_number = f'{individual[1]}' #The number in the hdf5 file is the same as the second number in the LP array for now, will change to first and second
    filepath = f"/cm/shared/databases/SMR_IPWR_DATABASE/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath): #Check here to make sure the file being looked up actually exists
        raise FileNotFoundError(f'The file Solutions_{file_number}.hdf5 does not exist. Make sure all assemblies in the LP are between 2-7.')
    
    with h5py.File(filepath,'r') as hdf5_file:
        if assembly_name in hdf5_file:

            assembly = hdf5_file[assembly_name]

            objectives = assembly["Objectives"][:] #Store all objectives for that LP
            #cost = assembly["cost"] #Store cost of LP #TODO!: Add cost back in once new dataase is generated
            BU = assembly["BU"][:] #Store burnup for each assembly in the LP
            return objectives, BU
        else:
            raise KeyError(f'The assembly "{assembly_name}" was not found in the file Solutions_{file_number}.hdf5. Please check that the assembly map was typed correctly.')
