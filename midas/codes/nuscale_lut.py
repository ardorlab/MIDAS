import h5py
import os

def evaluate(solution, input): #!TODO: Put parameters in docstring
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations

    Written by Cole Howard. 10/29/2024
    """
    for soln in solution:
        pass

    return solution
# parameters 
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
    file_number = individual[1]
    filepath = f"/cm/shared/databases/SMR_IPWR_DATABASE/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath):
        raise FileNotFoundError(f'The file Solutions_{file_number}.hdf5 does not exist. Make sure the second number in the loading pattern is between 2-7.')
    
    with h5py.File(filepath,'r') as hdf5_file:
        if assembly_name in hdf5_file:

            assembly = hdf5_file[assembly_name]

            objectives = assembly["Objectives"][:]
            BU = assembly["BU"][:]
            return objectives, BU
        else:
            raise KeyError(f'The assembly "{assembly_name}" was not found in the file Solutions_{file_number}.hdf5. Please check that the assembly map was typed correctly.')
