import h5py
# import pyarrow.parquet as pq
import os
import logging
from midas_data import __ipwr_lut__

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):
    """
    This function will find the output parameters for the NuScale SMR based on the SMR database of loading pattern calculations
    and write them into a dictionary in the solution.parameters object

    Written by Cole Howard. 10/29/2024
    updated by Jake Mikouchi. 1/3/25
    updated by Jake Mikouchi. 3/25/2026
    """
    if input.code_interface == 'ipwr_database_legacy':
        #Each objective is stored as one index in a single array within the hdf5 file, so I am getting each specific value
        objectives, BU= read_hdf5(solution.chromosome, input)
        #Create separate dictionary with parameters
        new_dict = {}
        new_dict["cycle_length"] = objectives[0]
        new_dict["fdeltah"] = objectives[1]
        new_dict["pinpowerpeaking"] = objectives[2]
        new_dict["max_boron"] = objectives[3]

        # Adding in burnup parameters
        new_dict["assembly_burnup"] = BU

    if input.code_interface == 'ipwr_database':
        #extract data from parquet file and return a dictionary 
        new_dict = read_parquet(solution.chromosome, input)

    #Only give the optimizer the parameters which were included in the input file
    for key in new_dict:
        if key in input.objectives:
            solution.parameters[key]["value"] = new_dict[key]

    return solution

def read_hdf5(soln, input): #TODO!: Comment code better
    """
    This function will search the hdf5 files containing the NuScale SMR database of loading patterns for the user-defined loading pattern,
    then return the parameters for that solution

    Parameters:
        individual: array
            array containing loading pattern for NuScale SMR (8 assemblies)

    Written by Cole Howard. 10/29/2024
    Updated by Jake Mikouchi. 12/16/25
    """
    individual = [input.fa_options['fuel'][sol]['type'] for sol in soln if sol in input.fa_options['fuel']]
    assembly_name = ''.join(map(str,individual))
    file_number = f'{individual[0]}{individual[1]}' #The number in the hdf5 file is the same as the first and second number of the loading pattern
    filepath = f"{__ipwr_lut__}/Solutions_{file_number}.hdf5"

    if not os.path.exists(filepath): #Check here to make sure the file being looked up actually exists
        raise FileNotFoundError(f'The file {filepath} does not exist. Make sure all assemblies in the LP are between 2-7.')
    #Open hdf5 file and index to the correct value, output the parameters
    with h5py.File(filepath,'r') as hdf5_file:

        assembly = hdf5_file[assembly_name]
        objectives = assembly["Objectives"][:] #Store all objectives for that LP
        BU = assembly["BU"][:] #Store burnup for each assembly in the LP
        hdf5_file.close()
        return objectives, BU
    
def read_parquet(soln, input): #TODO!: Comment code better
    """
    This function will extract the LP parameters from the parquet database.
    Unlike the hdf5, all LPs are saved in a single file making extraction easier.

    Written by Jake Mikouchi. 3/25/2026 
    """
    individual = [input.fa_options['fuel'][sol]['type'] for sol in soln if sol in input.fa_options['fuel']]
    core_name = ''.join(map(str,individual))
    filepath = f"{__ipwr_lut__}"

    if not os.path.exists(filepath): #ensures file exists
        raise FileNotFoundError(f'The file {filepath} does not exist.')

    data = pq.read_table(filepath, filters=[("Group", "=", core_name)]).to_pandas().set_index("Group")

    new_dict = {}
    new_dict["cycle_length"] = data.loc[core_name]["CL"]
    new_dict["fdeltah"] = data.loc[core_name]["FH"]
    new_dict["pinpowerpeaking"] = data.loc[core_name]["FQ"]
    new_dict["max_boron"] = data.loc[core_name]["BC"]
    new_dict["assembly_burnup"] = data.loc[core_name]["BU"]

    return new_dict