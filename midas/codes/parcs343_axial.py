## Import Block ##
import os
import gc
import logging
import shutil
import numpy as np
from copy import deepcopy
from pathlib import Path
import subprocess
from subprocess import STDOUT
from midas.utils import optimizer_tools as optools
from midas_data import __parcs343exe__

from midas.codes.parcs343 import get_results

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

## Functions ##
def evaluate(solution, input):
    """
    Interface used to run axial assembly optimizations using PARCSv343 calculations.
    
    evaluate function creates working directory and prepares depletion file.
    For the time being, axial assembly optimizations must be run using a template file.

    Written by Jake Mikouchi. 07/20/26
    """
    
## Create and move to unique directory for PARCS execution
    cwd = Path(os.getcwd())
    indv_dir = cwd.joinpath(input.results_dir_name / Path(solution.name))
    if not indv_dir.exists():
        logger.debug(f"Creating new results directory: {indv_dir}")
        os.mkdir(indv_dir)
    logger.debug(f"Changing to new working directory: {indv_dir}")
    os.chdir(indv_dir)

## Prepare depletion file template
    with open('boc_exp.dep',"w") as depfile:
        depfile.write("\n BEGIN STEP\n\n EXP 3D MAP 1.0E+00\n\n")
        columncount = 0
        for i in range(1,input.num_assemblies+1):
            ## write column headers
            if columncount == 0:
                depfile.write(" k lb ")
            depfile.write(str(i).ljust(8))
            columncount += 1
            ## write rows for every 10 columns
            if columncount == 10:
                depfile.write('\n')
                for j in range(input.number_axial-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
                    depfile.write(' '+str(j).ljust(3))
                    for k in range(columncount):
                        depfile.write('{:.3f}'.format(input.boc_exposure).rjust(8))
                    depfile.write('\n')
                depfile.write('\n')
                columncount = 0
        ## write rows for leftover columns
        if columncount!= 0:
            depfile.write('\n')
            for j in range(input.number_axial-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
                depfile.write(' '+str(j).ljust(3))
                for k in range(columncount):
                    depfile.write('{:.3f}'.format(input.boc_exposure).rjust(8))
                depfile.write('\n')
            depfile.write('\n')
        depfile.write(' END STEP\n')
    
        filename = solution.name + '.inp'
        # create input file based on if an input template is given
        if input.input_template['apply']:
            with_template(solution, input, cwd, filename)
        else: 
            raise ValueError("Axial Assembly optimization requires a template input file.")

    ## Run PARCS INPUT DECK #!TODO: separate the input writing and execution into two different functions that are called in sequence.
        parcscmd = __parcs343exe__
        try:
            output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
        ## Get Results
            if 'Finished' in str(output): #job completed
                logger.debug(f"Job {solution.name} completed successfully in PARCSv343.")
                solution.parameters = get_results(solution.parameters, solution.name, solution.chromosome, input)
            
            else: #job failed
                if input.calculation_type in ['eq_cycle']:
                    try:
                        solution.parameters = eq_cycle_convergence(input, solution, filename, parcscmd, input.code_walltime) #iteratively try to find an intial guess that will converge
                    except Exception as e:
                        logger.error(f"Job {solution.name} has failed to converge with the following exception: {e}")
                        solution.parameters = get_results(solution.parameters, solution.name, solution.chromosome, input, job_failed=True)
                        
                else: #standard execution pathway
                    logger.warning(f"Job {solution.name} has failed!")
                    solution.parameters = get_results(solution.parameters, solution.name, solution.chromosome, input, job_failed=True)
        
        except subprocess.TimeoutExpired: #job timed out
            os.system('rm -f {}.parcs_pin*'.format(solution.name))
            logger.error(f"Job {solution.name} has timed out!")
            solution.parameters = get_results(solution.parameters, solution.name, solution.chromosome, input, job_failed=True)
        except subprocess.CalledProcessError as e: #PARCS returned an abort signal
            logger.error(f"Job {solution.name} has failed with the following exception: {e}")
            solution.parameters = get_results(solution.parameters, solution.name, solution.chromosome, input, job_failed=True)
        
        logger.debug(f"Returning to original working directory: {cwd}")
        os.chdir(cwd)
        gc.collect()

        return solution

def with_template(solution, input, cwd, filename): 
    """
    Function to create parcs input file from template and populate input file with chromosome

    Written by Jake Mikouchi. 07/20/26
    """

## copy input file from template
    inp_template = str(cwd.joinpath(cwd / input.input_template['loc']))
    shutil.copy(inp_template, filename)

    with open(filename, "r") as file:
        lines = file.readlines()  

    with open(filename, "w") as ofile:
        realized_chromosome = prepare_chromosome(input, solution.chromosome)
        num_zones = len(input.genome.keys())
        write_assemb_flag = False 
        for line in lines:
            ## change CaseID ##
            if "caseid" in line.lower():
                ofile.write('CASEID {}  \n'.format(solution.name))  
            ## apply th coupling
            elif 'int_th' in line.lower():
                if input.th_fdbk['apply']:
                    if input.th_fdbk['loc'] is None:
                        ofile.write("      TH_FDBK    T\n")
                        ofile.write("      INT_TH     T -1\n")
                    else: 
                        ofile.write(f"      INT_TH     T 1 '{input.th_fdbk['loc']}'\n")
                else:
                    ofile.write("      TH_FDBK    F\n")

            elif "assy_type    10" in line:
                write_assemb_flag = True 
                ofile.write(line) 
                assmb_type = 100
                temp_list = [{key:input.zone_discretize[key]['order']} for key in input.zone_discretize.keys()]
                zone_list = [0 for x in temp_list]
                for zone in temp_list:
                    zone_list[zone[list(zone.keys())[0]] - 1] = list(zone.keys())[0]
                for i in range(int(len(solution.chromosome)/num_zones)):
                    line = f"    assy_type    {assmb_type}"
                    index = 0 
                    for zone in zone_list:
                        if zone not in input.zone_options['constant_xs'].keys():
                            line += f"  {input.zone_discretize[zone]['num_nodes']}*{realized_chromosome[(i*num_zones+index)]}"
                            index += 1
                        else: 
                            line += f"  {input.zone_discretize[zone]['num_nodes']}*{input.zone_options['constant_xs'][zone]['designation']}"
                    line += "  FUEL \n"
                    ofile.write(line) 
                    assmb_type += 10
            else:
                ofile.write(line)   

def prepare_chromosome(input_obj, chromosome):
    """
    Translates chromosome from representing the names of the cross-sections to the numeric 
    designation which parcs uses to map cross-sections to assemblies.
    
    Written by Jake Mikouchi. 07/21/2026
    """
    realized_chromosome = []
    
    for gene in chromosome:
        realized_chromosome.append(input_obj.zone_options['zone_xs'][gene]['designation'])

    return realized_chromosome
