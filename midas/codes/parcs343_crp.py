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
    Interface used to run control rod sequence optimizations using PARCSv343 calculations.
    
    evaluate function creates working directory and prepares depletion file.
    For the time being, control rod sequence optimizations must be run using a template file.

    Written by Jake Mikouchi. 08/12/26
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
            if input.crp_partial['apply']:
                # partial sequence optimization
                with_template_partial_opt(solution, input, cwd, filename)
            else:
                # full sequence optimization
                with_template_full_opt(solution, input, cwd, filename)
        else: 
            raise ValueError("Axial Assembly optimization requires a template input file.")

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

def with_template_full_opt(solution, input, cwd, filename): 
    """
    Function to create parcs input file from template and populate input file with chromosome

    Written by Jake Mikouchi. 08/12/26
    """

## copy input file from template
    inp_template = str(cwd.joinpath(cwd / input.input_template['loc']))
    shutil.copy(inp_template, filename)

    with open(filename, "r") as file:
        lines = file.readlines()  

    with open(filename, "w") as ofile:
        realized_chromosome = prepare_chromosome_crp(input, solution.chromosome)
        num_zones = len(input.genome.keys())
        write_pattern_flag = False 
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

            elif "bank_def" in line.lower():
                write_pattern_flag = True 
                sequence_step = 1
                banks_step = 0
                num_banks = len(input.control_bank_conf)
                for patterns in range(len(input.bank_pattern)):
                    line = f"    bank_def   {sequence_step}   "
                    for rod in realized_chromosome[banks_step:banks_step+num_banks]:
                        line += f"{rod}   "
                    line += "\n"
                    ofile.write(line) 
                    banks_step += num_banks
                    sequence_step += 1

            else:
                ofile.write(line)   

def with_template_partial_opt(solution, input, cwd, filename): 
    """
    Function to create parcs input file from template and populate input file with chromosome
    this function allows to a partial sequence to be optimized rather than a full sequence 
    to reduce the strain on the optimizers.

    Written by Jake Mikouchi. 08/20/26
    """

## copy input file from template
    inp_template = str(cwd.joinpath(cwd / input.input_template['loc']))
    shutil.copy(inp_template, filename)

    with open(filename, "r") as file:
        lines = file.readlines()  

    with open(filename, "w") as ofile:
        realized_chromosome = prepare_chromosome_crp(input, solution.chromosome)
        num_zones = len(input.genome.keys())
        write_pattern_flag = False 
        banks_step = 0
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

            elif "bank_def" in line.lower():
                pattern_step = int(line.lower().split()[1])
                if pattern_step >= input.crp_partial['steps'][0] and pattern_step <= input.crp_partial['steps'][1]:
                    num_banks = len(input.control_bank_conf)
                    line = f"    bank_def   {pattern_step}   "
                    for rod in realized_chromosome[banks_step:banks_step+num_banks]:
                        line += f"{rod}   "
                    line += "\n"
                    ofile.write(line) 
                    banks_step += num_banks
                else: 
                    ofile.write(line)  
            else:
                ofile.write(line)   

def calc_weighted_notches(input, chromosome, parameters, output_index):
    """
    Uses the concave power law to compute a scalar objctive which prioritizes
    fewer rod insertions at deeper depths.
    Written by. Jake Mikouchi 10/14/2026
    """
    for parameter_key in input.objectives.keys():
        if parameters[parameter_key]['output_index'] == output_index:
            p = parameters[parameter_key]['factor']
            break
    realized_chromosome = prepare_chromosome_crp(input, solution.chromosome)
    notches = [input.cr_bounds["fully_removed"]-cr for cr in realized_chromosome]
    weighted_notches = 0 
    for i in notches:
        weighted_notches += i**p

def prepare_chromosome_crp(input_obj, chromosome):
    """
    Translates chromosome from normalized values representataions to actual control rod heights.
    This function also orders genes in the correct sequence that is required by parcs.
    
    Written by Jake Mikouchi. 08/12/2026
    """
    realized_chromosome = optools.Solution.crp_chromosome_realization(input_obj, chromosome)

    genome = input_obj.genome
    bank_patterns = input_obj.bank_pattern
    control_rod_banks = input_obj.control_bank_conf
    crb_bounds = input_obj.cr_bounds

    prepared_chromosome = []
    chromosome_index = 0 
    for sequence in bank_patterns:
        map = genome[sequence]['map'] 
        for rod in map: 
            if rod > 0: 
                prepared_chromosome.append(realized_chromosome[chromosome_index])
                chromosome_index += 1
            elif rod == 0: 
                prepared_chromosome.append(crb_bounds['fully_removed'])

    arranged_chromosome = []
    num_banks = len(control_rod_banks)
    sequence_step = 0 
    for sequence in bank_patterns:
        step_pattern = prepared_chromosome[sequence_step:sequence_step+num_banks]
        sequence_step += num_banks

        intermitent = [None for i in range(num_banks)] 
        for indx in control_rod_banks: 
            intermitent[indx-1] = step_pattern[control_rod_banks.index(indx)]

        arranged_chromosome.extend(intermitent)

    return arranged_chromosome
