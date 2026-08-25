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

from midas.codes.parcs343 import calc_cycle_length
from midas.codes.parcs343 import calc_cpr
from midas.codes.parcs343 import calc_lhgr
from midas.codes.parcs343 import calc_aplhgr
from midas.codes.parcs343 import eq_cycle_convergence
from midas.codes.parcs343 import next_binary_search

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
                solution.parameters = get_results(solution.parameters, solution.name, input)
            
            else: #job failed
                if input.calculation_type in ['eq_cycle']:
                    try:
                        solution.parameters = eq_cycle_convergence(input, solution, filename, parcscmd, input.code_walltime) #iteratively try to find an intial guess that will converge
                    except Exception as e:
                        logger.error(f"Job {solution.name} has failed to converge with the following exception: {e}")
                        solution.parameters = get_results(solution.parameters, solution.name, input, job_failed=True)
                        
                else: #standard execution pathway
                    logger.warning(f"Job {solution.name} has failed!")
                    solution.parameters = get_results(solution.parameters, solution.name, input, job_failed=True)
        
        except subprocess.TimeoutExpired: #job timed out
            os.system('rm -f {}.parcs_pin*'.format(solution.name))
            logger.error(f"Job {solution.name} has timed out!")
            solution.parameters = get_results(solution.parameters, solution.name, input, job_failed=True)
        except subprocess.CalledProcessError as e: #PARCS returned an abort signal
            logger.error(f"Job {solution.name} has failed with the following exception: {e}")
            solution.parameters = get_results(solution.parameters, solution.name, input, job_failed=True)
        
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
        realized_chromosome = prepare_chromosome(input, solution.chromosome)
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
        realized_chromosome = prepare_chromosome(input, solution.chromosome)
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

def get_results(parameters, filename, input, job_failed=False):
    """
    function used to retrieve all potential output parameters that may be used in an optimization

    Parameters: 
        1 - cycle length: 
        2 - pxy: assembly level planar peaking
        3 - pxyz: assembly level  volumetric peaking
        4 - pin power peaking: Fq peaking factor
        5 - f delta h: enthalpy rise peaking factor
        6 - max_boron: 
        7 - chfr: critical heat flux ratio
        8 - keff_min: 
        9 - keff_max: 
        10 - keff_diff: difference between minimum and maximum keff
        11 - cpr: critical power ratio
        12 - lhgr: linear heat generation rate
        13 - aplhgr: average planar lhgr
        14 - axoff_boc: axial offset at begining of cycle
        15 - axoff_eoc: axial offset at end of cycle
        16 - axoff_moc: axial offset at middle of cycle

    Updated by Jake Mikouchi. 07/22/2026
    """
    ## Prepare container for results
    results_dict = {}
        
    if not job_failed:
        ## Read file for parsing
        with open(filename + ".parcs_dpl", "r") as ofile:
            filestr = ofile.read()
        
        ## Split file by section
        res_str = filestr.split('===============================================================================')
        res_str = res_str[-1].split('_______________________________________________________________________________')
        res_str = res_str[0].split('\n')
        
        ## Parse raw values by timestep
        efpd_list = []; boron_list = []; keff_list = []; pxy_list = []; pxyz_list = []; fq_list = []; fdh_list = []; chfr = []; axoff = []
        for i in range(2, len(res_str)-1):
            res_val=res_str[i].split()
            
            efpd_list.append(float(res_val[9]))
            boron_list.append(float(res_val[14]))
            keff_list.append(float(res_val[2]))            
            pxyz_list.append(float(res_val[7]))
            pxy_list.append(float(res_val[6]))
            fdh_list.append(float(res_val[21]))
            fq_list.append(float(res_val[22]))
            chfr.append(float(res_val[23]))
            axoff.append(float(res_val[4]))
        
        del filestr, res_str, res_val #unload file contents to clean up memory

        if input.crp_partial['apply']: 
            efpd_list = efpd_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            boron_list = boron_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            keff_list = keff_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            pxyz_list = pxyz_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            pxy_list = pxy_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            fdh_list = fdh_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            fq_list = fq_list[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            chfr = chfr[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]
            axoff = axoff[input.crp_partial['steps'][0]-1:input.crp_partial['steps'][1]]

        results_dict["cycle_length"] = {'value': calc_cycle_length(efpd_list,boron_list,keff_list,input.eoc_extrapolate), 'output_index':1}
        results_dict["pxy"] = {'value':max(pxy_list), 'output_index':2}
        results_dict["pxyz"]= {'value':max(pxyz_list), 'output_index':3}
        results_dict["pinpowerpeaking"] = {'value':max(fq_list), 'output_index':4}
        results_dict["fdeltah"] = {'value':max(fdh_list), 'output_index':5}
        results_dict["max_boron"] = {'value':max(boron_list), 'output_index':6}
        results_dict["chfr"] = {'value':min(chfr), 'output_index':7}
        
        results_dict["keff_min"] = {'value':min(keff_list), 'output_index':8}
        results_dict["keff_max"] = {'value':max(keff_list), 'output_index':9}
        results_dict["keff_diff"] = {'value':max(keff_list) - min(keff_list), 'output_index':10}

        if 11 in [parameters[key]['output_index'] for key in parameters.keys()]: 
            results_dict["cpr"] = {'value':calc_cpr(filename, parameters), 'output_index':11}
        if 12 in [parameters[key]['output_index'] for key in parameters.keys()]: 
            results_dict["lhgr"] = {'value':calc_lhgr(fq_list, parameters), 'output_index':12}
        if 13 in [parameters[key]['output_index'] for key in parameters.keys()]: 
            results_dict["aplhgr"] = {'value':calc_aplhgr(filename, parameters), 'output_index':13}
        results_dict["axoff_boc"] = {'value':axoff[0], 'output_index':14}
        results_dict["axoff_eoc"] = {'value':axoff[-1], 'output_index':15}
        results_dict["axoff_moc"] = {'value':axoff[int(len(axoff)/2)], 'output_index':16} #TODO parameterize this


        ## Correct Boron value if non-critical
        sorted_boron = sorted(boron_list,reverse=True)
        if sorted_boron[0] == sorted_boron[1]: #multiple values exceed maximum value in XS library, which is not reported by PARCS.
            new_max_boron = sorted_boron[0] #initialize variable
            for i in range(len(boron_list)): #!TODO: I think this serves to line up boron_list with keff_list. Could be replaced by index()
                if boron_list[i]== sorted_boron[0]:
                    boron_worth = 10.0 #pcm/ppm; assumed value.
                    excess_rho = (keff_list[i] - 1.0)*10**5 #pcm; excess reactivity
                    excess_boron = excess_rho/boron_worth #ppm
                    max_boron_corrected = sorted_boron[0] + excess_boron
                    if max_boron_corrected > new_max_boron:
                        new_max_boron = max_boron_corrected
            results_dict["max_boron"] = {'value':new_max_boron, 'output_index':14}
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["cycle_length"] = {'value':0.0, 'output_index':1}
        results_dict["pxy"] = {'value':10.0, 'output_index':2}
        results_dict["pxyz"] = {'value':10.0, 'output_index':3}
        results_dict["pinpowerpeaking"] = {'value':10.0, 'output_index':4}
        results_dict["fdeltah"] = {'value':10.0, 'output_index':5}
        results_dict["max_boron"]= {'value':10000, 'output_index':6}
        results_dict["chfr"] = {'value':0.0, 'output_index':7}
        results_dict["keff_min"] = {'value':0.0, 'output_index':8}
        results_dict["keff_max"] = {'value':10.0, 'output_index':9}
        results_dict["keff_diff"] = {'value':10.0, 'output_index':10}
        results_dict["cpr"] = {'value':0.0, 'output_index':11}
        results_dict["lhgr"] = {'value':100.0, 'output_index':12}
        results_dict["aplhgr"] = {'value':100.0, 'output_index':13}
        results_dict["axoff_boc"] = {'value':100, 'output_index':14}
        results_dict["axoff_eoc"] = {'value':-100, 'output_index':15}
        results_dict["axoff_moc"] = {'value':100, 'output_index':16}
    
    # populate parameters based on output index 
    for key in results_dict:
        output_index = results_dict[key]['output_index']
        for parameter_key in input.objectives:
            if 'output_index' in input.objectives[parameter_key].keys():
                if output_index == input.objectives[parameter_key]['output_index']:
                    parameters[parameter_key]["value"] = results_dict[key]['value']

    return parameters

def prepare_chromosome(input_obj, chromosome):
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
