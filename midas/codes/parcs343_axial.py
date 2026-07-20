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

def with_template(solution, input, cwd, filename): 

## Prepare values for file writing
    ## Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)

## copy input file from template
    inp_template = str(cwd.joinpath(cwd / input.input_template['loc']))
    shutil.copy(inp_template, filename)

    with open(filename, "r") as file:
        lines = file.readlines()  

    with open(filename, "w") as ofile:
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
                for i in range(int(len(solution.chromosome)/4)):
                    line = f"    assy_type    {assmb_type}  1*1  1*56  8*{solution.chromosome[(i*4)]}  5*{solution.chromosome[(i*4)+1]}  2*{solution.chromosome[(i*4)+2]}  6*{solution.chromosome[(i*4)+3]}  1*57  1*58 1*2  FUEL \n"
                    ofile.write(line) 
                    assmb_type += 10
            else:
                ofile.write(line)   


def get_results(parameters, filename, input, job_failed=False):
    """
    Currently supports cycle length, F_q, F_dh, max boron, keff, critical power ratio,
    linear heat generation rate, average planar linear heat generation rate.
    
    Updated by Nicholas Rollins. 09/27/2024
    Updated by Jake Mikouchi. 03/25/2025
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
        efpd_list = []; boron_list = []; keff_list = []; pxy_list = []; pxyz_list = []; fq_list = []; fdh_list = []; chfr = []
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
        
        del filestr, res_str, res_val #unload file contents to clean up memory
        
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
    
    # populate parameters based on output index 
    for key in results_dict:
        output_index = results_dict[key]['output_index']
        for parameter_key in input.objectives:
            if 'output_index' in input.objectives[parameter_key].keys():
                if output_index == input.objectives[parameter_key]['output_index']:
                    parameters[parameter_key]["value"] = results_dict[key]['value']

    return parameters

def calc_cycle_length(efpd,boron,keff,extrapolate):
    if not extrapolate:
        eoc1_ind = 0
        if boron[-1]==0.1: #boron went to zero before end of cycle.
            for i in range(len(efpd)):
                if boron[i] > 0.1 and boron[i+1] == 0.1:
                    eoc1_ind = i
                    break
        else: 
            eoc1_ind = len(efpd)-1

        eoc = efpd[eoc1_ind]
    else:
        if boron[-1]==0.1: #boron went to zero before end of cycle.
            eoc1_ind = 0
            eco2_ind = len(efpd)-1
            for i in range(len(efpd)):
                if boron[i] > 0.1 and boron[i+1] == 0.1:
                    eoc1_ind = i
                    eco2_ind = i+1
                    break
            if eoc1_ind != 0:
                dbor = abs(boron[eoc1_ind]-boron[eoc1_ind-1])
                defpd = abs(efpd[eoc1_ind]-efpd[eoc1_ind-1])
            else:
                dbor = abs(boron[eco2_ind]-boron[eoc1_ind])
                defpd = abs(efpd[eco2_ind]-efpd[eoc1_ind])
            try:
                def_dbor = defpd/dbor
            except ZeroDivisionError:
                def_dbor = 0.0
            eoc = efpd[eoc1_ind] + def_dbor*(boron[eoc1_ind]-boron[eco2_ind]) #linear extrapolation to efpd at boron=0.1
        elif boron[-1]==boron[0]: #true boron exceeds initial guess
            drho_dcb=10 #pcm/ppm
            drho1 = (keff[-2]-1.0)*10**5 #pcm
            cb1= boron[-2] + drho1/drho_dcb #corrected boron concentration
            drho2 = (keff[-1]-1.0)*10**5 #pcm
            cb2= boron[-1] + drho2/drho_dcb #corrected boron concentration
            dbor = abs(cb1-cb2) #ppm
            defpd = abs(efpd[-2]-efpd[-1]) #efpd
            def_dbor = defpd/dbor #efpd/ppm
            eoc = efpd[-1] + def_dbor*(cb2-0.1)
        else: #EOC boron is greater than 0.1
            dbor = abs(boron[-2]-boron[-1])
            defpd = abs(efpd[-2]-efpd[-1])
            def_dbor = defpd/dbor #slope
            eoc = efpd[-1] + def_dbor*(boron[-1]-0.1) #linear extrapolation
    return eoc

def calc_cpr(filename, parameters): #TODO update this to use CHFR output by parcs
    # maximum relative power fraction at each depletion step
    with open(f"{filename}.parcs_out","r") as f:
        lines = f.readlines()
        active_read = False
        peak_assembly_power = []
        for line in lines:
            if " assembly power distribution" in line.lower():
                active_read = True

            if "(" in line.lower() and active_read:
                peak_assembly_power.append(float(line.split()[-1]))
                active_read = False

    # single highest assembly power throughout cycle
    max_ap = max(peak_assembly_power) 

    #determines critical power ratio
    try: 
        cpr =  parameters['cpr']['critical_power'] / max_ap 
    except:
        cpr = 1 / max_ap

    return cpr

def calc_lhgr(fq_list, parameters):
    # max pin power
    max_pp = max(fq_list)

    #determines maximum linear heat generation rate
    lhgr = parameters['lhgr']['linear_power'] * max_pp

    return lhgr

def calc_aplhgr(filename, parameters):

    direct = os.getcwd()
    peakval = 0
    # retrieve average planar lhgr from each assembly
    for i  in os.listdir(direct):
        if os.path.isfile(os.path.join(direct,i)) and "parcs_pin" in i:
            with open(i, "r") as pin_file:
                lines = pin_file.readlines()
                fq = 0
                step = []
                counts = 0
                read = False
                for line in lines:
                    if read:
                        step.append(line.split()[1:])

                    if "Case:" in line and counts == 2:
                        read = True
                        counts = 0  

                    if "Case:" in line and float(line.split()[-1]) == 0.0:
                        try:
                            if counts == 0 and step:
                                pin_count = 1
                                for j in range(len(step[1:-1])):
                                    for k in step[j]:
                                        if float(k) > 0.0:
                                            fq += float(k)
                                            pin_count += 1 

                                fq = fq / pin_count
                                if fq > peakval:
                                    peakval = fq

                        except: 
                            pass

                        step = []
                        read = False
                        counts += 1 

    aplhgr = parameters['aplhgr']['linear_power'] * peakval

    return aplhgr

def eq_cycle_convergence(input, solution, filename, parcscmd, walltime):
    boc_exp = input.boc_exposure
    conv_list = [[],[]] #track convergence
    skip_convwrite = False
## fetch best cycle from previous attempt
    depfiles_list = []
    for file in os.listdir('./'):
        if '.parcs_cyc-' in file:
            depfiles_list.append(file)
    if not depfiles_list:
        raise ValueError("PARCS failed without generating any depletion results.")
    if os.path.getsize(depfiles_list[-1]) < 20000: #if file is too small the cycle didn't initialize
        lastcycle_dep = depfiles_list[-2]
    else:
        lastcycle_dep = depfiles_list[-1]
## reattempt convergence
    convergence_attempts = 0
    while convergence_attempts < 8: # number of attempts to make
        convergence_attempts += 1
    ## fetch best cycle from previous attempt
        depfiles_list = []
        for file in os.listdir('./'):
            if '.parcs_cyc-' in file:
                depfiles_list.append(file)
        if os.path.getsize(depfiles_list[-1]) < 20000: #if file is too small the cycle didn't initialize
            lastcycle_dep = depfiles_list[-2]
            os.replace(depfiles_list[-2],'restart_exp.dep')
        else:
            lastcycle_dep = depfiles_list[-1]
            os.replace(depfiles_list[-1],'restart_exp.dep')
        os.system("rm *.parcs_cyc-*") #delete old cycle results in case they don't get overwritten
    ## Update convergence tracking list
        if not skip_convwrite:
            conv_list[0].append(float(boc_exp))
            conv_list[1].append(int(lastcycle_dep[-2:]))
        else:
            skip_convwrite = False
    ## decide restart pathway
        if int(lastcycle_dep[-2:]) >= 5:
            skip_convwrite = True
            #restart from the new file
            logger.debug(f"Job {solution.name} has failed to converge {convergence_attempts} time(s). Retrying from cycle {int(lastcycle_dep[-2:])}...")
        ## edit inp file
            with open(filename, 'r+') as file:
                lines = file.readlines()  # Read all lines
                for i, line in enumerate(lines):
                    if line.strip().startswith("INP_HST"):
                        line = line.replace("'./boc_exp.dep' -2","'./restart_exp.dep' 1")
                        lines[i] = line # Update the line in the list
                # Move back to the start of the file and truncate to overwrite
                file.seek(0)
                file.writelines(lines)
                file.truncate()  # Ensures any remaining old content is removed if file size decreases
        else:
            #restart with new boc exposure
            boc_exp = next_binary_search(conv_list) #try a new boc exposure
            logger.debug(f"Job {solution.name} has failed to converge {convergence_attempts} time(s). Retrying with a BOC exposure of {boc_exp} GWd/MTU...")
        ## rewrite boc_exp.dep file
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
                                depfile.write('{:.3f}'.format(boc_exp).rjust(8))
                            depfile.write('\n')
                        depfile.write('\n')
                        columncount = 0
                ## write rows for leftover columns
                if columncount!= 0:
                    depfile.write('\n')
                    for j in range(input.number_axial-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
                        depfile.write(' '+str(j).ljust(3))
                        for k in range(columncount):
                            depfile.write('{:.3f}'.format(boc_exp).rjust(8))
                        depfile.write('\n')
                    depfile.write('\n')
                depfile.write(' END STEP\n')
        ## edit inp file
            with open(filename, 'r+') as file:
                lines = file.readlines()  # Read all lines
                for i, line in enumerate(lines):
                    if line.strip().startswith("INP_HST"):
                        line = line.replace("'./restart_exp.dep' 1","'./boc_exp.dep' -2")
                        lines[i] = line # Update the line in the list
                # Move back to the start of the file and truncate to overwrite
                file.seek(0)
                file.writelines(lines)
                file.truncate()  # Ensures any remaining old content is removed if file size decreases
        
        #try again with new starting point
        output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=walltime) #wait until calculation finishes
    if 'Finished' in str(output): #job completed
        logger.debug(f"Job {solution.name} completed successfully in PARCSv343.")
        solution.parameters = get_results(solution.parameters, solution.name, input)
    else:
        logger.warning(f"Job {solution.name} has failed!")
        solution.parameters = get_results(solution.parameters, solution.name, input, job_failed=True)
    
    return solution.parameters

def next_binary_search(search_list):
    """
    Find the next value to try in a binary search.
    
    search_list = list; should contain 2 arrays, the first continuing the already tried
                x values and the second should try the resulting y values.
    
    Written by Nicholas Rollins. 11/11/2024
    """
    #lower and upper bounds
    lower_bound = 0.0
    upper_bound = 25.0
    
    # Sort the lists by the elements of the first list
    x_list, y_list = zip(*sorted(list(zip(search_list[0], search_list[1])), key=lambda a: a[0]))
    if int(lower_bound) not in [int(x) for x in x_list]:
        return lower_bound
    elif int(upper_bound) not in [int(x) for x in x_list]:
        return upper_bound

    best_index = y_list.index(max(y_list))
    x1 = x_list[best_index] # get best x from index of best y
    if best_index == 0:
        x2 = x_list[best_index+1]
    elif best_index == len(x_list)-1:
        x2 = x_list[best_index-1]
    else:
        x2 = x_list[best_index+1] if y_list[best_index+1] > y_list[best_index-1] else x_list[best_index-1]

    if x1 >= x2:
        return float(x1 - (x1-x2)/2)
    else:
        return float(x2 - (x2-x1)/2)