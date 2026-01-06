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
from midas_data import __polaris624exe__;
import random #!needed for debug


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

## Functions ##
def evaluate(solution, input):
    """
    Interface used to run POLARISv6.2.4 calculations.
    
    Written by Nicholas Rollins. 03/06/2025
    """

## Create and move to unique directory for PARCS execution
    cwd = Path(os.getcwd())
    indv_dir = cwd.joinpath(input.results_dir_name / Path(solution.name))
    if not indv_dir.exists():
        logger.debug(f"Creating new results directory: {indv_dir}")
        os.mkdir(indv_dir)
    logger.debug(f"Changing to new working directory: {indv_dir}")
    os.chdir(indv_dir)
    
## Generate Input File
    filename = solution.name + '.inp'
    
    ## General Block ##
    with open(filename,"w") as ofile:
        ofile.write("=polaris\n")
        ofile.write("%\ntitle \"MIDAS Lattice Physics Calculation\"\n%\n")
        ofile.write("% Cross Section Library alias\n")
        ofile.write("lib \"{}\"\n%\n".format(input.xs_library)) 
        ofile.write("system {}\n%\n".format(input.system_type)) 
    
    ## Geometry Block ##
    with open(filename,"a") as ofile:
        ofile.write("%\n% % % % % % % % % %\n% %   GEOMETRY  % %\n% % % % % % % % % %\n%\n")
        if input.map_size == 'DIAGONAL':
            ofile.write(f"geom FA_1 :  ASSM {input.nrow} {input.pin_pitch} \n\n%\n")
        else:
            ofile.write(f"geom FA_1 :  ASSM {input.nrow} {input.pin_pitch} sym={input.map_size}\n\n%\n")
        if input.box != '0':
            ofile.write(f"box {input.box} \n")
        if input.hgap != '0':
            ofile.write(f"hgap {input.hgap} \n")        
        for name, pin in input.pin_options['rod_geometries'].items():
            ofile.write(f"pin {pin['type']} {int(np.ceil((pin['radii'][-1]*2 ) / input.pin_pitch))}: {' '.join(map(str,pin['radii']))} : {' '.join(map(str,pin['materials']))} % {name}\n")
        if 'control_rods' in input.pin_options.keys():
            for name, pin in input.pin_options['control_rods'].items():
                ofile.write(f"pin {pin['type']} : {' '.join(map(str,pin['radii']))} : {' '.join(map(str,pin['materials']))} % {name} in {pin['guide_tube']}\n")
        ofile.write("%\n")
        
        #!TODO: this only works for octant solution symmetry and needs to be moved to problem_preparation.py
        ofile.write("pinmap\n") # fuel pin map
        if input.map_size == "SE":
            half_size = int(np.ceil(input.nrow/2))
            for y in range(half_size):
                ofile.write("   ")
                for x in range(half_size):
                    if y < x:
                        index = int(x*(x+1)/2 + y)
                    else:
                        index = int(y*(y+1)/2 + x)
                    ofile.write(" " + str(input.pin_options['rod_geometries'][solution.chromosome[index]]['type']))
                ofile.write("\n")
            ofile.write("\n")
        elif input.map_size == "DIAGONAL":
            extended_chromosome = [0 for i in range(input.ncol * input.nrow)]
            for i in range(input.nrow):          
                start = i * (i + 1) // 2
                for c in range(i + 1): 
                    val = solution.chromosome[start + c]
                    extended_chromosome[i * input.nrow + c] = val
                    extended_chromosome[c * input.nrow + i] = val 

            index = 0
            for y in range(input.nrow):
                ofile.write("   ")
                for x in range(input.ncol):
                    ofile.write(" " + str(input.pin_options['rod_geometries'][extended_chromosome[index]]['type']))
                    index += 1
                ofile.write("\n")
            ofile.write("\n")
        else: # assume full lattice
            for y in range(int(input.ncol)):
                ofile.write("   ")
                for x in range(int(input.nrow)):
                    index = int(y*input.nrow + x)
                    ofile.write(" " + str(input.pin_options['rod_geometries'][solution.chromosome[index]]['type']))
                ofile.write("\n")
            ofile.write("\n")

        if 'control_rods' in input.pin_options.keys(): # control rod bank map
            ofile.write("control BankA : RODLET\n")
            if input.map_size == "SE":
                half_size = int(np.ceil(input.nrow/2))
                for y in range(half_size):
                    ofile.write("   ")
                    for x in range(half_size):
                        if y < x:
                            index = int(x*(x+1)/2 + y)
                        else:
                            index = int(y*(y+1)/2 + x)
                        _searching = True
                        for rod in input.pin_options['control_rods']:
                            if solution.chromosome[index] == input.pin_options['control_rods'][rod]['guide_tube']:
                                ofile.write(" " + str(input.pin_options['control_rods'][rod]['type']))
                                _searching = False
                                break
                        if _searching:
                            ofile.write(" _")
                    ofile.write("\n")
                ofile.write("\n%\n")
            elif input.map_size == "DIAGONAL":
                extended_chromosome = [0 for i in range(input.ncol * input.nrow)]
                for i in range(input.nrow):          
                    start = i * (i + 1) // 2
                    for c in range(i + 1): 
                        val = solution.chromosome[start + c]
                        extended_chromosome[i * input.nrow + c] = val
                        extended_chromosome[c * input.nrow + i] = val 

                index = 0
                for y in range(input.nrow):
                    ofile.write("   ")
                    for x in range(input.ncol):
                        _searching = True
                        for rod in input.pin_options['control_rods']:
                            if extended_chromosome[index] == input.pin_options['control_rods'][rod]['guide_tube']:
                                ofile.write(" " + str(input.pin_options['control_rods'][rod]['type']))
                                _searching = False
                                break
                        if _searching:
                            ofile.write(" _")
                        index += 1
                    ofile.write("\n")
                ofile.write("\n")
            else: # assume full lattice
                for y in range(int(input.ncol)):
                    ofile.write("   ")
                    for x in range(int(input.nrow)):
                        index = int(y*input.nrow + x)
                        _searching = True
                        for rod in input.pin_options['control_rods']:
                            if solution.chromosome[index] == input.pin_options['control_rods'][rod]['guide_tube']:
                                ofile.write(" " + str(input.pin_options['control_rods'][rod]['type']))
                                _searching = False
                                break
                        if _searching:
                            ofile.write(" _")
                    ofile.write("\n")
                ofile.write("\n%\n")
    
    ## Materials Block ##
    with open(filename,"a") as ofile:
        ofile.write("%\n% % % % % % % % % % %\n% %   Materials   % %\n% % % % % % % % % % %\n%\n")
        ofile.write("% Material compositions; referenced by mat cards\n")
        for name, comp in input.pin_options['compositions'].items():
            ofile.write(f"comp {name} : {comp['type']} {' '.join(map(str,comp['values']))}\n")
        ofile.write("%\n")
        if input.boronmat:
            ofile.write("% Set properties for material classes\n")
            ofile.write(f"prop boron {input.boronmat[0].split('.')[0]} : SOLP B scale=PPM\n")
            ofile.write("%\n")
        ofile.write("% Material Classes\n")
        for name, mat in input.pin_options['materials'].items():
            ofile.write(f"mat {name} : {mat['comp']} dens={mat['dens']}")
            if "temp" in mat:
                ofile.write(f" temp={mat['temp']}")
            ofile.write("\n")
    
    ## Statepoints Block ##
    with open(filename,"a") as ofile:
        ofile.write("%\n% % % % % % % % % % % %\n% %   STATEPOINTS   % %\n% % % % % % % % % % % %\n%\n")
        ofile.write("% Properties of stated materials\n")
        ofile.write(f"state ALL : temp={input.bulk_temps} ")
        for name, mat in input.pin_options['materials'].items():
            if mat['fueltype']:
                ofile.write(f" {name} : temp={input.fuel_temps} ")
        if 'control_rods' in input.pin_options.keys():
            ofile.write(f" BankA : in={str(input.cr_inserted).lower()} ")
        if input.boronmat:
            ofile.write(f" {input.boronmat[0]} : boron={input.boronmat[1]} \n")
        ofile.write("%\n")
        ofile.write(f"power {input.powdens} %W/gIHM\n%\n")
        ofile.write("deplete")
        for name, mat in input.pin_options['materials'].items():
            if mat['fueltype']:
                ofile.write(f" {name}=True")
        ofile.write("\n")
        ofile.write(f"bu")
        for step in input.depl_steps:
            ofile.write(' ' + str(step))
        ofile.write(" % burnup steps (GWd/MTU)\n%\n")
        ofile.write(f"opt GEOM MeshNumRings={input.num_meshrings} MeshNumSectors=1\n")
        
        ofile.write("end")
    
## Run POLARIS INPUT DECK #!TODO: separate the input writing and execution into two different functions that are called in sequence.
    polariscmd = __polaris624exe__
    try:
        completed_process = subprocess.run([polariscmd,filename], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
        calc_success = False
        with open(solution.name+'.msg','r')as ofile:
            for line in ofile:
                if 'Polaris execution completed with zero errors' in line:
                    calc_success = True

    ## Get Results
        if calc_success: #job completed
            logger.debug(f"Job {solution.name} completed successfully in POLARISv624.")
            solution.parameters = get_results(solution.parameters, solution.name)
        
        else: #job failed
            logger.warning(f"Job {solution.name} has failed!")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
    except subprocess.TimeoutExpired: #job timed out
        logger.error(f"Job {solution.name} has timed out!")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    except subprocess.CalledProcessError as e: #POLARIS returned an abort signal
        logger.error(f"Job {solution.name} has failed with the following exception: {e}")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
    logger.debug(f"Returning to original working directory: {cwd}")
    os.chdir(cwd)
    gc.collect()
    
    return solution

def get_results(parameters, solnname, job_failed=False):
    """
    Currently supports max Power Form Factor (pin power peaking), peak reactivity (k-inf),
    and max critical exposure (in GWd/MTU).
    
    Written by Nicholas Rollins. 03/10/2025
    """
    ## Prepare container for results
    results_dict = {}
    for res in ["pinpowerpeaking", "peak_reactivity", "max_critical_exposure"]:
        results_dict[res] = {}
        results_dict[res]['value'] = 0.0
    
    if not job_failed:
        exposure_list = []
        keff_list = []
        peak_power = 0.0
        ## Parse values at each statepoint
        with open(solnname+'.out', 'r') as fileread:
            valid = False
            for line in fileread:
                if not valid:
                    if "Initial burnup =" in line:
                        exposure_list.append(float(line.split('Initial burnup = ')[-1].split(',')[0]))
                        continue
                    elif "Transport: k-eff =  " in line:
                        keff_list.append(float(line.split('Transport: k-eff =  ')[-1].strip()))
                        continue
                    elif " Power Form Factors" in line and "Group" not in line:
                        valid = True
                        continue
                else:
                    if "-" in line:
                        continue #skip line
                    elif not line.strip():
                        valid = False
                        continue
                    else:
                        peak_power = max(peak_power,max([float(x) for x in line.split()]))
        
        results_dict["pinpowerpeaking"]["value"] = peak_power
        results_dict["peak_reactivity"]["value"] = max(keff_list)
        subcritical = False
        for i in range(1,len(keff_list)):
            if keff_list[i] >= 1.0:
                continue
            else:
                dbu_dkeff = (exposure_list[i] - exposure_list[i-1])/(keff_list[i] - keff_list[i-1])
                results_dict["max_critical_exposure"]["value"] = exposure_list[i-1] + dbu_dkeff*(1.0 - keff_list[i-1])#GWd/MTU; linear interpolation to keff = 1.0
                subcritical = True
                break
        if not subcritical:
            dbu_dkeff = (exposure_list[-1] - exposure_list[-2])/(keff_list[-1] - keff_list[-2])
            results_dict["max_critical_exposure"]["value"] = exposure_list[-2] + dbu_dkeff*(1.0 - keff_list[-2])#GWd/MTU; linear extrapolation to keff = 1.0
    
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["pinpowerpeaking"]["value"] = 10.0
        results_dict["peak_reactivity"]["value"] = 0.0
        results_dict["max_critical_exposure"]["value"] = 0.0 #GWd/MTU
    
    ## Save user-requested values from results
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
        else:
            if param not in ['cost_fuelcycle','av_fuelenrichment']: #check whitelist
                logger.warning(f"Parameter '{param}' not supported in POLARISv6.2.4 results parsing.")
    
    return parameters