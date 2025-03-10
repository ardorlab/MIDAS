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
from midas.utils.optimizer_tools import Constrain_Input
from midas_data import __polaris624exe__


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
        ofile.write("lib \"{}\"\n%\n".format("fine_n")) #!TODO: this value needs to be parameterized.
        ofile.write("system {}\n%\n".format("PWR")) #!TODO: this value needs to be parameterized.
    
    ## Geometry Block ##
    with open(filename,"a") as ofile:
        ofile.write("%\n% % % % % % % % % %\n% %   GEOMETRY  % %\n% % % % % % % % % %\n%\n")
        ofile.write(f"geom FA_1 :  ASSM {input.nrow} 1.26 sym={input.map_size}\n\n%\n")
        for name, pin in input.pin_options['rod_geometries'].items():
            ofile.write(f"pin {pin['type']} : {' '.join(map(str,pin['radii']))} : {' '.join(map(str,pin['materials']))} % {name}\n")
        for name, pin in input.pin_options['control_rods'].items():
            ofile.write(f"pin {pin['type']} : {' '.join(map(str,pin['radii']))} : {' '.join(map(str,pin['materials']))} % {name} in {pin['guide_tube']}\n")
        ofile.write("%\n")
        
        #!TODO: this only works for octant solution symmetry and needs to be moved to problem_preparation.py
        ofile.write("pinmap\n")
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
        else: # assume full lattice
            pass #!TODO: repeat above in all quadrants
        
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
        else: # assume full lattice
            pass #!TODO: repeat above in all quadrants
    
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
        ofile.write(f"state ALL : temp={input.bulk_temps} % K\n")
        for name, mat in input.pin_options['materials'].items():
            if mat['fueltype']:
                ofile.write(f"state {name.split('.')[0]} : temp={input.fuel_temps} % K\n")
        ofile.write(f"state BankA : in={str(input.cr_inserted).lower()} % insert rod cluster\n")
        if input.boronmat:
            ofile.write(f"state {input.boronmat[0]} : boron={input.boronmat[1]} % ppm\n")
        ofile.write("%\n")
        ofile.write(f"power {input.power} %W/gIHM\n%\n")
        ofile.write("deplete")
        for name, mat in input.pin_options['materials'].items():
            if mat['fueltype']:
                ofile.write(f" {name.split('.')[0]}=True")
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
        output = subprocess.check_output([polariscmd,filename,'>',''.join(filename.split('.')[:-1])+'.out'],\
                                            stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
    
    ## Get Results
        if 'Finished' in str(output): #job completed
            logger.debug(f"Job {solution.name} completed successfully in PARCSv343.")
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

def get_results(parameters, filename, job_failed=False):
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
        with open(''.join(filename.split('.')[:-1])+'.out', 'r') as fileread:
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
        for i in range(len(keff_list)):
            if keff_list[i] >= 1.0:
                continue
            else:
                results_dict["max_critical_exposure"]["value"] = exposure_list[i] #GWd/MTU
                subcritical = True
        if not subcritical: #!TODO: a more sophisticated way to handle this would be to extrapolate the rate of depletion to estimate the BU value at which k_eff crosses 1.0 (becomes subcritical).
            results_dict["max_critical_exposure"]["value"] = exposure_list[-1] #GWd/MTU
    
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["pinpowerpeaking"]["value"] = 10.0
        results_dict["peak_reactivity"]["value"] = 0.0
        results_dict["max_critical_exposure"]["value"] = 0.0 #GWd/MTU
    
    ## Save user-requested values from results
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
            if param not in ['cost_fuelcycle','av_fuelenrichment']: #check whitelist
                logger.warning(f"Parameter '{param}' not supported in POLARISv6.2.4 results parsing.")
    
    return parameters