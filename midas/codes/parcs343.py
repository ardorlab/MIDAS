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
from midas_data import __parcs343exe__


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

## Functions ##
def evaluate(solution, input):
    """
    #!TODO: write docstring.
    
    Updated by Nicholas Rollins. 10/03/2024
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
    
## Prepare values for file writing
    list_unique_xs = np.concatenate([value if isinstance(value,list) else np.concatenate(list(value.values()))\
                                    for value in input.xs_list.values()])

    ## Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)
    fuel_locations = [loc for loc in input.core_dict.keys() if 2 < len(loc) <  5]
    soln_fuel_locations = {}
    if input.calculation_type in ['eq_cycle']:
        soln_FAs = Constrain_Input.SS_decoder(solution.chromosome)
        for i in range(len(solution.chromosome)):
            soln_fuel_locations[fuel_locations[i]] = soln_FAs[i]
    else:
        for i in range(len(solution.chromosome)):
            soln_fuel_locations[fuel_locations[i]] = solution.chromosome[i]
    
    soln_core_dict = deepcopy(input.core_dict)
    for loc, label in soln_fuel_locations.items():
        tag = None
        for fueltype in input.tag_list['fuel']:
            if fueltype[1] == label:
                tag = fueltype[0]
        if not tag:
            raise ValueError("FA label not found in tag_list.")
        soln_core_dict[loc]['Value'] = tag
    #!for loc, label in soln_refl_locations.items(): #!TODO: create a way to specify reflector locs for multiple radial refls.

    soln_core_lattice = deepcopy(input.core_lattice) # core lattice filled with chromosome
    for loc, vals in soln_core_dict.items():
        sym_locs = [loc] + vals['Symmetric_Assemblies']
        for j in range(len(soln_core_lattice)):
            for i in range(len(soln_core_lattice[j])):
                if soln_core_lattice[j][i] in sym_locs:
                    if soln_core_lattice[j][i][0] == "R" and len(soln_core_lattice[j][i]) >= 5: #reflector
                        soln_core_lattice[j][i] = "10" #!TODO: add support more multiple radial refls.
                    else:
                        soln_core_lattice[j][i] = vals['Value']
    
## Generate Input File
    filename = solution.name + '.inp'
    
    ## CaseID Block ##
    with open(filename,"w") as ofile:
        ofile.write("!******************************************************************************\n")
        ofile.write('CASEID {}  \n'.format(solution.name))
        ofile.write("!******************************************************************************\n\n")

    ## CNTL Block ##
    with open(filename,"a") as ofile:
        ofile.write("CNTL\n")
        ofile.write("      RUN_OPTS   F T F F\n")
        if input.th_fdbk:
            ofile.write("      TH_FDBK    T\n")
            ofile.write("      INT_TH     T -1\n")
        else:
            ofile.write("      TH_FDBK    F\n")
        ofile.write("      CORE_POWER 100.0\n")
        ofile.write("      CORE_TYPE  PWR\n")
        ofile.write("      PPM        1000 1.0 1800.0 10.0\n") #!TODO: this should be a parameterized boron guess value.
        ofile.write("      DEPLETION  T  1.0E-5 T\n")
        if input.calculation_type in ['eq_cycle']:
            ofile.write("      MULT_CYC   T  F\n") #v3.4.2 specific line to enable the MCYCLE block
        ofile.write("      TREE_XS    T  {}  T  T  F  F  T  F  F  F  T  F  T  T  T  F  F \n".format(int(len(list_unique_xs))))
        ofile.write("      BANK_POS   100 100 100 100 100 100\n")
        ofile.write("      XE_SM      1 1 1 1\n")
        ofile.write("      SEARCH     PPM\n")
        ofile.write("      XS_EXTRAP  1.0 0.3\n")
        if input.pin_power_recon:
            ofile.write("      PIN_POWER  T\n")
        else:
            ofile.write("      PIN_POWER  F\n")
        ofile.write("      PRINT_OPT  T T T T T F T T T T  T  T  T  T  F  T  T")
        #!ofile.write("      PLOT_OPTS 0 0 0 0 0 2\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        
    ## PARAM Block ##
    with open(filename,"a") as ofile:
        ofile.write("PARAM\n")
        ofile.write("      LSOLVER     1 1 20\n")
        ofile.write("      NODAL_KERN  NEMMG\n")
        ofile.write("      CMFD        2\n")
        ofile.write("      DECUSP      2\n")
        ofile.write("      INIT_GUESS  0\n")
        ofile.write("      CONV_SS     1.e-6 5.e-5 1.e-3 0.001\n")
        ofile.write("      EPS_ERF     0.010\n")
        ofile.write("      EPS_ANM     0.000001\n")
        ofile.write("      NLUPD_SS    5 5 1\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
    
    ## GEOM Block Inputs ##
    with open(filename,"a") as ofile:
        ofile.write("GEOM\n")
        if input.map_size == 'quarter':
            dim_size = [np.floor(input.nrow/2)+1, np.floor(input.ncol/2)+1]
        else: #assume full geometry if not quarter-core
            dim_size = [input.nrow, input.ncol]
        ofile.write(f"      GEO_DIM {dim_size[0]} {dim_size[1]} {input.number_axial} 1 1\n")
        ofile.write("      RAD_CONF\n\n")
        for x in range(soln_core_lattice.shape[0]):
            ofile.write("      ")
            for y in range(soln_core_lattice.shape[1]):
                ofile.write(soln_core_lattice[x,y])
                ofile.write("  ")
            ofile.write("\n")
        ofile.write("\n")
    
        assembly_width = 21.50 #!TODO: change this to an input with default.
        if input.map_size == 'quarter':
            ofile.write(f"      GRID_X      1*{assembly_width/2} {dim_size[0]-1}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  1*1 {dim_size[0]-1}*1\n")
            ofile.write(f"      GRID_Y      1*{assembly_width/2} {dim_size[0]-1}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  1*1 {dim_size[0]-1}*1\n")
        else: #assume full geometry if not quarter-core
            ofile.write(f"      GRID_X      {dim_size[0]}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  {dim_size[0]}*1\n")
            ofile.write(f"      GRID_Y      {dim_size[1]}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  {dim_size[1]}*1\n")
        ofile.write("      GRID_Z      {}\n".format('  '.join([str(x) for x in input.axial_nodes])))
        # Write radial reflectors
        xsnum_radtop = 2 + len(input.xs_list['reflectors']['radial'])
        rad_tags = [tag[0] for tag in input.tag_list['reflectors']]
        for i in range(len(input.xs_list['reflectors']['radial'])):
            tag = input.tag_list['reflectors'][rad_tags.index(input.tag_list['reflectors'][i][0])][0]
            ofile.write("      ASSY_TYPE   {}   1*1  {}*{}  1*{} REFL\n".format(tag,input.number_axial-2,2+i,xsnum_radtop))
        # Write fuel types
        if 'blankets' in input.fa_options:
            xsnum_fuel = xsnum_radtop + len(input.xs_list['blankets'])
        else:
            xsnum_fuel = xsnum_radtop
        for key in input.fa_options['fuel'].keys():
            fuel = input.fa_options['fuel'][key]
            xsnum_fuel += 1
            if 'blanket' in fuel:
                xsnum_blanket = xsnum_radtop + \
                                input.xs_list['blankets'].index(input.fa_options['blankets'][fuel['blanket']]['serial']) + 1
                ofile.write("      ASSY_TYPE   {}   1*1  1*{} {}*{}  1*{}  1*{} FUEL\n".format(fuel['type'],xsnum_blanket,\
                                                                                       input.number_axial-4,xsnum_fuel,\
                                                                                       xsnum_blanket,xsnum_radtop))
            else:
                ofile.write("      ASSY_TYPE   {}   1*1  {}*{}  1*{} FUEL\n".format(fuel['type'],input.number_axial-2,\
                                                                                  xsnum_fuel,xsnum_radtop))
        ofile.write("\n")

        if input.map_size == 'quarter':
            ofile.write("      BOUN_COND   0 2 0 2 2 2\n")
            ofile.write("      SYMMETRY 4\n")
        else: #assume full geometry if not quarter-core
            ofile.write("      BOUN_COND   2 2 2 2 2 2\n")
            ofile.write("      SYMMETRY 1\n")

        ofile.write("    PINCAL_LOC\n")
        for x in range(input.pincal_loc.shape[0]):
            ofile.write("      ")
            for y in range(input.pincal_loc.shape[1]):
                val = input.pincal_loc[x,y]
                try:
                    if not np.isnan(val):
                        ofile.write(str(input.pincal_loc[x,y]))
                        ofile.write("  ")
                except TypeError:
                    ofile.write(str(input.pincal_loc[x,y]))
                    ofile.write("  ")
            ofile.write("\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## FDBK Block ##
    with open(filename,"a") as ofile:
        ofile.write("FDBK\n")
        ofile.write("      FA_POWPIT     {} {}\n".format(np.round(input.power/input.num_assemblies,4),assembly_width))
        ofile.write("      GAMMA_FRAC    0.0208    0.0    0.0\n")
        ofile.write("      EFF_DOPLT   T  0.5556\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## TH Block ##
    with open(filename,"a") as ofile:
        ofile.write("TH\n")
        if input.th_fdbk:
            ofile.write("      FLU_TYP       0\n")
            ofile.write("      N_PINGT    264 25\n")
            ofile.write("      PIN_DIM      4.1 4.75 0.58 6.13\n")
            ofile.write("      FLOW_COND    {}  {}\n".format(np.round(input.inlet_temp-273.15,2),\
                                                             np.round(input.flow/input.num_assemblies,4)))
            ofile.write("      HGAP     11356.0\n") #!TODO:check this value, should it be parameterized?
            ofile.write("      N_RING   6\n")
            ofile.write("      THMESH_X       9*1\n")
            ofile.write("      THMESH_Y       9*1\n")
            ofile.write("      THMESH_Z       1 2 3 4 5 6 7 8 9 10 11 12\n")
        else:
            ofile.write("      UNIF_TH   0.740    626.85     {}\n".format(np.round(input.inlet_temp-273.15,2))) #!TODO: how to deal with av. fuel temp?
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## DEPL Block ##
    with open(filename,"a") as ofile:
        ofile.write("DEPL\n")
        if input.calculation_type == 'single_cycle':
            ofile.write(f"      TIME_STP  {str(input.depl_steps).strip('[]')}\n")
        ofile.write("      INP_HST   './boc_exp.dep' -2 1\n")
        ofile.write("      OUT_OPT   T  T  T  T  F\n")
        # Write reflector cross sections
        ofile.write("      PMAXS_F   1 '{}{}' 1\n".format(input.xs_lib / Path(input.xs_list['reflectors']['bot'][0]),\
                                                        input.xs_extension))
        for i in range(len(input.xs_list['reflectors']['radial'])):
            rxs_index = 2 + i
            radpath = input.xs_lib / Path(input.xs_list['reflectors']['radial'][i])
            ofile.write("      PMAXS_F   {} '{}{}' {}\n".format(rxs_index,radpath,input.xs_extension,rxs_index))
        ofile.write("      PMAXS_F   {} '{}{}' {}\n".format(rxs_index+1,\
                                                          input.xs_lib / Path(input.xs_list['reflectors']['top'][0]),\
                                                          input.xs_extension,rxs_index+1))
        nxs_index = rxs_index + 2
        # Write blankets cross sections
        if 'blankets' in input.fa_options:
            for i in range(len(input.xs_list['blankets'])):
                bxs_index = i + rxs_index + 2
                blanketpath = input.xs_lib / Path(input.xs_list['blankets'][i])
                ofile.write("      PMAXS_F   {} '{}{}' {}\n".format(bxs_index,blanketpath,input.xs_extension,bxs_index))
            nxs_index = bxs_index + 1
            
        # Write fuel types cross sections
        for i in range(len(input.xs_list['fuel'])):
            fxs_index = i + nxs_index
            ofile.write("      PMAXS_F   {} '{}{}' {}\n".format(fxs_index,\
                                                              input.xs_lib / Path(input.xs_list['fuel'][i]),\
                                                              input.xs_extension,fxs_index))
    
    ## MCYCLE Block ##
    if input.calculation_type in ['eq_cycle']:
        soln_full_core_lattice = prepare_shuffling_map(input, solution.chromosome)
        with open(filename,"a") as ofile:
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
            
            ofile.write("MCYCLE\n")
            ofile.write("    CYCLE_DEF   1\n")
            ofile.write(f"      DEPL_STEP {str(input.depl_steps).strip('[]')}\n")
            ofile.write(f"      POWER_LEV {len(input.depl_steps)+1}*100.0\n")
            ofile.write(f"      BANK_SEQ  {len(input.depl_steps)+1}*1\n\n")
            
            ofile.write("    LOCATION   0\n")
            for x in range(input.full_core_locs.shape[0]):
                for y in range(input.full_core_locs.shape[1]):
                    val = input.full_core_locs[x,y]
                    try:
                        if not np.isnan(val):
                            ofile.write(str(input.full_core_locs[x,y]))
                            ofile.write("  ")
                    except TypeError:
                        ofile.write(str(input.full_core_locs[x,y]))
                        ofile.write("  ")
                ofile.write("\n")
            ofile.write("\n")
            
            ofile.write("    SHUF_MAP   1   1\n")
            for x in range(soln_full_core_lattice.shape[0]):
                for y in range(soln_full_core_lattice.shape[1]):
                    val = soln_full_core_lattice[x,y]
                    ofile.write(str(soln_full_core_lattice[x,y]))
                    ofile.write("  ")
                ofile.write('\n')
            ofile.write('\n')
            
            ofile.write("    CYCLE_IND    1  0  1\n")
            max_convergence_cycles = 10 #!TODO: this max number of cycles could easily be a parameter.
            for i in range(2,max_convergence_cycles+1):
                ofile.write(f"    CYCLE_IND    {i}  1  1\n")
            ofile.write(f"    CONV_EC    0.1  {i}\n")
    
    ## Termination Character ##
    with open(filename,"a") as ofile:
        ofile.write(".")

## Run PARCS INPUT DECK #!TODO: separate the input writing and execution into two different functions that are called in sequence.
    parcscmd = __parcs343exe__
    
    if input.calculation_type in ['eq_cycle']:
        walltime = 1800 #sec
    else:
        walltime = 600 #sec
    try:
        output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=walltime) #wait until calculation finishes
    ## Get Results
        if 'Finished' in str(output): #job completed
            logger.debug(f"Job {solution.name} completed successfully.")
            solution.parameters = get_results(solution.parameters, solution.name)
        
        else: #job failed
            if input.calculation_type in ['eq_cycle']:
                solution.parameters = eq_cycle_convergenc(input, solution, filename, parcscmd, walltime) #iteratively try to find an intial guess that will converge
            else: #standard execution pathway
                logger.warning(f"Job {solution.name} has failed!")
                solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
    except subprocess.TimeoutExpired: #job timed out
        os.system('rm -f {}.parcs_pin*'.format(solution.name))
        logger.error(f"Job {solution.name} has timed out!")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
    logger.debug(f"Returning to original working directory: {cwd}")
    os.chdir(cwd)
    gc.collect()
    
    return solution

def get_results(parameters, filename, job_failed=False): #!TODO: implement pin power reconstruction.
    """
    Currently supports cycle length, F_q, F_dh, and max boron.
    
    Updated by Nicholas Rollins. 09/27/2024
    """
    ## Prepare container for results
    results_dict = {}
    for res in ["cycle_length", "pinpowerpeaking", "fdeltah", "max_boron"]:
        results_dict[res] = {}
        results_dict[res]['value'] = []
        
    if not job_failed:
        ## Read file for parsing
        with open(filename + ".parcs_dpl", "r") as ofile:
            filestr = ofile.read()
        
        ## Split file by section
        res_str = filestr.split('===============================================================================')
        res_str = res_str[-1].split('_______________________________________________________________________________')
        res_str = res_str[0].split('\n')
        
        ## Parse raw values by timestep
        efpd_list = []; boron_list = []; keff_list = []; fq_list = []; fdh_list = []
        for i in range(2, len(res_str)-1):
            res_val=res_str[i].split()
            
            efpd_list.append(float(res_val[9]))
            boron_list.append(float(res_val[14]))
            keff_list.append(float(res_val[2]))
            fq_list.append(float(res_val[7]))
            fdh_list.append(float(res_val[6]))
        
        del filestr, res_str, res_val #unload file contents to clean up memory
        
        results_dict["cycle_length"]["value"] = calc_cycle_length(efpd_list,boron_list,keff_list)
        results_dict["pinpowerpeaking"]["value"] = max(fq_list)
        results_dict["fdeltah"]["value"] = max(fdh_list)
        results_dict["max_boron"]["value"] = max(boron_list)
        
        ## Correct Boron value if non-critical
        if results_dict["max_boron"]["value"] == 1800.0: #!TODO: initial guess should be a variable. can this be read from output file?
            new_max_boron = 0
            for i in range(len(boron_list)): #!TODO: I think this serves to line up boron_list with keff_list. Could be replaced by index()
                if boron_list[i]== 1800.0:
                    boron_worth = 10.0 #pcm/ppm
                    excess_rho = (keff_list[i] - 1.0)*10**5 #pcm; excess reactivity
                    excess_boron = excess_rho/boron_worth #ppm
                    max_boron_corrected = 1800.0 + excess_boron
                    if max_boron_corrected > new_max_boron:
                        new_max_boron = max_boron_corrected
            results_dict["max_boron"]["value"] = new_max_boron
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["cycle_length"]["value"] = 0.0
        results_dict["pinpowerpeaking"]["value"] = 10.0
        results_dict["fdeltah"]["value"] = 10.0
        results_dict["max_boron"]["value"] = 10000
    
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
        else:
            if param not in ['cost_fuelcycle']: #check whitelist
                logger.warning(f"Parameter '{param}' not supported in PARCS343 results parsing.")
    
    return parameters

def calc_cycle_length(efpd,boron,keff):
    if boron[-1]==0.1: #boron went to zero before end of cycle.
        eoc1_ind = 0
        eco2_ind = len(efpd)
        for i in range(len(efpd)):
            if boron[i] > 0.1 and boron[i+1] == 0.1:
                eoc1_ind = i
                eco2_ind = i+1
        dbor = abs(boron[eoc1_ind-1]-boron[eoc1_ind])
        defpd = abs(efpd[eoc1_ind-1]-efpd[eoc1_ind])
        def_dbor = defpd/dbor
        eoc = efpd[eoc1_ind] + def_dbor*(boron[eoc1_ind]-0.1) #linear extrapolation to efpd at boron=0.1
    elif boron[-1]==boron[0]==1800.0: #true boron exceeds initial guess #!TODO: this should be a parameterized boron guess value.
        drho_dcb=10
        drho1 = (keff[-2]-1.0)*10**5
        dcb1 = drho1/drho_dcb
        cb1= boron[-2] + dcb1
        drho2 = (keff[-1]-1.0)*10**5
        dcb2 = drho2/drho_dcb
        cb2= boron[-1] + dcb2
        dbor = abs(cb1-cb2)
        defpd = abs(efpd[-2]-efpd[-1])
        def_dbor = defpd/dbor
        eoc = efpd[-1] + def_dbor*(cb2-0.1)
    else:
        dbor = abs(boron[-2]-boron[-1])
        defpd = abs(efpd[-2]-efpd[-1])
        def_dbor = defpd/dbor
        eoc = efpd[-1] + def_dbor*(boron[-1]-0.1)
    return eoc

def prepare_shuffling_map(input, chromosome):
    """
    Prepares the formatted full-core shuffling scheme for the MCYCL block.
    
    Written by Nicholas Rollins. 10/16/2024
    """
    # map with labels that much the format in input.core_dict
    full_core_lattice = deepcopy(input.full_core_locs)
    for i in range(len(full_core_lattice)):
        for j in range(len(full_core_lattice[i])):
            full_core_lattice[i][j] = full_core_lattice[i][j].replace('-','')

    # array to hold the final map for printing
    soln_full_core_lattice = deepcopy(full_core_lattice)
    soln_labels_list = []
    for key in input.core_dict.keys():
        if key[0] != "R":
            soln_labels_list.append(key)        

    # fill map with shuffling scheme
    for i in range(len(chromosome)):
        labels = [soln_labels_list[i]] #which loc/FA we're pulling from
        labels.extend(input.core_dict[labels[0]]['Symmetric_Assemblies'])
        try:
            source_labels = [soln_labels_list[chromosome[i][1]]]
            source_labels.extend(input.core_dict[source_labels[0]]['Symmetric_Assemblies'])
            tags = []
            for loc in source_labels:
                tags.append(loc[0]+'-'+loc[1:]) #value to print in loc
        except TypeError:
            for fueltype in input.tag_list['fuel']:
                if fueltype[1] == chromosome[i][1]:
                    tags = [' -'+fueltype[0]]*len(labels)
                    break

        for idx in range(len(labels)):
            for irow in range(len(full_core_lattice)):
                for jcol in range(len(full_core_lattice[irow])):
                    if full_core_lattice[irow][jcol] == labels[idx]:
                        soln_full_core_lattice[irow][jcol] = tags[idx]
    
    return soln_full_core_lattice

def eq_cycle_convergenc(input, solution, filename, parcscmd, walltime):
    boc_exp = input.boc_exposure
    conv_list = [[],[]] #track convergence
    skip_convwrite = False
## fetch best cycle from previous attempt
    depfiles_list = []
    for file in os.listdir('./'):
        if '.parcs_cyc-' in file:
            depfiles_list.append(file)
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
        logger.debug(f"Job {solution.name} completed successfully.")
        solution.parameters = get_results(solution.parameters, solution.name)
    else:
        logger.warning(f"Job {solution.name} has failed!")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
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