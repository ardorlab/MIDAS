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
## for surrogate model
import midas_data 
from pathlib import Path


### import surrogate model here 
import itertools
import pickle 
import joblib
import copy
import time as TT
import matplotlib.pyplot as plt
from surmodel.paralel_MLmodel import get_result
from surmodel.paralel_MLmodel import get_result_157


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input,xsdict):
    """
    Interface used to run surrogate model calculations.
    
    evaluate function creates working directory and prepares depletion file.
    if a parcs input file template is provided by the user in the yaml file, the parcs input will be created in with_template()
    if no template is provided (base case) then the files will be created using without_template()

    Updated by Nicholas Rollins. 10/03/2024
    Updated by Jake Mikouchi. 03/18/2025
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
    # with open('boc_exp.dep',"w") as depfile:
    #     depfile.write("\n BEGIN STEP\n\n EXP 3D MAP 1.0E+00\n\n")
    #     columncount = 0
    #     for i in range(1,input.num_assemblies+1):
    #         ## write column headers
    #         if columncount == 0:
    #             depfile.write(" k lb ")
    #         depfile.write(str(i).ljust(8))
    #         columncount += 1
    #         ## write rows for every 10 columns
    #         if columncount == 10:
    #             depfile.write('\n')
    #             for j in range(input.number_axial-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
    #                 depfile.write(' '+str(j).ljust(3))
    #                 for k in range(columncount):
    #                     depfile.write('{:.3f}'.format(input.boc_exposure).rjust(8))
    #                 depfile.write('\n')
    #             depfile.write('\n')
    #             columncount = 0
    #     ## write rows for leftover columns
    #     if columncount!= 0:
    #         depfile.write('\n')
    #         for j in range(input.number_axial-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
    #             depfile.write(' '+str(j).ljust(3))
    #             for k in range(columncount):
    #                 depfile.write('{:.3f}'.format(input.boc_exposure).rjust(8))
    #             depfile.write('\n')
    #         depfile.write('\n')
    #     depfile.write(' END STEP\n')
    
    filename = solution.name + '.inp'
    # create input file based on if an input template is given
    if input.calculation_type not in ['multi_cycle']:
        if not input.input_template['apply']: 
            LPs = without_template(solution, input, cwd, filename)
        elif input.input_template['apply']:
            with_template(solution, input, cwd, filename)
    else:
        LPs, initburnup = multi_cycle_input_generation(solution, input, cwd, filename)

    ### surrogate depletion here 
    LPs = "".join(LPs).strip().split()
    ## initializing 
    # xsdict = pickle.load(open(midas_data.__path_xs_pickle__,'rb'))
    # xsdict = joblib.load(midas_data.__path_xs_pickle__,mmap_mode='r')
    fabulist = np.array([0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 
                         17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 
                         60.0, 65.0, 70.0, 75.0, 80.0])
    faaxial = np.array([15.24, 10.16, 5.08, 30.48, 30.48, 30.48, 30.48, 30.48,
                        30.48, 30.48, 30.48, 30.48, 30.48, 5.08, 10.16, 15.24])
    total_height = np.sum(faaxial)
    corebulist = [0.1, 0.4, 0.5, 1.0, 1.0] + [1.0]*28 # modify later to get from input parser
    idx11 = [0]
    idx22 = [1,2,3,4,5,6,7,9,18,27,36,45,54,63]
    if input.calculation_type not in ['multi_cycle']:
        initbumap = np.zeros((1,81,16))
    else:
        initbumap = initburnup
    if input.num_assemblies==157:
        Fd_all, Fq_all, maxboron, cycle_length = get_result_157(LPs, corebulist, xsdict, fabulist, faaxial, total_height, idx11, idx22, initbumap)
    elif input.num_assemblies==193:
        Fd_all, Fq_all, maxboron, cycle_length = get_result(LPs, corebulist, xsdict, fabulist, faaxial, total_height, idx11, idx22, initbumap)

        
    solution.parameters["cycle_length"]['value'] = cycle_length
    solution.parameters["pinpowerpeaking"]['value'] = Fq_all
    solution.parameters["fdeltah"]['value'] = Fd_all
    solution.parameters["max_boron"]['value'] = maxboron
    
    logger.debug(f"Returning to original working directory: {cwd}")
    os.chdir(cwd)
    ## clean xsdict 
    # del xsdict 
    # gc.collect()

    return solution

def extract_data(solution, input, xsdict):
    """
    Optimized data extraction for retraining models using vectorization.
    """
    # 1. Configuration and Path Setup
    retrain_path = midas_data.__path_to_store_retrain_data__
    os.makedirs(retrain_path, exist_ok=True)
    
    cwd = Path(os.getcwd())
    fabulist = np.array([0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 
                         17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 
                         60.0, 65.0, 70.0, 75.0, 80.0])

    # Pre-allocate master lists for this generation
    list_branch_data = [[] for _ in range(7)]
    list_outputRPF = []
    list_outputcore = []
    list_outputdata = [] # for storing true data to perform validation 
    if input.num_assemblies==193:
        flag17 = True
    else:
        flag17 = False

    # 2. Process each individual directory
    for dir_name in solution:
        folder_path = cwd.joinpath(input.results_dir_name / Path(dir_name))
        if input.calculation_type in ['multi_cycle']:
            output_path = os.path.join(folder_path, f"{dir_name}.parcs_cyc-02")## only run n+1 cycle for now
        else:
            output_path = os.path.join(folder_path, f"{dir_name}.parcs_dep")
        input_path = os.path.join(folder_path, f"{dir_name}.parcs_out")
        dplfile_path = os.path.join(folder_path, f"{dir_name}.parcs_dpl")

        if not os.path.exists(output_path):
            continue

        # Load data from PARCS files
        # Converting to numpy arrays immediately is key for speed
        if flag17:
            bu3d = np.array(get_burnup_17(output_path))     # Shape: (Steps, Active_FA, Nodes)
            rpf3d = np.array(get_rpf_17(output_path))       # Shape: (Steps, Active_FA, Nodes)
        else:
            bu3d = np.array(get_burnup_15(output_path))     # Shape: (Steps, Active_FA, Nodes)
            rpf3d = np.array(get_rpf_15(output_path))       # Shape: (Steps, Active_FA, Nodes)
        ## check if data is empty?
        if np.all(bu3d==0): ## if all zeros then skip extracting 
            continue
        coreparam = get_boron_cycle(output_path)     # List of (boron, cycle) for model 2 
        core_data = get_core_data(dplfile_path) ## array of [cyclelength, CBC,Fq,Fd] for retrain purpose
        if input.calculation_type in ['multi_cycle']:
            LPs=getLPs_cyc(output_path)
        else:
            LPs = getLP(input_path)              
        

        n_steps, n_fa_active, nz = bu3d.shape
        n_locations = len(LPs)

        # Pre-allocate temporary arrays for this specific PARCS run
        # Using (Steps, 81, Nodes) to maintain spatial structure
        temp_gen_rpf = np.zeros((n_steps, n_locations, nz))
        temp_gen_xs = [np.zeros((n_steps, n_locations, nz)) for _ in range(7)]

        active_fa_idx = 0
        for loc, fa in enumerate(LPs):
            # Skip reflectors or empty slots
            if fa == '10' or fa == '00':
                continue

            # Vectorization starts here: Iterate nodes, but vectorize burnup steps
            for node in range(nz):
                # 1. Pull XS library data for all 25 burnup points in the library once per node
                # This removes the "float(xs[0])" list comprehension from the inner-most loop
                raw_difc = np.array(getXSlist(xsdict, fa, 'difc', str(node)), dtype=float)
                raw_nufis = np.array(getXSlist(xsdict, fa, 'nufission', str(node)), dtype=float)
                raw_abs = np.array(getXSlist(xsdict, fa, 'absorption', str(node)), dtype=float)
                raw_rem = np.array(getXSlist(xsdict, fa, 'removal', str(node)), dtype=float)

                # 2. Get burnup values for ALL steps at this specific node
                bu_values = bu3d[:, active_fa_idx, node]

                # 3. Vectorized Interpolation: np.interp handles the whole array at once
                temp_gen_xs[0][:, loc, node] = np.interp(bu_values, fabulist, raw_difc[:, 0])
                temp_gen_xs[1][:, loc, node] = np.interp(bu_values, fabulist, raw_difc[:, 1])
                temp_gen_xs[2][:, loc, node] = np.interp(bu_values, fabulist, raw_nufis[:, 0])
                temp_gen_xs[3][:, loc, node] = np.interp(bu_values, fabulist, raw_nufis[:, 1])
                temp_gen_xs[4][:, loc, node] = np.interp(bu_values, fabulist, raw_abs[:, 0])
                temp_gen_xs[5][:, loc, node] = np.interp(bu_values, fabulist, raw_abs[:, 1])
                temp_gen_xs[6][:, loc, node] = np.interp(bu_values, fabulist, raw_rem[:, 0])

                # Store RPF values across all steps
                temp_gen_rpf[:, loc, node] = rpf3d[:, active_fa_idx, node]

            active_fa_idx += 1
        # Collect successful data blocks
        list_outputcore.extend(coreparam)
        list_outputdata.append(core_data)
        list_outputRPF.append(temp_gen_rpf)
        for i in range(7):
            list_branch_data[i].append(temp_gen_xs[i])
    # 3. Final Aggregation and Disk Saving
    if not list_outputRPF:
        raise ValueError("Extraction failed: No valid PARCS data found. Check your .parcs_dep files.")
    # Combine new data with existing .npy files on disk
    # We load them, stack them, and save them.
    # Note: load_or_empty is your helper that returns an empty array if file missing.
    try:
        # Process RPF and Core parameters
        outputRPF = np.vstack((load_or_empty(os.path.join(retrain_path, 'outputRFP.npy'), (0, 81, 16)), 
                               np.vstack(list_outputRPF)))
        outputcore = np.vstack((load_or_empty(os.path.join(retrain_path, 'outputcore.npy'), (0, 2)), 
                                np.array(list_outputcore)))
        outputdata = np.vstack((load_or_empty(os.path.join(retrain_path, 'outputdata.npy'), (0, 4)), 
                                np.array(list_outputdata)))


        np.save(os.path.join(retrain_path, 'outputRFP.npy'), outputRPF)
        np.save(os.path.join(retrain_path, 'outputcore.npy'), outputcore)
        np.save(os.path.join(retrain_path, 'outputdata.npy'), outputdata)

        # Process the 7 XS branches
        for i in range(7):
            filename = f'branchXS_{i+1}.npy'
            current_stacked = np.vstack(list_branch_data[i])
            existing_data = load_or_empty(os.path.join(retrain_path, filename), (0, 81, 16))
            
            final_branch = np.vstack((existing_data, current_stacked))
            np.save(os.path.join(retrain_path, filename), final_branch)
            
            # Explicit cleanup to keep RAM low during large stacks
            del existing_data, final_branch, current_stacked

        print(f"Data extraction complete. Current RPF shape: {outputRPF.shape}")

    except Exception as e:
        print(f"Error during final data save: {e}")
        raise

def load_or_empty(file_path, shape=(0,)):
    """Loads a numpy file if it exists, otherwise returns an empty array."""
    if os.path.exists(file_path):
        return np.load(file_path)
    else:
        return np.empty(shape)

def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def multi_cycle_input_generation(solution, input, cwd, filename):
    '''
    Generate PARCS343 input for multi-cycle run 
    '''
    list_unique_xs = np.concatenate([value if isinstance(value,list) else np.concatenate(list(value.values()))\
                                    for value in input.xs_list.values()])

    ## Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)
    fuel_locations = [loc for loc in input.core_dict.keys() if 2 < len(loc) <  5]
    soln_fuel_locations = {}
    if input.calculation_type in ['eq_cycle']:
        soln_FAs = optools.Solution.SS_decoder(solution.chromosome)
        for i in range(len(solution.chromosome)):
            soln_fuel_locations[fuel_locations[i]] = soln_FAs[i]
    else:
        for i in range(len(solution.chromosome)):
            soln_fuel_locations[fuel_locations[i]] = solution.chromosome[i]
    
    soln_core_dict = deepcopy(input.core_dict)
    soln_core_dict_tag = deepcopy(input.core_dict)
    for loc, label in soln_fuel_locations.items():
        tag = None
        for fueltype in input.tag_list['fuel']:
            if fueltype[1] == label:
                tag = fueltype[0]
        if not tag:
            raise ValueError(f"FA label '{label}' not found in fuel types ({input.tag_list['fuel']}).")
        soln_core_dict[loc]['Value'] = label ## use label rather than tag here
        soln_core_dict_tag[loc]['Value'] = [label,tag] ## use tag for assigning burnup value
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

    burnuplattice = deepcopy(input.core_lattice) # core lattice filled with chromosome
    for loc, vals in soln_core_dict_tag.items():
        sym_locs = [loc] + vals['Symmetric_Assemblies']
        for j in range(len(burnuplattice)):
            for i in range(len(burnuplattice[j])):
                if burnuplattice[j][i] in sym_locs:
                    if burnuplattice[j][i][0] == "R" and len(soln_core_lattice[j][i]) >= 5: #reflector
                        burnuplattice[j][i] = "10" #!TODO: add support more multiple radial refls.
                    elif vals['Value'] is None:
                        burnuplattice[j][i] = '00'
                    elif not is_integer(vals['Value'][0]):
                        burnuplattice[j][i] = vals['Value'][1]
                    else:
                        burnuplattice[j][i] = '00'
    ## replace '00' with just space 
    for i in range(len(soln_core_lattice)):
        for j in range(len(soln_core_lattice[i])):
            if soln_core_lattice[i][j] == '00':
                soln_core_lattice[i][j] = ' '
            elif soln_core_lattice[i][j] == '10':
                soln_core_lattice[i][j] = '0'
    ## get burnup map 
    boc_dep = (Path(cwd)/str(input.depletion_file)).resolve()

    if input.num_assemblies==193:
        bu3dmap = np.array(get_burnup_17(boc_dep)) ## (n_bu,n_fa,n_axial)
    elif input.num_assemblies==157:
        bu3dmap = np.array(get_burnup_15(boc_dep))
    else:
        raise ValueError(f"{input.num_assemblies} core is not supported in current version.")
    initbu3dmap = bu3dmap[input.depletion_restart_point-1]
    initburnup = np.zeros((1,9,9,16))
    for row in range(burnuplattice.shape[0]):
        for col in range(burnuplattice.shape[1]):
            faidx = int(burnuplattice[row][col])
            if faidx !=0:
                initburnup[:,row,col,:]=initbu3dmap[faidx-1,:]
    ## reshape 
    initburnup = initburnup.reshape((1,81,16))
    ## convert LPs
    locationmap = input.core_lattice.T ## transpose 
    depleted_map = np.array(input.depleted_map).reshape(soln_core_lattice.shape)
    LPs = deepcopy(depleted_map)
    for x in range(soln_core_lattice.shape[0]):
        for y in range(soln_core_lattice.shape[1]):
            if is_integer(soln_core_lattice[x,y]):
                if int(soln_core_lattice[x,y])!=0:
                    LPs[x,y] = soln_core_lattice[x,y].lstrip('-') ## take the FA id
            else:
                idx =  np.where(locationmap == soln_core_lattice[x,y].replace('-', ''))
                if len(idx[0]) >0:
                    LPs[x,y] = depleted_map[idx[0][0], idx[1][0]]
    LP_text = ""
    for x in range(LPs.shape[0]):
        for y in range(LPs.shape[1]):
            LP_text +=  LPs[x,y]+" "
    ## write input 
    ## CaseID Block ##
    with open(filename,"w") as ofile:
        ofile.write("!******************************************************************************\n")
        ofile.write('CASEID {}  \n'.format(solution.name))
        ofile.write("!******************************************************************************\n\n")

    ## CNTL Block ##
    with open(filename,"a") as ofile:
        ofile.write("CNTL\n")
        ofile.write("      RUN_OPTS   F T F F\n")
        if input.th_fdbk['apply']:
            if input.th_fdbk['loc'] is None:
                ofile.write("      TH_FDBK    T\n")
                ofile.write("      INT_TH     T -1\n")
            else: 
                ofile.write("      TH_FDBK    T\n")
                ofile.write(f"      INT_TH     T 1 '{input.th_fdbk['loc']}'\n")
        else:
            ofile.write("      TH_FDBK    F\n")
        ofile.write("      CORE_POWER 100.0\n")
        ofile.write(f"      CORE_TYPE  {input.core_type}\n")
        if input.core_type == "PWR":
            ofile.write("      PPM        1000 1.0 1800.0 10.0\n") #!TODO: this should be a parameterized boron guess value.
        ofile.write("      DEPLETION  T  1.0E-5 T\n")
        if input.calculation_type in ['eq_cycle']:
            ofile.write("      MULT_CYC   T  F\n") #v3.4.2 specific line to enable the MCYCLE block
        ofile.write("      TREE_XS    T  {}  T  T  F  F  T  F  F  F  T  F  T  T  T  F  F \n".format(int(len(list_unique_xs))))
        if input.core_type == "PWR":
            ofile.write("      BANK_POS   100 100 100 100 100 100\n")
        ofile.write("      XE_SM      1 1 1 1\n")
        if input.core_type == "PWR":
            ofile.write("      SEARCH     PPM 1.0 1800 10.0\n")
        elif input.core_type == "BWR":
            ofile.write("      SEARCH     KEFF 0.997\n")
        ofile.write("      MULT_CYC   T\n")
        ofile.write("      XS_EXTRAP  1.0 0.3\n")
        if input.pin_power_recon:
            ofile.write("      PIN_POWER  T\n")
        else:
            ofile.write("      PIN_POWER  F\n")
        ofile.write("      PRINT_OPT  T T T T T F T T T T  T  T  T  T  F  T  T\n")
        # ofile.write("      PLOT_OPTS 0 0 0 0 0 2\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        
    ## PARAM Block ##
    with open(filename,"a") as ofile:
        ofile.write("PARAM\n")
        ofile.write("      LSOLVER     1 1 20\n")
        if input.th_fdbk: #!TODO: temporary solution. This should be replaced with an actual parameter for the kernal.
            # ofile.write("      NODAL_KERN  HYBRID\n")
            ofile.write("      NODAL_KERN  NEMMG\n")
        else:
            ofile.write("      NODAL_KERN  NEMMG\n")
        ofile.write("      CMFD        2\n")
        ofile.write("      DECUSP      2\n")
        ofile.write("      INIT_GUESS  0\n")
        ofile.write("      CONV_SS     1.e-6 5.e-5 1.e-3\n")
        ofile.write("      EPS_ERF     0.010\n")
        ofile.write("      EPS_ANM     0.000001\n")
        ofile.write("      NLUPD_SS    5 5 1\n")
        #!if input.th_fdbk:
        #!    ofile.write("      N_ITERS_SS  4 1000\n")
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
        ## insert the core depleted core map
        for x in range(soln_core_lattice.shape[0]):
            ofile.write("      ")
            for y in range(soln_core_lattice.shape[1]):
                ofile.write(depleted_map[x,y])
                ofile.write("  ")
            ofile.write("\n")
        ofile.write("\n")
    
        assembly_width = 21.50 #!TODO: change this to an input with default.
        if input.map_size == 'quarter':
            ofile.write(f"      GRID_X      1*{assembly_width/2} {dim_size[0]-1:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  1*1 {dim_size[0]-1:.0f}*1\n")
            ofile.write(f"      GRID_Y      1*{assembly_width/2} {dim_size[0]-1:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  1*1 {dim_size[0]-1:.0f}*1\n")
        else: #assume full geometry if not quarter-core
            ofile.write(f"      GRID_X      {dim_size[0]:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  {dim_size[0]:.0f}*1\n")
            ofile.write(f"      GRID_Y      {dim_size[1]:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  {dim_size[1]:.0f}*1\n")
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
                # ofile.write("      ASSY_TYPE   {}   1*1  1*{} {}*{}  1*{}  1*{} FUEL\n".format(fuel['type'],xsnum_blanket,\
                #                                                                     input.number_axial-4,xsnum_fuel,\
                #                                                                     xsnum_blanket,xsnum_radtop))
                ofile.write("      ASSY_TYPE   {}   1*1  1*{} 1*{} {}*{}  1*{} 1*{}  1*{} FUEL\n".format(fuel['type'],xsnum_blanket,
                                                                                                         xsnum_blanket,
                                                                                                         input.number_axial-6,
                                                                                                         xsnum_fuel,
                                                                                                         xsnum_blanket,
                                                                                                         xsnum_blanket,
                                                                                                         xsnum_radtop))
            else:
                if fuel['serial']!='False':
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
                    if str(input.pincal_loc[x,y]) not in ['nan']:
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
        if input.th_fdbk['apply'] and not input.th_fdbk['loc']:
            ofile.write("      FLU_TYP       0\n")
            ofile.write("      N_PINGT    264 25\n")
            ofile.write("      PIN_DIM      4.1 4.75 0.58 6.13\n")
            ofile.write("      FLOW_COND    {}  {}\n".format(np.round(input.inlet_temp-273.15,2),\
                                                            np.round(input.flow/input.num_assemblies,4)))
            ofile.write("      STATE_CORE   {}  1301.86  1.5789E7\n".format(np.round(input.flow)))
            ofile.write("      HGAP     11356.0\n") #!TODO:check this value, should it be parameterized?
            ofile.write("      N_RING   6\n")
            ofile.write(f"      THMESH_X       {dim_size[0]}*1\n")
            ofile.write(f"      THMESH_Y       {dim_size[1]}*1\n")
            ofile.write(f"      THMESH_Z       {str([x+1 for x in range(input.number_axial)]).strip('[]').replace(',','')}\n")
        elif not input.th_fdbk['apply']:
            ofile.write("      UNIF_TH   0.740    626.85     {}\n".format(np.round(input.inlet_temp-273.15,2))) #!TODO: how to deal with av. fuel temp?
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        ## DEPL Block ##
    with open(filename,"a") as ofile:
        ofile.write("DEPL\n")
        ofile.write(f"      TIME_STP  {str(input.depl_steps).strip('[]')}\n")
        ofile.write(f"      INP_HST   '{str(boc_dep)}' 1 {str(input.depletion_restart_point)}\n") 
        # ofile.write("      OUT_OPT   T  T  T  T  F\n")
        # Write reflector cross sections
        ofile.write("      PMAXS_F   1 '{}{}' 1\n".format(input.xs_lib / Path(input.xs_list['reflectors']['bottom'][0]),\
                                                        input.xs_extension))
        # ofile.write("      PMAXS_F   1 '{}{}' 1\n".format(input.xs_lib / Path(input.xs_list['reflectors']['bot'][0]),\
        #                                                 input.xs_extension))
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
            # print(input.xs_list['fuel'][i])
            if input.xs_list['fuel'][i]!='False':
                ofile.write("      PMAXS_F   {} '{}{}' {}\n".format(fxs_index,\
                                                                input.xs_lib / Path(input.xs_list['fuel'][i]),\
                                                                input.xs_extension,fxs_index))
    ## MCYCLE Block ##
    with open(filename,"a") as ofile:
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        
        ofile.write("MCYCLE\n")
        ofile.write("    CYCLE_DEF   1\n")
        ofile.write(f"      DEPL_STEP 0\n")
        ofile.write(f"      POWER_LEV 2*100.0\n")
        ofile.write(f"      BANK_SEQ  2*1\n\n")
        ofile.write("    CYCLE_DEF   2\n")
        ofile.write(f"      DEPL_STEP {str(input.depl_steps).strip('[]')}\n")
        ofile.write(f"      POWER_LEV {len(input.depl_steps)+1}*100.0\n")
        ofile.write(f"      BANK_SEQ  {len(input.depl_steps)+1}*1\n\n")
        
        ofile.write("    LOCATION   0\n")
        for x in range(input.full_core_locs.shape[0]):
            ofile.write("      ")
            for y in range(input.full_core_locs.shape[1]):
                # val = input.full_core_locs[x,y]
                val = input.full_core_locs[y,x]
                try:
                    if not np.isnan(val):
                        # ofile.write(str(input.full_core_locs[x,y]))
                        ofile.write(str(input.full_core_locs[y,x]))
                        ofile.write("  ")
                except TypeError:
                    # ofile.write(str(input.full_core_locs[x,y]))
                    ofile.write(str(input.full_core_locs[y,x]))
                    ofile.write("  ")
            ofile.write("\n")
        ofile.write("\n")
        
        ofile.write("    SHUF_MAP   1   1\n")
        for x in range(soln_core_lattice.shape[0]):
            ofile.write("      ")
            for y in range(soln_core_lattice.shape[1]):
                val = soln_core_lattice[x,y]
                ofile.write(str(soln_core_lattice[x,y]))
                ofile.write("  ")
            ofile.write('\n')
        ofile.write('\n')
        
        ofile.write("    CYCLE_IND    1  0  1\n")
        ofile.write(f"    CYCLE_IND   2  1  2\n")
        ofile.write(f"    CONV_EC    0.1  2\n")
    
    ## Termination Character ##
    with open(filename,"a") as ofile:
        ofile.write(".")
    return LP_text, initburnup

def without_template(solution, input, cwd, filename):    
    ## Prepare values for file writing
    list_unique_xs = np.concatenate([value if isinstance(value,list) else np.concatenate(list(value.values()))\
                                    for value in input.xs_list.values()])

    ## Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)
    fuel_locations = [loc for loc in input.core_dict.keys() if 2 < len(loc) <  5]
    soln_fuel_locations = {}
    if input.calculation_type in ['eq_cycle']:
        soln_FAs = optools.Solution.SS_decoder(solution.chromosome)
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
            raise ValueError(f"FA label '{label}' not found in fuel types ({input.tag_list['fuel']}).")
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
    ## CaseID Block ##
    with open(filename,"w") as ofile:
        ofile.write("!******************************************************************************\n")
        ofile.write('CASEID {}  \n'.format(solution.name))
        ofile.write("!******************************************************************************\n\n")

    ## CNTL Block ##
    with open(filename,"a") as ofile:
        ofile.write("CNTL\n")
        ofile.write("      RUN_OPTS   F T F F\n")
        if input.th_fdbk['apply']:
            if input.th_fdbk['loc'] is None:
                ofile.write("      TH_FDBK    T\n")
                ofile.write("      INT_TH     T -1\n")
            else: 
                ofile.write("      TH_FDBK    T\n")
                ofile.write(f"      INT_TH     T 1 '{input.th_fdbk['loc']}'\n")
        else:
            ofile.write("      TH_FDBK    F\n")
        ofile.write("      CORE_POWER 100.0\n")
        ofile.write(f"      CORE_TYPE  {input.core_type}\n")
        if input.core_type == "PWR":
            ofile.write("      PPM        1000 1.0 1800.0 10.0\n") #!TODO: this should be a parameterized boron guess value.
        ofile.write("      DEPLETION  T  1.0E-5 T\n")
        if input.calculation_type in ['eq_cycle']:
            ofile.write("      MULT_CYC   T  F\n") #v3.4.2 specific line to enable the MCYCLE block
        ofile.write("      TREE_XS    T  {}  T  T  F  F  T  F  F  F  T  F  T  T  T  F  F \n".format(int(len(list_unique_xs))))
        if input.core_type == "PWR":
            ofile.write("      BANK_POS   100 100 100 100 100 100\n")
        ofile.write("      XE_SM      1 1 1 1\n")
        if input.core_type == "PWR":
            ofile.write("      SEARCH     PPM\n")
        elif input.core_type == "BWR":
            ofile.write("      SEARCH     KEFF 0.997\n")
        ofile.write("      XS_EXTRAP  1.0 0.3\n")
        if input.pin_power_recon:
            ofile.write("      PIN_POWER  T\n")
        else:
            ofile.write("      PIN_POWER  F\n")
        ofile.write("      PRINT_OPT  T T T T T F T T T T  T  T  T  T  F  T  T")
        # ofile.write("      PLOT_OPTS 0 0 0 0 0 2\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        
    ## PARAM Block ##
    with open(filename,"a") as ofile:
        ofile.write("PARAM\n")
        ofile.write("      LSOLVER     1 1 20\n")
        if input.th_fdbk: #!TODO: temporary solution. This should be replaced with an actual parameter for the kernal.
            ofile.write("      NODAL_KERN  HYBRID\n")
        else:
            ofile.write("      NODAL_KERN  NEMMG\n")
        ofile.write("      CMFD        2\n")
        ofile.write("      DECUSP      2\n")
        ofile.write("      INIT_GUESS  0\n")
        ofile.write("      CONV_SS     1.e-6 5.e-5 1.e-3\n")
        #!ofile.write("      EPS_ERF     0.010\n")
        ofile.write("      EPS_ANM     0.000001\n")
        ofile.write("      NLUPD_SS    3 5 1\n")
        #!if input.th_fdbk:
        #!    ofile.write("      N_ITERS_SS  4 1000\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
    
    ## GEOM Block Inputs ##
    LP_text = ""
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
                LP_text +=  soln_core_lattice[x,y]+" "
                ofile.write("  ")
            ofile.write("\n")
            LP_text += "\n"
        ofile.write("\n")
    
        assembly_width = 21.50 #!TODO: change this to an input with default.
        if input.map_size == 'quarter':
            ofile.write(f"      GRID_X      1*{assembly_width/2} {dim_size[0]-1:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  1*1 {dim_size[0]-1:.0f}*1\n")
            ofile.write(f"      GRID_Y      1*{assembly_width/2} {dim_size[0]-1:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  1*1 {dim_size[0]-1:.0f}*1\n")
        else: #assume full geometry if not quarter-core
            ofile.write(f"      GRID_X      {dim_size[0]:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_X  {dim_size[0]:.0f}*1\n")
            ofile.write(f"      GRID_Y      {dim_size[1]:.0f}*{assembly_width}\n")
            ofile.write(f"      NEUTMESH_Y  {dim_size[1]:.0f}*1\n")
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
                # ofile.write("      ASSY_TYPE   {}   1*1  1*{} {}*{}  1*{}  1*{} FUEL\n".format(fuel['type'],xsnum_blanket,\
                #                                                                     input.number_axial-4,xsnum_fuel,\
                #                                                                     xsnum_blanket,xsnum_radtop))
                ofile.write("      ASSY_TYPE   {}   1*1  1*{} 1*{} {}*{}  1*{} 1*{}  1*{} FUEL\n".format(fuel['type'],xsnum_blanket,
                                                                                                         xsnum_blanket,
                                                                                                         input.number_axial-6,
                                                                                                         xsnum_fuel,
                                                                                                         xsnum_blanket,
                                                                                                         xsnum_blanket,
                                                                                                         xsnum_radtop))
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
                    if str(input.pincal_loc[x,y]) not in ['nan']:
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
        if input.th_fdbk['apply'] and not input.th_fdbk['loc']:
            ofile.write("      FLU_TYP       0\n")
            ofile.write("      N_PINGT    264 25\n")
            ofile.write("      PIN_DIM      4.1 4.75 0.58 6.13\n")
            ofile.write("      FLOW_COND    {}  {}\n".format(np.round(input.inlet_temp-273.15,2),\
                                                            np.round(input.flow/input.num_assemblies,4)))
            ofile.write("      STATE_CORE   {}  1301.86  1.5789E7\n".format(np.round(input.flow)))
            ofile.write("      HGAP     11356.0\n") #!TODO:check this value, should it be parameterized?
            ofile.write("      N_RING   6\n")
            ofile.write(f"      THMESH_X       {dim_size[0]}*1\n")
            ofile.write(f"      THMESH_Y       {dim_size[1]}*1\n")
            ofile.write(f"      THMESH_Z       {str([x+1 for x in range(input.number_axial)]).strip('[]').replace(',','')}\n")
        elif not input.th_fdbk['apply']:
            ofile.write("      UNIF_TH   0.740    626.85     {}\n".format(np.round(input.inlet_temp-273.15,2))) #!TODO: how to deal with av. fuel temp?
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
    ## DEPL Block ##
    with open(filename,"a") as ofile:
        ofile.write("DEPL\n")
        if input.calculation_type == 'single_cycle':
            ofile.write(f"      TIME_STP  {str(input.depl_steps).strip('[]')}\n")
        # ofile.write("      INP_HST   './boc_exp.dep' -2 1\n") TODO add functionality in through input yaml files
        ofile.write("      OUT_OPT   T  T  T  T  F\n")
        # Write reflector cross sections
        ofile.write("      PMAXS_F   1 '{}{}' 1\n".format(input.xs_lib / Path(input.xs_list['reflectors']['bottom'][0]),\
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
    
    return LP_text

def with_template(solution, input, cwd, filename): 

    ## Prepare values for file writing
    ## Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)
    fuel_locations = [loc for loc in input.core_dict.keys() if 2 < len(loc) <  5]
    soln_fuel_locations = {}
    if input.calculation_type in ['eq_cycle']:
        soln_FAs = optools.Solution.SS_decoder(solution.chromosome)
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
            raise ValueError(f"FA label '{label}' not found in fuel types ({input.tag_list['fuel']}).")
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
    ## copy input file from template
    inp_template = str(cwd.joinpath(cwd / input.input_template['loc']))
    shutil.copy(inp_template, filename)
    soln_full_core_lattice = prepare_shuffling_map(input, solution.chromosome)

    with open(filename, "r") as file:
        lines = file.readlines()  

    with open(filename, "w") as ofile:
        for line in lines:
            ## change CaseID ##
            if "caseid" in line.lower():
                ofile.write('CASEID {}  \n'.format(solution.name))  
            ## add LP ##
            elif "rad_conf" in line.lower() and "!" not in line.lower():
                ofile.write("      RAD_CONF\n\n")
                for x in range(soln_core_lattice.shape[0]):
                    ofile.write("      ")
                    for y in range(soln_core_lattice.shape[1]):
                        ofile.write(soln_core_lattice[x,y])
                        ofile.write("  ")
                    ofile.write("\n")
                ofile.write("\n")

            elif 'int_th' in line.lower():
                if input.th_fdbk['apply']:
                    if input.th_fdbk['loc'] is None:
                        ofile.write("      TH_FDBK    T\n")
                        ofile.write("      INT_TH     T -1\n")
                    else: 
                        ofile.write(f"      INT_TH     T 1 '{input.th_fdbk['loc']}'\n")
                else:
                    ofile.write("      TH_FDBK    F\n")
            # add shuffling map
            elif ("location" in line.lower() and "!" not in line.lower()) and input.calculation_type in ['eq_cycle']:
                ofile.write("      LOCATION  0\n")
                for x in range(input.full_core_locs.shape[0]):
                    ofile.write("      ")
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
            
            elif ("shuf_map" in line.lower() and "!" not in line.lower()) and input.calculation_type in ['eq_cycle']:
                ofile.write("      SHUF_MAP 1 1\n")
                for x in range(soln_full_core_lattice.shape[0]):
                    ofile.write("      ")
                    for y in range(soln_full_core_lattice.shape[1]):
                        val = soln_full_core_lattice[x,y]
                        ofile.write(str(soln_full_core_lattice[x,y]))
                        ofile.write("  ")
                    ofile.write('\n')
                ofile.write('\n')

            else:
                ofile.write(line)   

def get_results(parameters, filename, job_failed=False): #!TODO: implement pin power reconstruction.
    """
    Currently supports cycle length, F_q, F_dh, max boron, keff, critical power ratio,
    linear heat generation rate, average planar linear heat generation rate.
    
    Updated by Nicholas Rollins. 09/27/2024
    Updated by Jake Mikouchi. 03/25/2025
    """
    ## Prepare container for results
    results_dict = {}
    for res in ["cycle_length", "pinpowerpeaking", "fdeltah",'pxyz', 'pxy', "max_boron", 
                        "keff_min", "keff_max", "keff_diff", "cpr", "lhgr", "aplhgr", "chfr"]:
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
        
        results_dict["cycle_length"]["value"] = calc_cycle_length(efpd_list,boron_list,keff_list)
        results_dict["pxy"]["value"] = max(pxy_list)
        results_dict["pxyz"]["value"] = max(pxyz_list)
        results_dict["pinpowerpeaking"]["value"] = max(fq_list)
        results_dict["fdeltah"]["value"] = max(fdh_list)
        results_dict["max_boron"]["value"] = max(boron_list)
        results_dict["chfr"]["value"] = min(chfr)
        
        results_dict["keff_min"]["value"] = min(keff_list)
        results_dict["keff_max"]["value"] = max(keff_list)
        results_dict["keff_diff"]["value"] = max(keff_list) - min(keff_list)
        if "cpr" in parameters.keys(): 
            results_dict["cpr"]["value"] = calc_cpr(filename, parameters)
        if "lhgr" in parameters.keys(): 
            results_dict["lhgr"]["value"] = calc_lhgr(fq_list, parameters)
        if "aplhgr" in parameters.keys(): 
            results_dict["aplhgr"]["value"] = calc_aplhgr(filename, parameters)


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
            results_dict["max_boron"]["value"] = new_max_boron
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["cycle_length"]["value"] = 0.0
        results_dict["pinpowerpeaking"]["value"] = 10.0
        results_dict["fdeltah"]["value"] = 10.0
        results_dict["pxy"]["value"] = 10.0
        results_dict["pxyz"]["value"] = 10.0
        results_dict["max_boron"]["value"] = 10000
        results_dict["chfr"]["value"] = 0.0
        results_dict["keff_min"]["value"] = 0.0
        results_dict["keff_max"]["value"] = 10.0
        results_dict["keff_diff"]["value"] = 10.0
        results_dict["cpr"]["value"] = 0.0
        results_dict["lhgr"]["value"] = 100.0
        results_dict["aplhgr"]["value"] = 100.0
    
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
        #!else: #!TODO: is this practical? It would need to have all the TRACE parameters whitelisted as well.
        #!    if param not in ['cost_fuelcycle','av_fuelenrichment']: #check whitelist
        #!        logger.warning(f"Parameter '{param}' not supported in PARCS343 results parsing.")
    
    return parameters

def calc_cycle_length(efpd,boron,keff):
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
    

## function to get EXPOSURE each step for surrogate 

def getLP(inpfile):
    """
    Extract LPs from input
    
    :param inpfile: inputfile path
    """
    with open(inpfile,'r') as f:
        txt = f.readlines()
    for i,line in enumerate(txt):
        if line.find("GEO_DIM")>=0:
            coresize = int(float(line.split()[1]))
        elif line.find('RAD_CONF')>=0:
            st = i+1
            break
    LPs=''
    for ii in range(st, st+coresize):
        LPs+=txt[ii]
    LPs = "".join(LPs).strip().split()

    return LPs

def getLPs_cyc(ofile):
    """
    Extract LPs from cycle file in PARCS in case of multi-cycle calculation
    
    :param ofile: outputfile path
    """
    with open(ofile,'r') as f:
        txt = f.readlines()
    for i,line in enumerate(txt):
        if line.find("Assembly Type Layout")>=0:
            st=i+2
            coresize = len(txt[st].strip().split())
            break
    LPs=''
    for ii in range(st, st+coresize):
        ## add '00' if not equal to xore size
        temp = txt[ii]
        parts = temp.strip().split()
        linenew = ' '.join(parts + ['00'] * (9 - len(parts)))
        LPs+=linenew+'\n'
    LPs = "".join(LPs).strip().split()

    return LPs

def get_boron_cycle(ofile):
    """
    Get boron and cycle at all burn up step 
    
    :param ofile: depletion file 
    """
    with open(ofile,'r') as f:
        txt = f.readlines()
    bo_cyc = []
    key = 'PT   RE     Keff '
    for i,line in enumerate(txt):
        if line.find(key)>=0:
            temp = line.split()
            bor_idx = temp.index("ppm")
            cyc_idx = temp.index("Days")
            val = txt[i+1].split()
            bo_cyc.append([float(val[bor_idx]), float(val[cyc_idx])])
    return bo_cyc

def get_burnup_17(ofile,FULL_CORE=False):
    '''
    3D burnup from PARCS .dep output files
    Some geometry predefined parameters are required.
    '''
    nfa=56
    z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    nz=len(z_id)
    refl_id=[9,18,27,36,44,45,53,60,61,66,67,68,69,70,71,72,73]
    if FULL_CORE:
        refl_id=[1,2,3,4,5,6,7,8,9,10,11,12,
                    20,21,22,23,24,
                    36,37,38,
                    52,53,54,
                    68,69,70,
                    86,87,
                    103,104,
                    120,121,
                    137,138,
                    154,155,
                    171,172,
                    188,189,190,
                    204,205,206,
                    220,221,222,
                    234,235,236,237,238,
                    246,247,248,249,250,251,252,253,254,255,256,257]
    with open(ofile,'r') as f:
        txt = f.readlines()
    alldeplines = []
    k_sta = 'EXP 3D MAP'
    k_end = ' I_D 2D MAP'
    for i,line in enumerate(txt):
        if line.find(k_sta)>=0:
            i_sta = i+1
        if line.find(k_end)>=0:
            i_end = i-1
            depline = [txt[ii] for ii in range(i_sta, i_end+1)]
            depline = "".join(depline)
            alldeplines.append(depline)
    nbu=len(alldeplines)
    bu_3d=np.zeros((nbu,nfa,nz))
    for bu,step in enumerate(alldeplines):
        txt_dep=step.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        asb_counter=0
        fasb_counter=0
        ifass = 0
        iass = 0
        for i in range(1,len(txt_dep)):
            line_dep = txt_dep[i].split()
            if line_dep[0]=='k':
                asb_counter=copy.deepcopy(iass)
                fasb_counter=copy.deepcopy(ifass)
            elif int(line_dep[0]) not in z_id: 
                pass
            else:
                iz = z_id[-1] - int(line_dep[0])
                ifass = copy.deepcopy(fasb_counter)
                iass = copy.deepcopy(asb_counter)
                for j in range(1,len(line_dep)):
                    iass +=1
                    if iass in refl_id:
                        pass
                    else:
                        ifass +=1
                        val = float(line_dep[j])
                        bu_3d[bu, ifass-1, iz]=val
    
    return bu_3d

def get_rpf_17(ofile,FULL_CORE=False):
    '''
    3D RPF from PARCS .dep output files
    Some geometry predefined parameters are required.
    '''
    nfa=56
    z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    nz=len(z_id)
    refl_id=[9,18,27,36,44,45,53,60,61,66,67,68,69,70,71,72,73]
    if FULL_CORE:
        refl_id=[1,2,3,4,5,6,7,8,9,10,11,12,
                    20,21,22,23,24,
                    36,37,38,
                    52,53,54,
                    68,69,70,
                    86,87,
                    103,104,
                    120,121,
                    137,138,
                    154,155,
                    171,172,
                    188,189,190,
                    204,205,206,
                    220,221,222,
                    234,235,236,237,238,
                    246,247,248,249,250,251,252,253,254,255,256,257]
    with open(ofile,'r') as f:
        txt = f.readlines()
    alldeplines = []
    k_sta = ' RPF 3D MAP'
    k_end = ' EXP 2D MAP'
    for i,line in enumerate(txt):
        if line.find(k_sta)>=0:
            i_sta = i+1
        if line.find(k_end)>=0:
            i_end = i-1
            depline = [txt[ii] for ii in range(i_sta, i_end+1)]
            depline = "".join(depline)
            alldeplines.append(depline)
    nbu=len(alldeplines)
    rpf_3d=np.zeros((nbu,nfa,nz))
    for bu,step in enumerate(alldeplines):
        txt_dep=step.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        asb_counter=0
        fasb_counter=0
        ifass = 0
        iass = 0
        for i in range(1,len(txt_dep)):
            line_dep = txt_dep[i].split()
            if line_dep[0]=='k':
                asb_counter=copy.deepcopy(iass)
                fasb_counter=copy.deepcopy(ifass)
            elif int(line_dep[0]) not in z_id: 
                pass
            else:
                iz = z_id[-1] - int(line_dep[0])
                ifass = copy.deepcopy(fasb_counter)
                iass = copy.deepcopy(asb_counter)
                for j in range(1,len(line_dep)):
                    iass +=1
                    if iass in refl_id:
                        pass
                    else:
                        ifass +=1
                        val = float(line_dep[j])
                        rpf_3d[bu, ifass-1, iz]=val
    return rpf_3d

def get_burnup_15(ofile,FULL_CORE=False):
    '''
    3D burnup from PARCS .dep output files
    Some geometry predefined parameters are required.
    '''
    nfa=47
    z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    nz=len(z_id)
    refl_id=[9, 18, 26, 27, 35, 42, 43, 49, 50, 55, 56, 59, 60, 61, 62, 63, 64]
    if FULL_CORE:
        refl_id=[1,2,3,4,5,6,7,8,9,10,11,12,
                    20,21,22,23,24,
                    36,37,38,
                    52,53,54,
                    68,69,70,
                    86,87,
                    103,104,
                    120,121,
                    137,138,
                    154,155,
                    171,172,
                    188,189,190,
                    204,205,206,
                    220,221,222,
                    234,235,236,237,238,
                    246,247,248,249,250,251,252,253,254,255,256,257]
    with open(ofile,'r') as f:
        txt = f.readlines()
    alldeplines = []
    k_sta = 'EXP 3D MAP'
    k_end = ' I_D 2D MAP'
    for i,line in enumerate(txt):
        if line.find(k_sta)>=0:
            i_sta = i+1
        if line.find(k_end)>=0:
            i_end = i-1
            depline = [txt[ii] for ii in range(i_sta, i_end+1)]
            depline = "".join(depline)
            alldeplines.append(depline)
    nbu=len(alldeplines)
    bu_3d=np.zeros((nbu,nfa,nz))
    for bu,step in enumerate(alldeplines):
        txt_dep=step.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        asb_counter=0
        fasb_counter=0
        ifass = 0
        iass = 0
        for i in range(1,len(txt_dep)):
            line_dep = txt_dep[i].split()
            if line_dep[0]=='k':
                asb_counter=copy.deepcopy(iass)
                fasb_counter=copy.deepcopy(ifass)
            elif int(line_dep[0]) not in z_id: 
                pass
            else:
                iz = z_id[-1] - int(line_dep[0])
                ifass = copy.deepcopy(fasb_counter)
                iass = copy.deepcopy(asb_counter)
                for j in range(1,len(line_dep)):
                    iass +=1
                    if iass in refl_id:
                        pass
                    else:
                        ifass +=1
                        val = float(line_dep[j])
                        bu_3d[bu, ifass-1, iz]=val
    
    return bu_3d

def get_rpf_15(ofile,FULL_CORE=False):
    '''
    3D RPF from PARCS .dep output files
    Some geometry predefined parameters are required.
    '''
    nfa=47
    z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
    nz=len(z_id)
    refl_id=[9, 18, 26, 27, 35, 42, 43, 49, 50, 55, 56, 59, 60, 61, 62, 63, 64]
    if FULL_CORE:
        refl_id=[1,2,3,4,5,6,7,8,9,10,11,12,
                    20,21,22,23,24,
                    36,37,38,
                    52,53,54,
                    68,69,70,
                    86,87,
                    103,104,
                    120,121,
                    137,138,
                    154,155,
                    171,172,
                    188,189,190,
                    204,205,206,
                    220,221,222,
                    234,235,236,237,238,
                    246,247,248,249,250,251,252,253,254,255,256,257]
    with open(ofile,'r') as f:
        txt = f.readlines()
    alldeplines = []
    k_sta = ' RPF 3D MAP'
    k_end = ' EXP 2D MAP'
    for i,line in enumerate(txt):
        if line.find(k_sta)>=0:
            i_sta = i+1
        if line.find(k_end)>=0:
            i_end = i-1
            depline = [txt[ii] for ii in range(i_sta, i_end+1)]
            depline = "".join(depline)
            alldeplines.append(depline)
    nbu=len(alldeplines)
    rpf_3d=np.zeros((nbu,nfa,nz))
    for bu,step in enumerate(alldeplines):
        txt_dep=step.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        asb_counter=0
        fasb_counter=0
        ifass = 0
        iass = 0
        for i in range(1,len(txt_dep)):
            line_dep = txt_dep[i].split()
            if line_dep[0]=='k':
                asb_counter=copy.deepcopy(iass)
                fasb_counter=copy.deepcopy(ifass)
            elif int(line_dep[0]) not in z_id: 
                pass
            else:
                iz = z_id[-1] - int(line_dep[0])
                ifass = copy.deepcopy(fasb_counter)
                iass = copy.deepcopy(asb_counter)
                for j in range(1,len(line_dep)):
                    iass +=1
                    if iass in refl_id:
                        pass
                    else:
                        ifass +=1
                        val = float(line_dep[j])
                        rpf_3d[bu, ifass-1, iz]=val
    return rpf_3d

def getXSlist(xsdict,fatype,xstype,axialnode):
    bukeys=xsdict[fatype].keys()
    xsval = [xsdict[fatype][key][xstype][axialnode] for key in bukeys]
    return xsval

def get_core_data(dplfile):
    '''
    Get cycle length, CBC, Fq,Fd for retrain purpose
    '''
    with open(dplfile) as ofile:
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
    cycle = calc_cycle_length(efpd_list,boron_list,keff_list)
    cbc = max(boron_list)
    fq = max(fq_list)
    fd = max(fdh_list)
    return np.array([cycle,cbc,fq,fd])

