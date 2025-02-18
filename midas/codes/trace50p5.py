## Import Block ##
import os
import gc
import logging
import shutil
import numpy as np
from copy import deepcopy
from pathlib import Path
from time import sleep
import subprocess
from subprocess import STDOUT
from midas.utils.optimizer_tools import Constrain_Input
from midas_data import __trace50p5exe__
from midas.codes import parcs343


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

## Functions ##
def evaluate(solution, input): #!TODO: this doesn't include the initial PARCS calculation needed to get the EOC.dep file.
    """
    Interface used to run multiphysics TRACE-PARCS calculations, including transients. This assumes TRACEv50p5 and PARCSv32m21co.
    #!TODO: In its current configuration, this interface has an exagurated reliance on user-provided template files.
    
    Written by Nicholas Rollins. 2/12/2025
    """
## Run Initial Solution
    if input.init_code == "parcs343":
        solution = parcs343.evaluate(solution, input)
    
    # using the "evaluate" function of another code will automatically create the job directory, so this step is skipped.
    cwd = Path(os.getcwd())
    indv_dir = cwd.joinpath(input.results_dir_name / Path(solution.name))
    logger.debug(f"Continuing to working directory: {indv_dir}")
    os.chdir(indv_dir)

## TRACE SS input file
    shutil.copy(cwd / input.inp_template_ss, indv_dir / Path("TRACE_" + solution.name + "_ss.inp"))
    #!TODO: edit template file to edit inputs to be parameterized such as nhtstr and ncomp.

## PARCS SS input file

    ## Write EOC Exposure Map
    writeDEPfromOUT("eoc_exp.dep", solution.name, input.pincal_loc, input.number_axial)
    
    ## Prepare values for file writing
    list_unique_xs = np.concatenate([value if isinstance(value,list) else np.concatenate(list(value.values()))\
                                    for value in input.xs_list.values()])

    # Fill loading pattern with chromosome (core_dict from Prepare_Problem_Values.prepare_cycle)
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
    
    filename = "TRACE_" + solution.name + '_ss.parcs_inp'
    
    ## CaseID Block ##
    with open(filename,"w") as ofile:
        ofile.write("!******************************************************************************\n")
        ofile.write('CASEID {}               TRACE-PARCS SS\n'.format(solution.name))
        ofile.write("!******************************************************************************\n\n")

    ## CNTL Block ##
    with open(filename,"a") as ofile:
        ofile.write("CNTL\n")
        ofile.write("      ext_th     T     '{}'   TRAC   10\n".format(cwd / input.inp_maptabfile))
        ofile.write("      CORE_POWER {}\n".format(input.ss_powerfraction))
        ofile.write("      CORE_TYPE  PWR\n")
        ofile.write("      DEPLETION  T  1.0E-3 T\n")
        ofile.write("      TREE_XS    T  {}  T  T  F  F  T  F  T  F  T  F  F  T  T  F \n".format(int(len(list_unique_xs))))
        ofile.write("      BANK_POS   100 100 100 100 100 100 100 100\n") #!TODO: parameterize control rods
        ofile.write("      XE_SM      1 1\n")
        ofile.write("      XS_EXTRAP  1.0 0.3\n")
        if input.pin_power_recon:
            ofile.write("      PIN_POWER  T\n")
        else:
            ofile.write("      PIN_POWER  F\n")
        ofile.write("      PRINT_OPT  T T T T T F T T T T  T  T  T  T  F  T  T")
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
        
        if input.pin_power_recon:
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
        #!TODO: this needs to be parameterized. Not sure the best way to do that for bank_conf, though.
        ofile.write("      cr_axinfo  36.2255 3.427055\n")
        ofile.write("      bank_conf\n")
        ofile.write("                        0   0   0   0   0   0   0   0   0               \n")
        ofile.write("                0   0   0   0   0   0   0   0   0   0   0   0   0       \n")
        ofile.write("            0   0   0   0   0   1   0   6   0   1   0   0   0   0   0   \n")
        ofile.write("            0   0   0   0   3   0   5   0   5   0   3   0   0   0   0   \n")
        ofile.write("        0   0   0   0   7   0   8   0   7   0   8   0   7   0   0   0   0\n")
        ofile.write("        0   0   0   3   0   5   0   4   0   4   0   5   0   3   0   0   0\n")
        ofile.write("        0   0   1   0   8   0   6   0   2   0   6   0   8   0   1   0   0\n")
        ofile.write("        0   0   0   5   0   4   0   2   0   2   0   4   0   5   0   0   0\n")
        ofile.write("        0   0   6   0   7   0   2   0   7   0   2   0   7   0   6   0   0\n")
        ofile.write("        0   0   0   5   0   4   0   2   0   2   0   4   0   5   0   0   0\n")
        ofile.write("        0   0   1   0   8   0   6   0   2   0   6   0   8   0   1   0   0\n")
        ofile.write("        0   0   0   3   0   5   0   4   0   4   0   5   0   3   0   0   0\n")
        ofile.write("        0   0   0   0   7   0   8   0   7   0   8   0   7   0   0   0   0\n")
        ofile.write("            0   0   0   0   3   0   5   0   5   0   3   0   0   0   0   \n")
        ofile.write("            0   0   0   0   0   1   0   6   0   1   0   0   0   0   0   \n")
        ofile.write("                0   0   0   0   0   0   0   0   0   0   0   0   0       \n")
        ofile.write("                        0   0   0   0   0   0   0   0   0\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## FDBK Block ##
    with open(filename,"a") as ofile:
        ofile.write("FDBK\n")
        ofile.write("      FA_POWPIT     {} {}\n".format(np.round(input.power/input.num_assemblies,4),assembly_width))
        ofile.write("      GAMMA_FRAC    0.0208    0.0    0.0\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

   ## DEPL Block ##
    with open(filename,"a") as ofile:
        ofile.write("DEPL\n")
        ofile.write("      TIME_STP   0.001\n")
        ofile.write("      INP_HST   './{}' -2 1\n".format("eoc_exp.dep"))
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

    ## Termination Character ##
    with open(filename,"a") as ofile:
        ofile.write(".")

## Run TRACE-PARCS SS Coupled Calc
    try:
        output = subprocess.check_output([__trace50p5exe__,"-p","TRACE_" + solution.name + "_ss.inp"], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
    ## Continue only upon successfull completion
        if 'steady state converged' in str(output): #job completed
            if not input.inp_template_tr: #steady-state calc ONLY.
                logger.debug(f"Job {solution.name} completed successfully in TRACEv50p5.")
                solution.parameters = get_results(solution.parameters, solution.name)
                logger.debug(f"Returning to original working directory: {cwd}")
                os.chdir(cwd)
                gc.collect()
                return solution
            else:
                ss_completed = True #continue on to transient calc.
        else: #job failed
            ss_completed = False
            logger.warning(f"Job {solution.name} has failed!")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
            logger.debug(f"Returning to original working directory: {cwd}")
            os.chdir(cwd)
            gc.collect()
            return solution
    except subprocess.TimeoutExpired: #job timed out
        ss_completed = False
        logger.error(f"Job {solution.name} has timed out!")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        logger.debug(f"Returning to original working directory: {cwd}")
        os.chdir(cwd)
        gc.collect()
        return solution
    except subprocess.CalledProcessError as e: #TRACE returned an abort signal
        ss_completed = False
        logger.error(f"Job {solution.name} has failed with the following exception: {e}")
        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        logger.debug(f"Returning to original working directory: {cwd}")
        os.chdir(cwd)
        gc.collect()
        return solution
    
    if ss_completed:
## rename restart file for TR
        os.rename("./TRACE_" + solution.name + "_ss.tpr","./TRACE_" + solution.name + "_tr.rst")
    
## TRACE TR input file
        shutil.copy(cwd / input.inp_template_tr, indv_dir / Path("TRACE_" + solution.name + "_tr.inp"))
        #!TODO: edit template file to edit inputs to be parameterized such as nhtstr and ncomp.
    
## PARCS TR input file
        filename = "TRACE_" + solution.name + '_tr.parcs_inp'
        
    ## CaseID Block ##
        with open(filename,"w") as ofile:
            ofile.write("!******************************************************************************\n")
            ofile.write('CASEID {}               TRACE-PARCS TR\n'.format(solution.name))
            ofile.write("!******************************************************************************\n\n")

    ## CNTL Block ##
        with open(filename,"a") as ofile:
            ofile.write("CNTL\n")
            ofile.write("      ext_th     T     '{}'   TRAC   10\n".format(cwd / input.inp_maptabfile))
            ofile.write("      CORE_POWER 100\n")
            ofile.write("      CORE_TYPE  PWR\n")
            ofile.write("      DEPLETION  T  1.0E-3 T\n")
            ofile.write("      decay_heat T\n")
            ofile.write("      transient  T\n")
            ofile.write("      restart    T './{}.parcs_rst'  -1\n".format("TRACE_" + solution.name + '_ss'))
            ofile.write("      TREE_XS    T  {}  T  T  F  F  T  F  T  F  T  F  F  T  T  F \n".format(int(len(list_unique_xs))))
            ofile.write("      BANK_POS   100 100 100 100 100 100 100 100\n") #!TODO: parameterize control rods
            ofile.write("      XE_SM      1 1\n")
            ofile.write("      XS_EXTRAP  1.0 0.3\n")
            if input.pin_power_recon:
                ofile.write("      PIN_POWER  T\n")
            else:
                ofile.write("      PIN_POWER  F\n")
            ofile.write("      PRINT_OPT  T T T T T F T T T T  T  T  T  T  F  T  T")
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
            
            if input.pin_power_recon:
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
            #!TODO: this needs to be parameterized. Not sure the best way to do that for bank_conf, though.
            ofile.write("      cr_axinfo  36.2255 3.427055\n")
            ofile.write("      bank_conf\n")
            ofile.write("                        0   0   0   0   0   0   0   0   0               \n")
            ofile.write("                0   0   0   0   0   0   0   0   0   0   0   0   0       \n")
            ofile.write("            0   0   0   0   0   1   0   6   0   1   0   0   0   0   0   \n")
            ofile.write("            0   0   0   0   3   0   5   0   5   0   3   0   0   0   0   \n")
            ofile.write("        0   0   0   0   7   0   8   0   7   0   8   0   7   0   0   0   0\n")
            ofile.write("        0   0   0   3   0   5   0   4   0   4   0   5   0   3   0   0   0\n")
            ofile.write("        0   0   1   0   8   0   6   0   2   0   6   0   8   0   1   0   0\n")
            ofile.write("        0   0   0   5   0   4   0   2   0   2   0   4   0   5   0   0   0\n")
            ofile.write("        0   0   6   0   7   0   2   0   7   0   2   0   7   0   6   0   0\n")
            ofile.write("        0   0   0   5   0   4   0   2   0   2   0   4   0   5   0   0   0\n")
            ofile.write("        0   0   1   0   8   0   6   0   2   0   6   0   8   0   1   0   0\n")
            ofile.write("        0   0   0   3   0   5   0   4   0   4   0   5   0   3   0   0   0\n")
            ofile.write("        0   0   0   0   7   0   8   0   7   0   8   0   7   0   0   0   0\n")
            ofile.write("            0   0   0   0   3   0   5   0   5   0   3   0   0   0   0   \n")
            ofile.write("            0   0   0   0   0   1   0   6   0   1   0   0   0   0   0   \n")
            ofile.write("                0   0   0   0   0   0   0   0   0   0   0   0   0       \n")
            ofile.write("                        0   0   0   0   0   0   0   0   0\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")

    ## FDBK Block ##
        with open(filename,"a") as ofile:
            ofile.write("FDBK\n")
            ofile.write("      FA_POWPIT     {} {}\n".format(np.round(input.power/input.num_assemblies,4),assembly_width))
            ofile.write("      GAMMA_FRAC    0.0208    0.0    0.0\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")

   ## DEPL Block ##
        with open(filename,"a") as ofile:
            ofile.write("DEPL\n")
            ofile.write("      TIME_STP   0.001\n")
            ofile.write("      INP_HST   './{}' -2 1\n".format("eoc_exp.dep"))
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
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")

    ## TRAN Block ##
        with open(filename,"a") as ofile: #!TODO: this needs to correspond with the TRACE file settings and should be parameterized.
            ofile.write("TRAN\n")
            ofile.write("     Time_Step  100.0 0.1 100.0  10.0\n")
            ofile.write("     Scram      T 114.0 0.0 2.2\n")
            ofile.write("     nlupd_tr   5 1 5 10\n")
            ofile.write("     conv_tr    0.0001\n")
            ofile.write("     eps_xsec   0.005\n")

    ## Termination Character ##
        with open(filename,"a") as ofile:
            ofile.write(".")
    
## Run TRACE-PARCS TR Coupled Calc
        try:
            output = subprocess.check_output([__trace50p5exe__,"-p","TRACE_" + solution.name + "_tr.inp"], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
        ## Get Results
            if not '## Fatal Error ##' in str(output): #job completed
                logger.debug(f"Job {solution.name} completed successfully in TRACEv50p5.")
                solution.parameters = get_results(solution.parameters, solution.name)
            
            else: #job failed
                logger.warning(f"Job {solution.name} has failed!")
                solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        
        except subprocess.TimeoutExpired: #job timed out
            logger.error(f"Job {solution.name} has timed out!")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        except subprocess.CalledProcessError as e: #TRACE returned an abort signal
            logger.error(f"Job {solution.name} has failed with the following exception: {e}")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
    
    logger.debug(f"Returning to original working directory: {cwd}")
    os.chdir(cwd)
    gc.collect()
    
    return solution

def writeDEPfromOUT(depFileName, soln_name, fuel_locs_inp, num_axial_nodes):
    """
    Method for automatically parsing and generating the 3D nodel exposure map data and '.dep' file
    needed to initialize the coupled TRACE-PARCS calculation. Assumes EOC conditions are desired.
    
    parameters:
        depFileName: str, name to be used for the '*exp.dep' 3D exposure map file.
        soln_name: str, designation for current generation and individual.
        fuel_locs_inp: This is expected to be the "pincal_loc" attribute from the "input_parser" object.
        num_axial_nodes: int, the number of axial nodes in the PARCS model, including top and bottom reflectors.
    
    Written by Nicholas Rollins. 2/14/2025
    """
    ## Identify PARCS output file containing the 3D exposure map data
    if Path(soln_name + ".parcs_dpl").exists: #multiple PARCS output files present; EOC of last cycle will be used.
        cycles = []
        with open(soln_name + ".parcs_dpl", 'r') as dpl_read:
            for line in dpl_read:
                if " CYCLE   " in line:
                    sline = line.split()
                    cycles.append(int(sline[-1]))
        if cycles:
            lastcycle = max(cycles)
            filename = Path(soln_name + ".parcs_cyc-" + str(lastcycle).rjust(2,'0'))
        else:
            filename = Path(soln_name + ".parcs_dep")
    else:
        filename = Path(soln_name + ".dep") #this is necessary for earlier PARCS versions.
    
    fuel_locs_list = []
    for line in fuel_locs_inp:
        for col in line:
            try:
                fuel_locs_list.append(int(col))
            except:
                pass
    num_assemblies = sum(fuel_locs_list)

    exp_map = np.zeros((num_axial_nodes-2,num_assemblies)) #(axial node, FA index)

    ## Parse exposure map from PARCS output file
    with open(filename, 'r') as fileread:
        loc_offset = -10
        FAi_offset = 0
        line_offset = 0
        invalid = True
        for line in fileread:
            if invalid:
                if " EXP 3D MAP 1.0E+00" in line:
                    invalid = False
                    continue
            else:
                if not line.strip():
                    continue
                elif " I_D 2D MAP" in line:
                    invalid = True
                    loc_offset = -10
                    FAi_offset = 0
                    line_offset = 0
                    continue
                elif " k lb   " in line:
                    FAi_offset += line_offset
                    loc_offset += 10
                else:
                    sline = line.split()
                    axindex = int(sline[0])
                    if 1 < axindex < num_axial_nodes: #omit top and bot refl
                        line_offset = 0
                        for i in range(1,len(sline)):
                            if not fuel_locs_list[i+loc_offset-1]:
                                continue
                            else:
                                exp_map[axindex-2][FAi_offset+line_offset] = float(sline[i])
                                line_offset += 1

    ## Prepare depletion file template
    with open(depFileName,"w") as depfile:
        depfile.write("\n BEGIN STEP\n\n EXP 3D MAP 1.0E+00\n\n")
        columncount = 0
        ioffset = -10
        for i in range(1,num_assemblies+1):
            ## write column headers
            if columncount == 0:
                depfile.write(" k lb ")
                ioffset += 10
            depfile.write(str(i).ljust(8))
            columncount += 1
            ## write rows for every 10 columns
            if columncount == 10:
                depfile.write('\n')
                for j in range(num_axial_nodes-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
                    depfile.write(' '+str(j).ljust(3))
                    for k in range(columncount):
                        depfile.write('{:.3f}'.format((exp_map[j-1][k+ioffset])).rjust(8))
                    depfile.write('\n')
                depfile.write('\n')
                columncount = 0
        ## write rows for leftover columns
        if columncount!= 0:
            depfile.write('\n')
            for j in range(num_axial_nodes-2,0,-1): #iterate in reverse; assume 1 node each top and bottom reflectors.
                depfile.write(' '+str(j).ljust(3))
                for k in range(columncount):
                    depfile.write('{:.3f}'.format(exp_map[j-1][k+ioffset]).rjust(8))
                depfile.write('\n')
            depfile.write('\n')
        depfile.write(' END STEP\n')
    return

def get_results(parameters, filename, job_failed=False): #!TODO: implement pin power reconstruction.
    """
    Currently supports cycle length, F_q, F_dh, and max boron from initial PARCS calc.
    Currently supports maximum inner cladding temperature, fuel centerline temperature, and gap heat transfer in TRACE.
    
    Written by Nicholas Rollins. 02/17/2025
    """
    ## Prepare container for results
    results_dict = {}
    for res in ["maxcladtemp", "maxfueltemp", "maxgapq"]:
        results_dict[res] = {}
        results_dict[res]['value'] = 0.0
        
    if not job_failed:
    ## Read TRACE transient file for parsing
        cold_temp = 562.95 #K; CZP temperature (such as in the radial reflectors) #!TODO: this should be parameterized.
        split_center_row = True #!TODO: this should be parameterized.
        core_shape = [
            [None, None, None, None, 909, 909, 909, 909, 909, 909, 909, 909, 909, None, None, None, None],
            [None, None, 909, 909, 909, 701, 702, 703, 704, 705, 706, 707, 909, 909, 909, None, None],
            [None, 909, 909, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 909, 909, None],
            [None, 909, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 909, None],
            [909, 909, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 909, 909],
            [909, 745, 746, 747, 748, 749, 750, 751, 752, 753, 754, 755, 756, 757, 758, 759, 909],
            [909, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 909],
            [909, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 909],
            [909, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 804, 909],
            [909, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 909],
            [909, 820, 821, 822, 823, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 909],
            [909, 835, 836, 837, 838, 839, 840, 841, 842, 843, 844, 845, 846, 847, 848, 849, 909],
            [909, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 864, 909],
            [909, 909, 865, 866, 867, 868, 869, 870, 871, 872, 873, 874, 875, 876, 877, 909, 909],
            [None, 909, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 909, None],
            [None, 909, 909, 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 909, 909, None],
            [None, None, 909, 909, 909, 902, 903, 904, 905, 906, 907, 908, 909, 909, 909, None, None],
            [None, None, None, None, 909, 909, 909, 909, 909, 909, 909, 909, 909, None, None, None, None]
        ] #STP; #!TODO: this should be parameterized.
        
        ## Parse FA-averaged maps at each timestep and axial node
            #dictionaries organized as {hsnum:[time:axial]}
        CladTemp_dict = {} #K
        CenterlineTemp_dict = {} #K
        GapQ_dict = {} #W/m^2
        with open("TRACE_" + filename + "_tr.out", 'r') as fileread:
            valid = False
            lineskip = 8
            for line in fileread:
                if valid:
                    if lineskip > 0:
                        lineskip -= 1
                        continue
                    else:
                        sline = line.split()
                        if len(sline) <= 2:
                            valid = False
                            continue
                        else:
                            CladTemp_dict[hsnum][-1].append(float(sline[6]))
                            CenterlineTemp_dict[hsnum][-1].append(float(sline[8]))
                            GapQ_dict[hsnum][-1].append(float(sline[3]))
                if " For HS num =" in line:
                    sline = line.split()
                    try:
                        hsnum = int(sline[-1])
                    except ValueError:
                        continue #skip HS GRAVS edits
                    valid = True
                    lineskip = 8

                    try:
                        CladTemp_dict[hsnum].append([])
                    except KeyError:
                        CladTemp_dict[hsnum] = [[]]
                    try:
                        CenterlineTemp_dict[hsnum].append([])
                    except KeyError:
                        CenterlineTemp_dict[hsnum] = [[]]
                    try:
                        GapQ_dict[hsnum].append([])
                    except KeyError:
                        GapQ_dict[hsnum] = [[]]
                    continue
        
        ## Fetch number of timesteps and axial nodes in data
        num_timesteps = len(CladTemp_dict[list(CladTemp_dict.keys())[0]])
        num_axialnodes = len(CladTemp_dict[list(CladTemp_dict.keys())[0]][0])
        
        ## Generate 2D map of values at each timestep and axial node.
            ## This is needed to apply split_center corrections and 
            ##    may be useful in the future if the formatted 2D map is desired for printing.
        for ts in range(num_timesteps):
            for axn in range(num_axialnodes):
                if split_center_row == True: #some calcs split the center row of FAs into different coolant channels.
                    split_center = int(len(core_shape)/2)-1
                clad_coremap = deepcopy(core_shape)
                fuel_coremap = deepcopy(core_shape)
                gapq_coremap = deepcopy(core_shape)
                for row in range(len(core_shape)):
                    for col in range(len(core_shape[row])):
                    ## Format 2D maps from parsed data
                        if core_shape[row][col]:
                            if split_center_row != True:
                                clad_coremap[row][col] = CladTemp_dict[core_shape[row][col]][ts][axn]
                                fuel_coremap[row][col] = CenterlineTemp_dict[core_shape[row][col]][ts][axn]
                                gapq_coremap[row][col] = GapQ_dict[core_shape[row][col]][ts][axn]
                            else:
                                if not row == split_center:
                                    clad_coremap[row][col] = CladTemp_dict[core_shape[row][col]][ts][axn]
                                    fuel_coremap[row][col] = CenterlineTemp_dict[core_shape[row][col]][ts][axn]
                                    gapq_coremap[row][col] = GapQ_dict[core_shape[row][col]][ts][axn]
                                else: #if split_center, recombine center row FA values.
                                    try:
                                        clad_TperQ = [(CladTemp_dict[core_shape[row][col]][ts][axn]-cold_temp)/GapQ_dict[core_shape[row][col]][ts][axn],
                                                 (CladTemp_dict[core_shape[row+1][col]][ts][axn]-cold_temp)/GapQ_dict[core_shape[row+1][col]][ts][axn]]
                                        clad_coremap[row][col] = cold_temp + (GapQ_dict[core_shape[row][col]][ts][axn] + \
                                                                       GapQ_dict[core_shape[row+1][col]][ts][axn])*(clad_TperQ[0] + clad_TperQ[1])/2
                                    except ZeroDivisionError:
                                        clad_coremap[row][col] = CladTemp_dict[core_shape[row][col]][ts][axn]
                                    try:
                                        fuel_TperQ = [(CenterlineTemp_dict[core_shape[row][col]][ts][axn]-cold_temp)/GapQ_dict[core_shape[row][col]][ts][axn],
                                                 (CenterlineTemp_dict[core_shape[row+1][col]][ts][axn]-cold_temp)/GapQ_dict[core_shape[row+1][col]][ts][axn]]
                                        fuel_coremap[row][col] = cold_temp + (GapQ_dict[core_shape[row][col]][ts][axn] + \
                                                                       GapQ_dict[core_shape[row+1][col]][ts][axn])*(fuel_TperQ[0] + fuel_TperQ[1])/2
                                    except ZeroDivisionError:
                                        fuel_coremap[row][col] = CenterlineTemp_dict[core_shape[row][col]][ts][axn]

                                    gapq_coremap[row][col] = GapQ_dict[core_shape[row][col]][ts][axn] + GapQ_dict[core_shape[row+1][col]][ts][axn]
                
                ## Skip second channel for center row of FAs 
                if split_center_row == True:
                    clad_coremap = clad_coremap[:split_center+1] + clad_coremap[split_center+2:]
                    fuel_coremap = fuel_coremap[:split_center+1] + fuel_coremap[split_center+2:]
                    gapq_coremap = gapq_coremap[:split_center+1] + gapq_coremap[split_center+2:]
            
            ## Calculate max value in 2D map and compare to global max value in results.
                cladtemps_list = []
                fueltemps_list = []
                gapqtemps_list = []
                for row in range(len(clad_coremap)):
                    for col in range(len(clad_coremap[row])):
                        try:
                            cladtemps_list.append(float(clad_coremap[row][col]))
                        except TypeError:
                            continue
                        try:
                            fueltemps_list.append(float(fuel_coremap[row][col]))
                        except TypeError:
                            continue
                        try:
                            gapqtemps_list.append(float(gapq_coremap[row][col]))
                        except TypeError:
                            continue
                results_dict["maxcladtemp"]["value"] = max(max(cladtemps_list),results_dict["maxcladtemp"]["value"])
                results_dict["maxfueltemp"]["value"] = max(max(fueltemps_list),results_dict["maxfueltemp"]["value"])
                results_dict["maxgapq"]["value"] = max(max(gapqtemps_list),results_dict["maxgapq"]["value"])
    
    
    else: #job has failed; fill parameters with absurdly negative values.
        results_dict["maxcladtemp"]["value"] = 2000 #K
        results_dict["maxfueltemp"]["value"] = 2000 #K
        results_dict["maxgapq"]["value"] = 1E6 #W/m^2
    
    ## Save user-requested values from results
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
        #!else: #!TODO: is this practical? It would need to have all the PARCS parameters whitelisted as well.
        #!    if param not in ['cost_fuelcycle','av_fuelenrichment']: #check whitelist
        #!        logger.warning(f"Parameter '{param}' not supported in TRACE50p5 results parsing.")
    
    return parameters