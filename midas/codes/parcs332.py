## Import Block ##
import os
import gc
from pathlib import Path


## Functions ##
def get_results(parameters, filename): #!TODO: implement pin power reconstruction.
    """
    Currently supports cycle length, F_q, F_dh, and max boron.
    
    Updated by Nicholas Rollins. 09/27/2024
    """
    ## Prepare container for results
    results_dict = {}
    for res in ["cycle_length", "PinPowerPeaking", "FDeltaH", "max_boron"]:
        results_dict[res]['value'] = []
    
    ## Read file for parsing
    with open(filename + ".parcs_dpl", "r") as ofile:
        filestr = ofile.read()
    
    ## Split file by section
    res_str = filestr.split('===============================================================================')
    res_str = res_str[1].split('-------------------------------------------------------------------------------')
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
    results_dict["PinPowerPeaking"]["value"] = max(fq_list)
    results_dict["FDeltaH"]["value"] = max(fdh_list)
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
                if mboron > new_max_boron:
                    new_max_boron = mboron
        results_dict["max_boron"]["value"] = new_max_boron
    
    for param in parameters.keys():
        if param in results_dict:
            parameters[param]['value'] = results_dict[param]["value"]
    
    return parameters

def evaluate(solution, input):
    """
    #!TODO: write docstring.
    
    Updated by Nicholas Rollins. 10/03/2024
    """
## Create and move to unique directory for PARCS execution
    cwd = Path(os.getcwd())
    indv_dir = cwd.joinpath(input.results_dir_name + '/' + solution.name)
    os.mkdir(indv_dir)
    os.chdir(indv_dir)

## Prepare depletion file template #!TODO: can this file be dynamically generated instead of copied?
    if input.map_size == 'quarter':
        if input.number_assemblies == 193:
            shutil.copyfile('/home/nkrollin/midas/MIDAS/xslib/' + 'boc_exp_quart193_18.dep', 'boc_exp.dep') #!TODO: change this path to global variable
    else: #assume full geometry if not quarter-core
        if input.number_assemblies == 193:
            shutil.copyfile('/home/nkrollin/midas/MIDAS/xslib/' + 'boc_exp_full193.dep', 'boc_exp.dep')
        elif input.number_assemblies == 157:
            shutil.copyfile('/home/nkrollin/midas/MIDAS/xslib/' + 'boc_exp_full157.dep', 'boc_exp.dep')
    
## Fill loading pattern with chromosome
    fuel_locations = list(input.core_dict['fuel'].keys())
    soln_loading_pattern = {}
    for i in range(len(solution.chromosome)):
        soln_loading_pattern[fuel_locations[i]] = solution.chromosome[i]
    
    raise InputError("DEBUG STOP.")#!
    #!input.core_lattice
    
## Generate Input File
    filename = solution.name + '.inp'
    
    ## CaseID Block
    with open(filename,"w") as ofile:             
        ofile.write("!******************************************************************************\n")
        ofile.write('CASEID {}  \n'.format(solution.name))
        ofile.write("!******************************************************************************\n\n")

    ## CNTL Block
    with open(filename,"a") as ofile:             
        ofile.write("CNTL\n")
        ofile.write("     RUN_OPTS F T F F\n")
        ofile.write("     TH_FDBK    T\n")
        ofile.write("     INT_TH     T -1\n")
        ofile.write("     CORE_POWER 100.0\n")
        ofile.write("     CORE_TYPE  PWR\n")
        ofile.write("     PPM        1000 1.0 1800.0 10.0\n")
        ofile.write("     DEPLETION  T  1.0E-5 T\n")
        ofile.write("     TREE_XS    T  {} T  T  F  F  T  F  T  F  T  F  T  T  T  F \n".format(int(len(xs_unique)+4)))
        ofile.write("     BANK_POS   100 100 100 100 100 100\n")
        ofile.write("     XE_SM      1 1 1 1\n")
        ofile.write("     SEARCH     PPM\n")
        ofile.write("     XS_EXTRAP  1.0 0.3\n")
        ofile.write("     PIN_POWER  T\n")
        ofile.write("     PLOT_OPTS 0 0 0 0 0 2\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
        
    ## PARAM Block
    with open(filename,"a") as ofile:             
        ofile.write("PARAM\n")
        ofile.write("     LSOLVER  1 1 20\n")
        ofile.write("     NODAL_KERN     NEMMG\n")
        ofile.write("     CMFD     2\n")
        ofile.write("     DECUSP   2\n")
        ofile.write("     INIT_GUESS 0\n")
        ofile.write("     CONV_SS   1.e-6 5.e-5 1.e-3 0.001\n")
        ofile.write("     EPS_ERF   0.010\n")
        ofile.write("     EPS_ANM   0.000001\n")
        ofile.write("     NLUPD_SS  5 5 1\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")
    
    ## GEOM Block
    with open(filename,"a") as ofile:             
        ofile.write("GEOM\n")
        ofile.write("     GEO_DIM 9 9 18 1 1\n")
        ofile.write("     RAD_CONF\n")
        for x in range(self.core_lattice.shape[0]):
            ofile.write("     ")
            for y in range(self.core_lattice.shape[1]):
                ofile.write(self.core_lattice[x,y])
                ofile.write("  ")
            ofile.write("\n")
    
        ofile.write("     GRID_X      1*10.75 8*21.50\n")
        ofile.write("     NEUTMESH_X  1*1 8*1\n")
        ofile.write("     GRID_Y      1*10.75 8*21.50\n")
        ofile.write("     NEUTMESH_Y  1*1 8*1\n")
        ofile.write("     GRID_Z      30.48 15.24 10.16 5.08 10*30.48 5.08 10.16 15.24 30.48\n")            
        ofile.write("     ASSY_TYPE   10   1*2   16*2    1*2 REFL\n")
        for i in range(xs_unique.shape[0]):
            if 'gd_0' in xs_unique[i]:
                ofile.write("     ASSY_TYPE   {}   1*1 1*4  14*{}  1*4  1*3 FUEL\n".format(tag_unique[i],xs_ref[i]))
            else:
                ofile.write("     ASSY_TYPE   {}   1*1 1*4  1*4 12*{} 1*4 1*4  1*3 FUEL\n".format(tag_unique[i],xs_ref[i]))
        ofile.write("\n")

        ofile.write("     boun_cond   0 2 0 2 2 2\n")
        ofile.write("     SYMMETRY 4\n")

        ofile.write("     PINCAL_LOC\n")
        for x in range(pincal_loc.shape[0]):
            ofile.write("      ")
            for y in range(pincal_loc.shape[1]):
                val = pincal_loc[x,y]
                if np.isnan(val):
                    pass
                else:
                    ofile.write(str(int(pincal_loc[x,y])))
                    ofile.write("  ")
            ofile.write("\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## FDBK Block
    with open(filename,"a") as ofile:             
        ofile.write("FDBK\n")
        ofile.write("     FA_POWPIT     {} 21.5\n".format(np.round(self.power/193,4)))
        ofile.write("     GAMMA_FRAC    0.0208    0.0    0.0\n")
        ofile.write("     EFF_DOPLT   T  0.5556\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## TH Block
    with open(filename,"a") as ofile:   
        ofile.write("TH\n")          
        ofile.write("     FLU_TYP       0\n")
        ofile.write("     N_PINGT    264 25\n")
        ofile.write("     PIN_DIM      4.1 4.75 0.58 6.13\n")
        ofile.write("     FLOW_COND    {}  {}\n".format(np.round(self.inlet_temperature-273.15,2),np.round(self.flow/193,4)))
        ofile.write("     HGAP     11356.0\n")
        ofile.write("     N_RING   6\n")
        ofile.write("     THMESH_X       9*1\n")
        ofile.write("     THMESH_Y       9*1\n")
        ofile.write("     THMESH_Z       1 2 3 4 5 6 7 8 9 10 11 12\n")
        ofile.write("\n")
        ofile.write("!******************************************************************************\n\n")

    ## DEPL Block
    with open(filename,"a") as ofile:             
        ofile.write("DEPL\n")
        ofile.write("     TIME_STP  1 1 14*30\n")
        ofile.write("     INP_HST   './boc_exp.dep' -2 1\n")
        ofile.write("     PMAXS_F   1 '{}' 1\n".format(cdir + '/' + 'xs_gbot'))
        ofile.write("     PMAXS_F   2 '{}' 2\n".format(cdir + '/' + 'xs_grad'))
        ofile.write("     PMAXS_F   3 '{}' 3\n".format(cdir + '/' + 'xs_gtop'))
        ofile.write("     PMAXS_F   4 '{}' 4\n".format(cdir + '/' + 'xs_g250_gd_0_wt_0'))
        for i in range(xs_unique.shape[0]):
            ofile.write("     PMAXS_F   {} '{}' {}\n".format(5+i,cdir + '/' + xs_unique[i],5+i))
        ofile.write("\n")
        ofile.write(".")

## Run PARCS INPUT DECK
    parcscmd = "/cm/shared/codes/TRACE51341_PARCS_332/PARCS-v332_Exe/Executables/Linux/parcs-v332-linux2-intel-x64-release.x" #!TODO: how to move this to a global or environmental variable?
    try:
        output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=50) #wait until calculation finishes
    ## Get Results
        if 'Finished' in str(output): #job completed
            ofile = solution.name + '.out'
            solution.parameters = get_results(solution.parameters, solution.name)
        else: #job failed
            pass #!TODO: add logging and failure step
    except subprocess.TimeoutExpired: #job timed out
        os.system('rm -f {}.parcs_pin*'.format(solution.name))
        #!TODO: add logging and failure step
    
    os.chdir(cwd)
    gc.collect()
    
    return solution.parameters

def calc_cycle_length(efpd,boron,keff):
    if boron[-1]==0.1:
        eoc1_ind = 0
        eco2_ind = len(efpd)
        for i in range(len(efpd)):
            if boron[i] > 0.1 and boron[i+1] == 0.1:
                eoc1_ind = i
                eco2_ind = i+1
        dbor = abs(boron[eoc1_ind-1]-boron[eoc1_ind])
        defpd = abs(efpd[eoc1_ind-1]-efpd[eoc1_ind])
        def_dbor = defpd/dbor
        eoc = efpd[eoc1_ind] + def_dbor*(boron[eoc1_ind]-0.1)
    elif boron[-1]==boron[0]==1800.0:
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