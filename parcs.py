import os
import gc
import sys
import copy
import h5py
import math
import time
import numpy as np
import pickle
import random
from solution_types import Solution
from pathlib import Path
import shutil 
from lcoe import LCOE

def get_clength(efpd,boron):
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
    else:
        dbor = abs(boron[-2]-boron[-1])
        defpd = abs(efpd[-2]-efpd[-1])
        def_dbor = defpd/dbor
        eoc = efpd[-1] + def_dbor*(boron[-1]-0.1)
    return(eoc)

def get_pin_power(filepath):
    start = time.time()
    print('Reading of Pin Powers')
    ofile = open(filepath+".pin", "r")
    filestr = ofile.read()
    ofile.close()
    bustr=filestr.split("At Time:")
    nbu= len(bustr)-1
    nz_tag = bustr[1][11:48]
    nass_tag = bustr[1][49:86]
    nass = len(bustr[1].split(nass_tag))-1
    nz = len(bustr[1].split(nz_tag))-2
    asb1_str = bustr[1].split(nz_tag)
    npx = int(asb1_str[1].split("\n")[1].split()[-1])
    npy = len(asb1_str[1].split("\n")) -4
    pp_mat = np.zeros((nbu,nass,nz,npx*npy))
    for ibu in range(nbu):
        iass=0
        iz=0
        ibustr = bustr[ibu+1]
        iasszstr = ibustr.split('Assembly Coordinate (i,j):')
        for j in range(1,len(iasszstr)):
            iassz = iasszstr[j].split('\n')
            ass_id = np.array([int(iassz[0][0:5]),int(iassz[0][5:9])])
            z_id = int(iassz[0][35:38])
            if z_id == 0:
                iz=0
                iass+=1
                continue
            else:
                pp_str = iassz[2:2+npy]
                for iy in range(npy):
                    for ix in range(npx):
                        pp_id = iy*npx + ix 
                        try:
                            pp_mat[ibu,iass,iz,pp_id] = float(pp_str[iy][(7*ix + 8):(7*ix + 14)])
                        except:
                            print("Non physical peaking factors")
                            pp_mat[ibu,iass,iz,pp_id] = 10.0
                iz+=1
    end = time.time()
    print('Pin Power Duration = {} s'.format(end-start))
    return(pp_mat)

def get_asb_power(filepath):
    start = time.time()
    print('Reading of Assembly Powers')
    ofile = open(filepath+".dep", "r")
    filestr = ofile.read()
    ofile.close()
    bustr=filestr.split(" RPF 3D MAP")
    nbu= len(bustr)-1
    nasb_str = filestr.split('===============================================================================')
    nass=int(nasb_str[0].split()[-1])
    nz_str = bustr[0].split(' RPF 1D MAP')[1].split('\n')
    nrefl=2
    nz = len(nz_str)-4-nrefl
    ztag = np.arange(2,nz+1 + 1)
    asb_mat = np.zeros((nbu,nass,nz))
    for ibu in range(nbu):
        ibustr = bustr[ibu+1].split(' EXP 2D MAP')[0]
        asb_str = ibustr.split(' k lth')
        iasb=0
        for ik in range(1,len(asb_str)):
            asb_line=asb_str[ik].split('\n')
            for iz in range(1,len(asb_line)):
                asb_val=asb_line[iz].split()
                if len(asb_val)>0:
                    if int(asb_val[0]) in ztag:
                        zid = int(asb_val[0])-2
                        asb_count = 0
                        for ia in range(1,len(asb_val)):
                            val = float(asb_val[ia])
                            if  val !=0.0:
                                asb_id = iasb + asb_count
                                asb_mat[ibu,asb_id,zid]=val
                                asb_count+=1
                            else:
                                continue
                else:
                    continue
            iasb += asb_count
    end = time.time()
    print('Assembly Power Duration = {} s'.format(end-start))
    return(asb_mat)

def get_lcoe(cycle_length,full_core,core_dict):

# Compute number of fuel assemblies

    nfuel=0
    nfuel_sym=0
    nrefl = 0
    nrefl_sym = 0
    for key, value in core_dict['core'].items():
        symmetry_multiplier = len(core_dict['core'][key]['Symmetric_Assemblies'])+1
        if key is None:
            continue
        elif key[0]=='R':
            nrefl_sym+=1
            nrefl+=symmetry_multiplier
        else:
            nfuel_sym+=1
            nfuel+=symmetry_multiplier
    
    cycle_param={'EFPD': cycle_length,
                'Batches': 3,
                'Thermal_Power': core_dict['Parameters']['Thermal_Power'],
                'Efficiency': 0.33,
                'Fuel_Assemblies': nfuel}

    lcoe_param={'Discount_Rate': 0.07,
                'Uranium_Ore_Price': 80,
                'Conversion_Price': 10,
                'Enrichment_Price': 160,
                'Fabrication_Price': 250,
                'Uranium_Ore_Loss': 0.002,
                'Conversion_Loss': 0.002,
                'Enrichment_Loss': 0.002,
                'Fabrication_Loss': 0.002,
                'Enrichment_Feed': 0.00711,
                'Enrichment_Tail': 0.003,
                'Storage_Price': 200,
                'Disposal_Price': 463,
                'Uranium_Ore_Time': -2.0,
                'Conversion_Time': -1.5,
                'Enrichment_Time': -1.0,
                'Fabrication_Time': -0.5,
                'Storage_Time': 5.0+cycle_param['EFPD']*cycle_param['Batches']/365.25,
                'Disposal_Time': cycle_param['EFPD']*cycle_param['Batches']/365.25}

    unique_fa =  np.unique(list(full_core.values()))
    asb_param = {}
    for i in range(len(unique_fa)):
        nfa = list(full_core.values()).count(unique_fa[i])
        enr = float(unique_fa[i][2:5])/10000
        asb_dict = {'Number': nfa,
                    'Fuel_Rods': 264,
                    'Fuel_Radius': 0.41,
                    'Fuel_Height': 365.76,
                    'Enrichment': enr,
                    'Fuel_Density': 10.23,
                    'Fabrication_Price': 250
                    }
        asb_param[unique_fa[i]]=asb_dict

    lcoe, bu, asb_cost = LCOE(cycle_param,lcoe_param, asb_param)
    asb_cost_dict = {}
    for i in range(len(unique_fa)):
        asb_cost_dict[unique_fa[i]]=asb_cost[i]

    return((lcoe, bu, asb_cost_dict))

def get_results(filepath,full_core,core_dict,pin_power=False):
    efpd=[]
    boron =[]
    fq=[]
    fdh=[]
    read_bool  = False
    with open(filepath + ".dep","r") as ifile:             
        for line in ifile:
            s=line.split()
            if(len(s)>0):
                if(s[0] =="summary:"):
                    read_bool=True
                    continue
            if(len(s)>1 and len(s)>20):
                if read_bool and s[1]=="1":
                    efpd.append(float(s[2]))
                    boron.append(float(s[12]))
                    fq.append(float(s[4].split('(')[0]))
                    fdh.append(float(s[7]))
            elif(len(s)>1 and len(s)==20):
                if read_bool and s[1]=="1":
                    efpd.append(float(s[2]))
                    boron.append(float(s[11]))
                    fq.append(float(s[4].split('(')[0]))
                    fdh.append(float(s[6]))
    
    res = {}
    res["Cycle_Length"] = get_clength(efpd,boron)
    res["Fq"] = max(fq)
    res["Fdh"] = max(fdh)
    res["Max_Boron"] = max(boron)
    lcoe, discharge_bu, asb_cost = get_lcoe(res["Cycle_Length"],full_core,core_dict)
    res["LCOE"] = lcoe
    res["Discharge_Burnup"]=discharge_bu
    res["Assemblies_Costs"] = asb_cost
    if pin_power:
        zh = np.array([15.24, 10.16, 5.08, 30.48, 30.48, 30.48, 30.48, 30.48,
                       30.48, 30.48, 30.48, 30.48, 30.48, 5.08, 10.16, 15.24])
        asb_mat=get_asb_power(filepath)
        pp_mat = get_pin_power(filepath)
        fq_asb = np.max(asb_mat)
        fdh_asb = 0
        for ibu in range(asb_mat.shape[0]):
            for iasb in range(asb_mat.shape[1]):
                fdh_i = np.dot(asb_mat[ibu,iasb,:],zh)/np.sum(zh)
                if fdh_i > fdh_asb:
                    fdh_asb = fdh_i
        fq_pp = 0
        fdh_pp = 0
        fq_id = np.array([0,0])
        fdh_id = np.array([0,0])
        for ibu in range(asb_mat.shape[0]):
            for iasb in range(asb_mat.shape[1]):
                iasb_mat = np.zeros((pp_mat.shape[2],pp_mat.shape[3]))
                for iz in range(pp_mat.shape[2]):
                    iasb_mat[iz,:]=pp_mat[ibu,iasb,iz,:]*asb_mat[ibu,iasb,iz]
                fq_i = np.max(iasb_mat)
                if fq_i > fq_pp:
                    fq_pp = fq_i
                    fq_id[0]=ibu 
                    fq_id[1]=iasb
                for ip in range(pp_mat.shape[3]):
                    fdh_i = np.dot(iasb_mat[:,ip],zh)/np.sum(zh)
                    if fdh_i > fdh_pp:
                        fdh_pp = fdh_i
                        fdh_id[0]=ibu 
                        fdh_id[1]=iasb
        res["Fq"] = fq_pp
        res["Fdh"] = fdh_pp
    return(res)

def get_full_core(core_dict):
        """
        Generates the 17x17 full fuel core from symmetry.

        Parameters: None
    
        Written by Gregory Delipei 7/24/2022
        """
        full_core  = {}
        for key, value in core_dict['fuel'].items():
            full_core[key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core[skey]=value['Value'] 
        return(full_core)

def get_quarter_lattice(core_dict):
        """
        Generates the 17x17 quarter core lattice.

        Parameters: None
    
        Written by Gregory Delipei 7/24/2022
        """
        core_map = core_dict['core_map']
        core_id = core_dict['core_id']
        nx = int(np.ceil(len(core_map[0])/2))
        ny = int(np.ceil(len(core_map)/2))
        quarter_core = np.zeros((ny,nx), dtype='<U8')
        for y in range(ny):
            for x in range(nx):
                val =core_map[ny-1+y,nx-1+x]
                if val is None:
                    val = "00"
                quarter_core[y,x] = val
        return(quarter_core)

def evaluate(fname,ldir,core_dict):
        """
        Creates the input deck, runs the calculation and retrieves the results and the cost.

        Parameters: 
        loc: String - Directory of execution
        fname: String - File name

        Written by Gregory Delipei 7/29/2022
        """

        # Create PARCS INPUT DECK

        pwd = Path(os.getcwd())

        if not os.path.exists(ldir):
            os.makedirs(ldir)
        else:
            shutil.rmtree(ldir, ignore_errors=True)
            os.makedirs(ldir)
        
        cdir = core_dict['Parameters']['xslib']
        shutil.copyfile(cdir + '/' + 'boc_exp_quart193.dep', ldir +"/" + 'boc_exp.dep')
        os.chdir(ldir)

        full_core = get_full_core(core_dict)
        qcore_lattice = get_quarter_lattice(core_dict)

        xs_array = np.zeros((qcore_lattice.shape[0],qcore_lattice.shape[1]), dtype='<U20')
        pincal_loc = np.zeros((qcore_lattice.shape[0],qcore_lattice.shape[1]))
        for x in range(qcore_lattice.shape[0]):
            for y in range(qcore_lattice.shape[1]):
                loc = qcore_lattice[x,y]
                if loc != "00" and loc[0] != "R":
                    qcore_lattice[x,y] = core_dict['Inventory'][full_core[loc]]['Tag']
                    xs_array[x,y] = core_dict['Inventory'][full_core[loc]]['Cross_Section']
                    pincal_loc[x,y]=1
                elif loc[0] == "R":
                    qcore_lattice[x,y] = "10 "
                    xs_array[x,y] = None
                    pincal_loc[x,y]=0
                elif loc == "00":
                    qcore_lattice[x,y] = "00 "
                    xs_array[x,y] = None
                    pincal_loc[x,y]=np.nan

        xs_unique = np.unique(xs_array)
        xs_unique = np.delete(xs_unique, np.argwhere(xs_unique == 'None'))
        xslib = core_dict["Parameters"]['xslib']
        tag_unique = copy.deepcopy(xs_unique)
        xs_ref = np.arange(5,5+len(xs_unique)) # 1-3 for reflectors and 4 for blankets
        for key,value in core_dict["Inventory"].items():
            for i in range(xs_unique.shape[0]):
                if value['Cross_Section']==xs_unique[i]:
                    tag_unique[i]=value['Tag']
       
        filename = fname+'.inp'
        with open(filename,"w") as ofile:             
            ofile.write("!******************************************************************************\n")
            ofile.write('CASEID pwr_193  \n')
            ofile.write("!******************************************************************************\n\n")

        with open(filename,"a") as ofile:             
            ofile.write("CNTL\n")
            ofile.write("     TH_FDBK    T\n")
            ofile.write("     CORE_POWER 100\n")
            ofile.write("     CORE_TYPE  PWR\n")
            ofile.write("     PPM        1000\n")
            ofile.write("     DEPLETION  T  1.0E-5 T\n")
            ofile.write("     TREE_XS    T  {} T  T  F  F  T  F  T  F  T  F  T  T  T  F \n".format(int(len(xs_unique)+4)))
            ofile.write("     BANK_POS   100 100 100 100 100 100\n")
            ofile.write("     XE_SM      1 1\n")
            ofile.write("     SEARCH     PPM\n")
            ofile.write("     XS_EXTRAP  1.0 0.3 0.8 0.2\n")
            ofile.write("     PIN_POWER  T\n")
            ofile.write("     PRINT_OPT T F T T F T F F F T  F  F  F  F  F\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
            
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
             
        with open(filename,"a") as ofile:             
            ofile.write("GEOM\n")
            ofile.write("     GEO_DIM 9 9 18 1 1\n")
            ofile.write("     COORDINATE  1\n")
            ofile.write("     SYMMETRY 4\n")
            ofile.write("     RAD_CONF\n")
            for x in range(qcore_lattice.shape[0]):
                ofile.write("     ")
                for y in range(qcore_lattice.shape[1]):
                    ofile.write(qcore_lattice[x,y])
                    ofile.write("  ")
                ofile.write("\n")
            ofile.write("     GRID_X      1*10.75 8*21.50\n")
            ofile.write("     NEUTMESH_X  1*1 8*1\n")
            ofile.write("     GRID_Y      1*10.75 8*21.50\n")
            ofile.write("     NEUTMESH_Y  1*1 8*1\n")
            ofile.write("     GRID_Z      30.48 15.24 10.16 5.08 10*30.48 5.08 10.16 15.24 30.48\n")
            ofile.write("     Boun_cond   0 2 0 2 2 2\n")
            
            ofile.write("     ASSY_TYPE   10   1*2   16*2    1*2 REFL\n")
            for i in range(xs_unique.shape[0]):
                if 'gd_0' in xs_unique[i]:
                    ofile.write("     ASSY_TYPE   {}   1*1 1*4  14*{}  1*4  1*3 FUEL\n".format(tag_unique[i],xs_ref[i]))
                else:
                    ofile.write("     ASSY_TYPE   {}   1*1 1*4  1*4 12*{} 1*4 1*4  1*3 FUEL\n".format(tag_unique[i],xs_ref[i]))
            ofile.write("\n")

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
        
        with open(filename,"a") as ofile:             
            ofile.write("TH\n")
            ofile.write("     FA_POWPIT     {} 21.5\n".format(np.round(core_dict['Parameters']['Thermal_Power']/193,4)))
            ofile.write("     GAMMA_FRAC    0.0208    0.0    0.0\n")
            ofile.write("     FLU_TYP       0\n")
            ofile.write("     N_PINGT    264 25\n")
            ofile.write("     PIN_DIM      4.1 4.75 0.58 6.13\n")
            ofile.write("     FLOW_COND    {}  {}\n".format(np.round(core_dict['Parameters']['Inlet_Temperature']-273.15,2),np.round(core_dict['Parameters']['Core_Flow']/193,4)))
            ofile.write("     HGAP     11356.0\n")
            ofile.write("     N_RING   6\n")
            ofile.write("     THMESH_X       9*1\n")
            ofile.write("     THMESH_Y       9*1\n")
            ofile.write("     THMESH_Z       1 2 3 4 5 6 7 8 9 10 11 12\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")

        with open(filename,"a") as ofile:             
            ofile.write("DEPL\n")
            ofile.write("     TIME_STP  1 1 14*30\n")
            ofile.write("     INP_HST   './boc_exp.dep' -2 1\n")
            ofile.write("     PMAXS_F   1 '{}' 1\n".format(xslib+ '/' + 'xs_gbot'))
            ofile.write("     PMAXS_F   2 '{}' 2\n".format(xslib+ '/' + 'xs_grad'))
            ofile.write("     PMAXS_F   3 '{}' 3\n".format(xslib+ '/' + 'xs_gtop'))
            ofile.write("     PMAXS_F   4 '{}' 4\n".format(xslib+ '/' + 'xs_g250_gd_0_wt_0'))
            for i in range(xs_unique.shape[0]):
                ofile.write("     PMAXS_F   {} '{}' {}\n".format(5+i,xslib+ '/' + xs_unique[i],5+i))
            ofile.write("\n")
            ofile.write(".")

        # Run PARCS INPUT DECK
        parcscmd = "/cm/shared/apps/ncsu/TRACE_PARCS/PARCS-v32m19co_src/run/parcs_ifr.x " + filename
        os.system(parcscmd)
        # Get Results
        ofile = fname + '.out'
        res=get_results('pwr_193',full_core,core_dict, pin_power=True)
        os.chdir(pwd)
        return(res)


