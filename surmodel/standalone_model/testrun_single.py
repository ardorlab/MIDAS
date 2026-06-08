import multiprocessing as mp
import os
import pickle
import numpy as np
import time as TT

from test import get_result, interpolatecycle
# os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

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
    # cycle = calc_cycle_length(efpd_list,boron_list,keff_list)
    ## get the index 
    arr = np.array(boron_list)
    idx = np.where(np.isclose(arr, 0.1))[0]
    index = idx[0] if idx.size > 0 else -1

    cycle = interpolatecycle(10,np.array(boron_list[:index]).reshape(-1,1),np.array(efpd_list[:index]).reshape(-1,1))
    cbc = max(boron_list)
    fq = max(fq_list)
    fd = max(fdh_list)
    # print(boron_list, efpd_list)
    # print(cycle)
    
    return index, np.array([cycle,cbc,fq,fd])

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



if __name__ == "__main__":
    # 1. Load shared data once in main process
    xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/updateXSMay26.pkl','rb'))
    listFA = [251,252,321,322,201,461,462,501,502,526,566,586,250,280,320,400,567,587]
    # combo_generator = itertools.product(listFA, repeat=5)
    corebulist = [0.1, 0.4, 0.5,1.0,1.0]
    initbumap = np.zeros((1,81,16))
    for i in range(30):
        corebulist.append(1*1.0)
    # corebulist.append(0.84)
    # Start timing
    start_time = TT.time()

    Fqmax = []
    Fdmax = []
    boronmax_pred = []
    cycle_length_pred = []
    boronmax_true = []
    cycle_length_true = []
    count = 0
    bustep = 0
    ## your  testing LPs
    LPs = ['252  252  201  201  201  251  252  322  10               \n', 
           '252  201  201  251  251  252  251  201  10               \n', 
           '201  201  251  321  322  322  322  10  10            \n', 
           '201  321  251  322  201  251  201  10  00            \n', 
           '201  252  252  322  252  321  10  10  00       \n',
           '251  201  321  321  251  10  10  00  00    \n',
           '252  321  321  321  10  10  00  00  00     \n',
           '322  252  10  10  10  00  00  00  00   \n',
           '10  10  10  00  00  00  00  00  00   \n']

    LPs_o = " ".join(LPs).strip().split()
    fd,fq,cbc,cycle=get_result(LPs_o,bustep=bustep, count=count, 
               cycle_length_pred=cycle_length_pred, 
               boronmax_pred = boronmax_pred,
               Fqmax=Fqmax, Fdmax= Fdmax, corebulist = corebulist, listFA=listFA, index=-1, initbumap=initbumap)
    print(f"Cycle length (EFPD): {cycle:.3f}")
    print(f"Peak boron concentration (ppm): {cbc:.3f}")
    print(f"Peak Fq: {fq:.3f}")
    print(f"Peak Fd: {fd:.3f}")
    
    print("Calculation complete (s):", TT.time() -start_time)
