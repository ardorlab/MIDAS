import multiprocessing as mp
import os
import pickle
import numpy as np
import time as TT
import copy
import re

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

def get_map(inputfile, key):
    with open(inputfile,'r') as f:
        txt = f.readlines()
    for i,line in enumerate(txt):
        if line.find(key)>=0:
            st = i+1 
            coresize = len(txt[st].strip().split())
            break
    nmap=''
    for ii in range(st, st+coresize):
        temp = txt[ii]
        parts = temp.strip().split()
        linenew = ' '.join(parts + ['00'] * (9 - len(parts)))
        nmap+=linenew+'\n'
    
    resultmap = re.sub(r'(?<![0-9-])0(?![0-9x])', '00', nmap)
    nmap = "".join(resultmap).strip().split()
    return nmap 

def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_initbumap(initdepfile,inputfile, outputfile):
    initbu3dmap = get_burnup_17(initdepfile)[-1] ## 
    initburnup = np.zeros((1,9,9,16))
    locationmap = np.array(get_map(inputfile, '    LOCATION   0')).reshape(9,9)
    rlmap = np.array(get_map(inputfile, '    SHUF_MAP   1   1')).reshape(9,9)
    ## depleted map 
    depleted_map = getLP(outputfile)
    depleted_map = np.array(depleted_map).reshape(9,9)
    ## indexing map
    idx_map = []
    idx=1
    for i in range(depleted_map.shape[0]):
        idx_map.append([])
        for j in range(depleted_map.shape[1]):
            if depleted_map[i,j] not in ['10','00']:
                idx_map[i].append(idx)
                idx+=1
            else:
                idx_map[i].append(0)
    idx_map = np.array(idx_map)
    shufmap = copy.deepcopy(rlmap)
    for x in range(depleted_map.shape[0]):
        for y in range(depleted_map.shape[1]):
            if not is_integer(rlmap[x,y]):
                idx =  np.where(locationmap == rlmap[x,y])
                if len(idx[0]) >0:
                    shufmap[x,y] = idx_map[idx[0][0], idx[1][0]]
    ## assign initbumap
    for row in range(shufmap.shape[0]):
        for col in range(shufmap.shape[1]):
            faidx = int(shufmap[row][col])
            # print(faidx)
            if faidx >0:
                initburnup[:,row,col,:]=initbu3dmap[faidx-1,:]
    ## reshape 
    initburnup = initburnup.reshape((1,81,16))
    return initburnup

if __name__ == "__main__":
    # 1. Load shared data once in main process
    xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/updateXS.pkl','rb'))
    listFA = [461,462,501,502,526,566,586,250,280,320,400,567,587]
    # combo_generator = itertools.product(listFA, repeat=5)
    corebulist = [0.1, 0.4, 0.5,1.0,1.0]
    for i in range(28):
        corebulist.append(1*1.0)
    # Start timing
    initbumap = np.zeros((1,81,16))
    start_time = TT.time()
    ## run all the loop 
    pathtotest = '/home/khnguy22/Deeponet-midas/MIDAS/corereload197-surrogate-nofb_new/output_files'
    initbufile = '/home/khnguy22/Deeponet-midas/MIDAS/corereload197-surrogate-nofb/BOC_exposure.dep'
    preddata = []
    truedata = []
    ## run the loop 
    for root, dirs, files in os.walk(pathtotest):
        for file in files:
            if file.endswith('.inp'):
                inputfile = os.path.join(root, file)
                ## get initbumap 
                # depfilename = file.replace('.inp', '.parcs_dep')
                depfilename = file.replace('.inp', '.parcs_cyc-02') # for multi_cycle case
                pinfilename = file.replace('.inp','.parcs_pin')
                dplfilename = file.replace('.inp','.parcs_dpl')
                outfilename = file.replace('.inp','.parcs_out')
                print(f"--- Processing for {inputfile} ---")
                folder_path = root
                outputdep_path = os.path.join(folder_path, depfilename)
                outputpin_path = os.path.join(folder_path, pinfilename)
                outputdpl_path = os.path.join(folder_path, dplfilename)
                output_path = os.path.join(folder_path, outfilename)
                initbumap = get_initbumap(initbufile,inputfile,output_path)
                index, array = get_core_data(outputdpl_path)
                truedata.append(array)
                # test = getLP(output_path)
                test = getLPs_cyc(outputdep_path)  # for multi_cycle
                Fqmax = []
                Fdmax = []
                boronmax_pred = []
                cycle_length_pred = []
                boronmax_true = []
                cycle_length_true = []
                count = 0
                bustep = 0
    
                # 2. Define list of Loading Patterns (LPs)
                # LPs = ['280	250	567	567	320	567	587	587	10 \n', 
                #        '250	250	567	567	567	567	587	587	10 \n', 
                #        '567	567	250	567	250	587	587	587	10 \n', 
                #        '567	567	567	567	567	587	587	587	10 \n', 
                #        '320	567	250	567	320	587	587	10	10 \n',
                #        '567	567	587	587	587	587	400	10	00 \n',
                #        '587	587	587	587	587	400	10	10	00\n',
                #        '587	587	587	587	10	10	10	00	00\n',
                #        '10	10	10	10	10	00	00	00	00\n']
            
                # LPs_o = " ".join(LPs).strip().split()
                # print(LPs_o)
                # list_of_LPs = [ LPs_o] 
                LPs_o = test
        
                fd,fq,cbc,cycle=get_result(test,bustep=bustep, count=count, 
                           cycle_length_pred=cycle_length_pred, 
                           boronmax_pred = boronmax_pred,
                           Fqmax=Fqmax, Fdmax= Fdmax, corebulist = corebulist, listFA=listFA, index=index, initbumap=initbumap)
                
                if cbc >2000:
                    print('waaaaaaa',fd,fq,cbc,cycle)
                    print('trueee', array)
                    stop
                preddata.append(np.array([cycle,cbc,fq,fd]))
                print('pred', cycle,cbc,fq,fd)
                print('true', array )
                
                print("Parallel processing complete:", inputfile)
np.save('preddata_51.npy',np.array(preddata))
np.save('truedata_51.npy',np.array(truedata))
print(f'average time cost {1/500*(TT.time()-start_time):.3f}')



# from test import get_result
# # # import Pool 

# LPs = ['526  526  526  566  501  526  526  501  10 \n', 
#        '526  526  462  502  526  566  501  501  10 \n', 
#        '462  526  566  501  526  502  526  501  10 \n', 
#        '566  462  526  526  501  526  462  501  10 \n', 
#        '526  501  566  526  501  526  566  10   10 \n',
#        '462  566  526  502  462  526  501  10   00 \n',
#        '526  462  501  462  526  526  10   10   00 \n',
#        '526  526  501  501  10   10   10   00   00 \n',
#        '10   10   10   10   10   00   00   00   00 \n']

# LPs_o = " ".join(LPs).strip().split()
# print(LPs)
# LPs = ["526","526","526" ,"566","501","526","526", "10", "10",
# "526","526","462" ,"502","526","566","501", "10", "10",
# "462","526","566" ,"501","526","502","526", "10", "10",
# "566","462","526" ,"526","501","526","462", "10", "10",
# "526","501","566" ,"526","501","526","566", "10", "10",
# "462","566","526" ,"502","462","526","501", "10", "00",
# "526","462","501" ,"462","526","526", "10", "10", "00",
# "526","462","462" ,"502", "10", "10", "10", "00", "00",
#  "10", "10", "10" , "10", "10", "00", "00", "00", "00",]

# for i in range(81):
#     print(LPs[i], LPs_o[i])


# stop
# print(get_result(LPs))