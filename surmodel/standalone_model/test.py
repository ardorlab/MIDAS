'''
Create depletion module to calculate the burn up 

'''

import os
import shutil
import random
import subprocess
import itertools
from pathlib import Path
import pickle 
import numpy as np
import copy
import time as TT
import matplotlib.pyplot as plt
st = TT.time()

from loadmodel_Nov25 import minmaxReverse, minmaxscale
from loadmodel_Nov25 import model
from loadmodel_Nov25 import scalerparam
from loadmodel_Nov25 import trunks
from loadmodelcore import trunks as trunk_core
from loadmodelcore import scalerparam as scaler_core
from loadmodelcore import model as modelcore
from loadmodelpin import best_model as pinmodel
# from loadpinmodel_new import best_model as pinmodel

import tensorflow as tf
import os

# 1. Get the number of physical cores (not logical threads/hyperthreading)
# If you know you have 32 cores, hardcode this to 32.
num_cores = os.cpu_count()  
# Often better to use physical cores only, so try os.cpu_count() // 2 if using Hyperthreading

# 2. Set environment variables (do this BEFORE any other TF code)
os.environ["OMP_NUM_THREADS"] = str(num_cores)
os.environ["TF_NUM_INTRAOP_THREADS"] = str(num_cores)
os.environ["TF_NUM_INTEROP_THREADS"] = "1" # Usually 1 or 2 is best for sequential inference

# 3. Configure TensorFlow
tf.config.threading.set_intra_op_parallelism_threads(num_cores)
tf.config.threading.set_inter_op_parallelism_threads(1)


## function to get EXPOSURE each step 
# testdatapath ='./rawdata12k_coreparam/'
# tempborcyc = np.load(testdatapath+'output.npy').astype(np.float32)
# bor_truth = tempborcyc[:,0]
# cyc_truth = tempborcyc[:,1]
## for comparing
def get_burnup(ofile,FULL_CORE=False):
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

def get_rpf(ofile,FULL_CORE=False):
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


def extract_data_old(solution,input,xsdict):
    """
    Extract data for retraining models 
    DO NOT USE!!!
    :param solution: list of solution 
    :param input: input object form yaml
    """
    cwd = Path(os.getcwd())
    listFA = [461,462,501,502,526,566,586,250,280,320,400,567,587] ## for now
    # xsdict = pickle.load(open(midas_data.__path_xs_pickle__,'rb'))
    # xsdict = joblib.load(midas_data.__path_xs_pickle__, mmap_mode='r')
    fabulist = np.array([0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 
                         17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 
                         60.0, 65.0, 70.0, 75.0, 80.0])
    retrain_path = midas_data.__path_to_store_retrain_data__
    os.makedirs(retrain_path, exist_ok=True)
    branchXS_1 = load_or_empty(os.path.join(retrain_path, 'branchXS_1.npy'),shape=(0,81,16))
    branchXS_2 = load_or_empty(os.path.join(retrain_path, 'branchXS_2.npy'),shape=(0,81,16))
    branchXS_3 = load_or_empty(os.path.join(retrain_path, 'branchXS_3.npy'),shape=(0,81,16))
    branchXS_4 = load_or_empty(os.path.join(retrain_path, 'branchXS_4.npy'),shape=(0,81,16))
    branchXS_5 = load_or_empty(os.path.join(retrain_path, 'branchXS_5.npy'),shape=(0,81,16))
    branchXS_6 = load_or_empty(os.path.join(retrain_path, 'branchXS_6.npy'),shape=(0,81,16))
    branchXS_7 = load_or_empty(os.path.join(retrain_path, 'branchXS_7.npy'),shape=(0,81,16))
    outputRPF  = load_or_empty(os.path.join(retrain_path, 'outputRFP.npy'), shape=(0,81,16))
    outputcore = load_or_empty(os.path.join(retrain_path, 'outputcore.npy'), shape=(0,2))
    list_branchXS_1 = [] 
    list_branchXS_2 = []
    list_branchXS_3 = []
    list_branchXS_4 = []
    list_branchXS_5 = []
    list_branchXS_6 = []
    list_branchXS_7 = []
    list_outputRPF  = []
    list_outputcore = []
    for dir in solution:
        folder_path = cwd.joinpath(input.results_dir_name / Path(dir))
        of_filename = f"{dir}.parcs_dep"
        input_filename = f"{dir}.parcs_out"
        dpl_filename = f"{dir}.parcs_dpl"
        output_path = os.path.join(folder_path, of_filename)
        if os.path.exists(output_path):
            bu3d = get_burnup_17(output_path)
            rpf3d = get_rpf(output_path)
            coreparam = get_boron_cycle(output_path)
            nbu = len(bu3d)
            LPs = getLP(os.path.join(folder_path, input_filename))
            if nbu !=0:
                nz = len(bu3d[0][0])
                for bu in range(nbu):
                    idx = 0 
                    tempbranchXS_1 = np.zeros((81,nz))
                    tempbranchXS_2 = np.zeros((81,nz))
                    tempbranchXS_3 = np.zeros((81,nz))
                    tempbranchXS_4 = np.zeros((81,nz))
                    tempbranchXS_5 = np.zeros((81,nz))
                    tempbranchXS_6 = np.zeros((81,nz))
                    tempbranchXS_7 = np.zeros((81,nz))
                    tempoutput = np.zeros((81,nz))
                    for loc,fa in enumerate(LPs):
                        if fa =='10':
                            # use 0 for reflector 
                            tempoutput[loc][:]=[0.0]*nz
                            tempbranchXS_1[loc][:]=[0.0]*nz
                            tempbranchXS_2[loc][:]=[0.0]*nz
                            tempbranchXS_3[loc][:]=[0.0]*nz
                            tempbranchXS_4[loc][:]=[0.0]*nz
                            tempbranchXS_5[loc][:]=[0.0]*nz
                            tempbranchXS_6[loc][:]=[0.0]*nz
                            tempbranchXS_7[loc][:]=[0.0]*nz
                        elif fa != '00':
                            for node in range(nz):
                                buval = bu3d[bu][idx][node]
                                ## Doing interpolation 
                                XSdifc = getXSlist(xsdict,fa,'difc',str(node))
                                XSnufis = getXSlist(xsdict,fa,'nufission',str(node))
                                XSabs = getXSlist(xsdict,fa,'absorption',str(node))
                                XSrem = getXSlist(xsdict,fa,'removal',str(node))
                                XS1 = [float(xs[0]) for xs in XSdifc]
                                XS2 = [float(xs[1]) for xs in XSdifc]
                                XS3 = [float(xs[0]) for xs in XSnufis]
                                XS4 = [float(xs[1]) for xs in XSnufis]
                                XS5 = [float(xs[0]) for xs in XSabs]
                                XS6 = [float(xs[1]) for xs in XSabs]
                                XS7 = [float(xs[0]) for xs in XSrem]
                                tempbranchXS_1[loc][node]=np.interp(buval, fabulist, XS1)
                                tempbranchXS_2[loc][node]=np.interp(buval, fabulist, XS2)
                                tempbranchXS_3[loc][node]=np.interp(buval, fabulist, XS3)
                                tempbranchXS_4[loc][node]=np.interp(buval, fabulist, XS4)
                                tempbranchXS_5[loc][node]=np.interp(buval, fabulist, XS5)
                                tempbranchXS_6[loc][node]=np.interp(buval, fabulist, XS6)
                                tempbranchXS_7[loc][node]=np.interp(buval, fabulist, XS7)
                                tempoutput[loc][node]=rpf3d[bu][idx][node]
                            idx = idx + 1
                    list_outputcore.append(coreparam[bu] )
                    list_outputRPF.append(tempoutput)        
                    list_branchXS_1.append(tempbranchXS_1)        
                    list_branchXS_2.append(tempbranchXS_2)        
                    list_branchXS_3.append(tempbranchXS_3)        
                    list_branchXS_4.append(tempbranchXS_4)        
                    list_branchXS_5.append(tempbranchXS_5)        
                    list_branchXS_6.append(tempbranchXS_6)        
                    list_branchXS_7.append(tempbranchXS_7)        
    
    if not list_branchXS_1:
        ## empty data
        raise ValueError("Failed to extract data due to empty list -- Check output PARCS.")

    branchXS_1 = np.vstack((branchXS_1,list_branchXS_1))
    branchXS_2 = np.vstack((branchXS_2,list_branchXS_2))
    branchXS_3 = np.vstack((branchXS_3,list_branchXS_3))
    branchXS_4 = np.vstack((branchXS_4,list_branchXS_4))
    branchXS_5 = np.vstack((branchXS_5,list_branchXS_5))
    branchXS_6 = np.vstack((branchXS_6,list_branchXS_6))
    branchXS_7 = np.vstack((branchXS_7,list_branchXS_7))
    outputRPF  = np.vstack((outputRPF ,list_outputRPF))
    outputcore = np.vstack((outputcore,list_outputcore))
    ## save to a temporary directory for retraining later 
    
    np.save(os.path.join(retrain_path, 'outputRFP.npy'),outputRPF)
    np.save(os.path.join(retrain_path, 'outputcore.npy'),outputcore)
    np.save(os.path.join(retrain_path, 'branchXS_1.npy'),branchXS_1)
    np.save(os.path.join(retrain_path, 'branchXS_2.npy'),branchXS_2)
    np.save(os.path.join(retrain_path, 'branchXS_3.npy'),branchXS_3)
    np.save(os.path.join(retrain_path, 'branchXS_4.npy'),branchXS_4)
    np.save(os.path.join(retrain_path, 'branchXS_5.npy'),branchXS_5)
    np.save(os.path.join(retrain_path, 'branchXS_6.npy'),branchXS_6)
    np.save(os.path.join(retrain_path, 'branchXS_7.npy'),branchXS_7)
    print('retrain data output shape',outputRPF.shape)
    # del xsdict 
    # gc.collect()
    # print(outputcore.shape)
    return



def getdata(bumap, xsdict, trunks, scalerparam, coremap,fabulist):
    '''
    Take the burnup map and FA XS data to generate branch data for surrogate model
    bumap : input (1,81,16)
    xsdict: input dictionary of FA XS data
    scaler: input for scaling
    coremap:input core FA layout [9x9]
    iterFAdata: output 7 branch and 1 trunk data 
    '''
    _,nfa,nz = bumap.shape ## _, 81, 16
    ##initializing
    tempbranchXS_1 = np.zeros((nfa,nz))
    tempbranchXS_2 = np.zeros((nfa,nz))
    tempbranchXS_3 = np.zeros((nfa,nz))
    tempbranchXS_4 = np.zeros((nfa,nz))
    tempbranchXS_5 = np.zeros((nfa,nz))
    tempbranchXS_6 = np.zeros((nfa,nz))
    tempbranchXS_7 = np.zeros((nfa,nz))
    for loc,fa in enumerate(coremap):
        if fa =='10':
            tempbranchXS_1[loc][:]=[-1]*nz
            tempbranchXS_2[loc][:]=[-1]*nz
            tempbranchXS_3[loc][:]=[-1]*nz
            tempbranchXS_4[loc][:]=[-1]*nz
            tempbranchXS_5[loc][:]=[-1]*nz
            tempbranchXS_6[loc][:]=[-1]*nz
            tempbranchXS_7[loc][:]=[-1]*nz
        elif fa != '00':
            for node in range(nz):
                buval = bumap[0,loc,node]
                ## Doing interpolation 
                XSdifc  = np.array(getXSlist(xsdict,fa,'difc',str(node)))
                XSnufis = np.array(getXSlist(xsdict,fa,'nufission',str(node)))
                XSabs   = np.array(getXSlist(xsdict,fa,'absorption',str(node)))
                XSrem   = np.array(getXSlist(xsdict,fa,'removal',str(node)))
                XSdifc_0, XSdifc_1 = XSdifc[:, 0], XSdifc[:, 1]
                XSnufis_0, XSnufis_1 = XSnufis[:, 0], XSnufis[:, 1]
                XSabs_0, XSabs_1 = XSabs[:, 0], XSabs[:, 1]
                XSrem_0 = XSrem[:, 0]
            
                tempbranchXS_1[loc, node] = np.interp(buval, fabulist, XSdifc_0)
                tempbranchXS_2[loc, node] = np.interp(buval, fabulist, XSdifc_1)
                tempbranchXS_3[loc, node] = np.interp(buval, fabulist, XSnufis_0)
                tempbranchXS_4[loc, node] = np.interp(buval, fabulist, XSnufis_1)
                tempbranchXS_5[loc, node] = np.interp(buval, fabulist, XSabs_0)
                tempbranchXS_6[loc, node] = np.interp(buval, fabulist, XSabs_1)
                tempbranchXS_7[loc, node] = np.interp(buval, fabulist, XSrem_0)

    ## scaling 
    tempbranchXS_1 = np.array(tempbranchXS_1).reshape(1,nfa,nz)
    tempbranchXS_2 = np.array(tempbranchXS_2).reshape(1,nfa,nz)
    tempbranchXS_3 = np.array(tempbranchXS_3).reshape(1,nfa,nz)
    tempbranchXS_4 = np.array(tempbranchXS_4).reshape(1,nfa,nz)
    tempbranchXS_5 = np.array(tempbranchXS_5).reshape(1,nfa,nz)
    tempbranchXS_6 = np.array(tempbranchXS_6).reshape(1,nfa,nz)
    tempbranchXS_7 = np.array(tempbranchXS_7).reshape(1,nfa,nz)
    X1_news  = minmaxscale(tempbranchXS_1,scalerparam['b1max'], scalerparam['b1min'])
    X2_news  = minmaxscale(tempbranchXS_2,scalerparam['b2max'], scalerparam['b2min'])
    X3_news  = minmaxscale(tempbranchXS_3,scalerparam['b3max'], scalerparam['b3min'])
    X4_news  = minmaxscale(tempbranchXS_4,scalerparam['b4max'], scalerparam['b4min'])
    X5_news  = minmaxscale(tempbranchXS_5,scalerparam['b5max'], scalerparam['b5min'])
    X6_news  = minmaxscale(tempbranchXS_6,scalerparam['b6max'], scalerparam['b6min'])
    X7_news  = minmaxscale(tempbranchXS_7,scalerparam['b7max'], scalerparam['b7min'])
    iterFAdata =(X1_news,
                 X2_news,
                 X3_news, 
                 X4_news, 
                 X5_news, 
                 X6_news, 
                 X7_news, trunks)
    return iterFAdata

## new version
def getXSlist(xsdict, fatype, xstype, axialnode):
    return [xsdict[fatype][key][xstype][axialnode] for key in xsdict[fatype]]

## test vectorization for interpolation

def getdata_optimized(bumap, xsdict, trunks, scalerparam, coremap, fabulist):
    '''
    Optimized version of the interpolation process
    '''
    _, nfa, nz = bumap.shape  # _, 81, 16
    
    # Pre-allocate all arrays at once
    tempbranchXS = np.zeros((7, nfa, nz))
    
    # Convert fabulist to numpy array once
    fabulist_np = np.array(fabulist)
    
    # Create mask for valid FAs (not '10' or '00')
    valid_mask = np.array([(fa != '10' and fa != '00') for fa in coremap])
    invalid_mask = np.array([fa == '10' for fa in coremap])
    
    # Set invalid locations to -1 all at once
    tempbranchXS[:, invalid_mask, :] = 0.0
    
    # Get unique FA types to reduce redundant dictionary lookups
    unique_fas = list(set([fa for fa in coremap if fa != '10' and fa != '00']))
    
    # Pre-compute XS data for all unique FA types
    xs_cache = {}
    for fa in unique_fas:
        xs_cache[fa] = {
            'difc': getXSlist_vectorized(xsdict, fa, 'difc'),
            'nufission': getXSlist_vectorized(xsdict, fa, 'nufission'), 
            'absorption': getXSlist_vectorized(xsdict, fa, 'absorption'),
            'removal': getXSlist_vectorized(xsdict, fa, 'removal')
        }
    
    # Process each unique FA type
    for fa in unique_fas:
        # Find all locations with this FA type
        fa_locations = np.where(np.array(coremap) == fa)[0]
        
        if len(fa_locations) == 0:
            continue
            
        # Get burnup values for all locations of this FA type
        buvals = bumap[0, fa_locations, :]  # shape: (n_locations, nz)
        
        # Get XS data for this FA type
        bukeys = list(xs_cache[fa]['difc'].keys())
        
        # Convert XS data to arrays for vectorized operations
        xs_arrays = {}
        for xs_type in ['difc', 'nufission', 'absorption', 'removal']:
            xs_data = xs_cache[fa][xs_type]
            # Stack all burnup points for all axial nodes
            xs_arrays[xs_type] = np.array([[xs_data[key][str(node)] for key in bukeys] 
                                          for node in range(nz)])  # shape: (nz, n_burnup_points, n_groups)
        
        # Vectorized interpolation for all locations and nodes at once
        for i, loc in enumerate(fa_locations):
            buval_loc = buvals[i, :]  # burnup for this location, all nodes
            
            for node in range(nz):
                buval = buval_loc[node]
                
                # Get XS values for this axial node
                difc_vals = xs_arrays['difc'][node]    # shape: (n_burnup_points, 2)
                nufis_vals = xs_arrays['nufission'][node]  # shape: (n_burnup_points, 2)  
                abs_vals = xs_arrays['absorption'][node]   # shape: (n_burnup_points, 2)
                rem_vals = xs_arrays['removal'][node]      # shape: (n_burnup_points, 1)
                
                # Vectorized interpolation
                tempbranchXS[0, loc, node] = np.interp(buval, fabulist_np, difc_vals[:, 0])
                tempbranchXS[1, loc, node] = np.interp(buval, fabulist_np, difc_vals[:, 1])
                tempbranchXS[2, loc, node] = np.interp(buval, fabulist_np, nufis_vals[:, 0])
                tempbranchXS[3, loc, node] = np.interp(buval, fabulist_np, nufis_vals[:, 1])
                tempbranchXS[4, loc, node] = np.interp(buval, fabulist_np, abs_vals[:, 0])
                tempbranchXS[5, loc, node] = np.interp(buval, fabulist_np, abs_vals[:, 1])
                tempbranchXS[6, loc, node] = np.interp(buval, fabulist_np, rem_vals[:, 0])
    
    # Reshape and scale - vectorized operations
    tempbranchXS = tempbranchXS.reshape(7, 1, nfa, nz)
    
    # Vectorized scaling
    scaled_results = []
    scale_keys = ['b1max', 'b2max', 'b3max', 'b4max', 'b5max', 'b6max', 'b7max']
    scale_mins = ['b1min', 'b2min', 'b3min', 'b4min', 'b5min', 'b6min', 'b7min']
    
    for i in range(7):
        scaled = minmaxscale(tempbranchXS[i], scalerparam[scale_keys[i]], scalerparam[scale_mins[i]])
        scaled_results.append(scaled)
    # Return as tuple
    return scaled_results 

def depletion_power(model, modelinput, scaler,idx1, idx2):
    '''
    Predict the new powdistribution based on new burnup distribution 
    input: bustep: bu step on the burnup distribution (1)
    input: model: DeepOnet model 
    input: modelinput: (7FA data, trunk) to generate input for model 
    input: scaler for scaling input output of model 
    ouput: new power distribution (1,81,16)
    '''
    # t0 = TT.time()
    newpow = model.predict(modelinput)
    # print('pred time, ',TT.time() -t0)
    tmp_pred = newpow.reshape(-1,81, 16)
    ### assign 0 to reflector and void region 
    tmp_pred[:,8, :] = 0
    tmp_pred[:,17, :] = 0
    tmp_pred[:,26, :] = 0
    tmp_pred[:,35, :] = 0
    tmp_pred[:, 43:45, :] = 0
    tmp_pred[:, 52:54, :] = 0
    tmp_pred[:, 60:63, :] = 0
    tmp_pred[:, 67:73, :] = 0
    tmp_pred[:, 73:, :] = 0
    # tmp_pred[abs(tmp_pred)<1e-4] = 0
    tmp_pred = normalize(tmp_pred,idx2, idx1,total_height, faaxial)
    # y_predr = minmaxReverse(tmp_pred.reshape(tmp_pred.shape[0],newpow.shape[-1]),scaler['omax'],scaler['omin'])
    y_predr = tmp_pred.reshape(newpow.shape)
    return y_predr

def depletion_power_157(model, modelinput, scaler,idx1, idx2):
    '''
    Predict the new powdistribution based on new burnup distribution for 157 core
    input: bustep: bu step on the burnup distribution (1)
    input: model: DeepOnet model 
    input: modelinput: (7FA data, trunk) to generate input for model 
    input: scaler for scaling input output of model 
    ouput: new power distribution (1,81,16)
    '''
    # t0 = TT.time()
    newpow = model.predict(modelinput)
    # print('pred time, ',TT.time() -t0)
    tmp_pred = newpow.reshape(-1,81, 16)
    ### assign 0 to reflector and void region 
    idx = np.r_[8, 17, 25, 26, 34, 35, 42:45, 50:54, 58:63, 65:tmp_pred.shape[1]]

    tmp_pred[:, idx, :] = 0
    # tmp_pred[abs(tmp_pred)<1e-4] = 0
    tmp_pred = normalize(tmp_pred,idx2, idx1,total_height, faaxial,157)
    # y_predr = minmaxReverse(tmp_pred.reshape(tmp_pred.shape[0],newpow.shape[-1]),scaler['omax'],scaler['omin'])
    y_predr = tmp_pred.reshape(newpow.shape)
    return y_predr

def depletion_boroncycle(model, modelinput, scaler):
    '''
    Predict the new boron concentration based on new burnup distribution 
    input: bustep: bu step on the burnup distribution (1)
    input: model: DeepOnet model 
    input: modelinput: (7FA data, trunk) to generate input for model 
    input: scaler for scaling input output of model 
    ouput: boron. cycle (1,2)
    '''
    # t0 = TT.time()
    temp = model.predict(modelinput)
    # print('predict boron takes', TT.time() -t0)
    y_predr = minmaxReverse(temp,scaler['omax'],scaler['omin'])
    # bor = y_predr[:,0]
    # cyc = y_predr[:,1]
    return y_predr[:,0], y_predr[:,1]

def getXSlist_vectorized(xsdict, fatype, xstype):
    """Vectorized version that gets all axial nodes at once"""
    bukeys = list(xsdict[fatype].keys())
    # Get all axial nodes for this FA type and XS type at once
    xsvals = {}
    for key in bukeys:
        xsvals[key] = xsdict[fatype][key][xstype]
    return xsvals

def interp_nd(fabulist, pp, buval):
    """
    Interpolates pp values at query points buval based on fabulist.

    Parameters
    ----------
    fabulist : array-like, shape (N,)
        Monotonic 1D array of x-values.
    pp : array-like, shape (N, F)
        2D array of y-values, where F is number of features.
    buval : float or array-like
        Query point(s) to interpolate.

    Returns
    -------
    pp_interp : ndarray
        Interpolated values.
        Shape (F,) if buval is scalar,
        Shape (M, F) if buval has M elements.
    """
    fabulist = np.asarray(fabulist)
    pp = np.asarray(pp)
    buval = np.atleast_1d(buval)  # ensure array

    # Find interval indices
    idx = np.searchsorted(fabulist, buval) - 1
    idx = np.clip(idx, 0, len(fabulist) - 2)

    x0, x1 = fabulist[idx], fabulist[idx + 1]    # (M,)
    y0, y1 = pp[idx, :], pp[idx + 1, :]          # (M, F)

    # Linear interpolation
    t = (buval - x0) / (x1 - x0)                 # (M,)
    pp_interp = y0 + (y1 - y0) * t[:, None]      # (M, F)

    # If input was scalar → return 1D (F,)
    if np.isscalar(buval) or np.ndim(buval) == 0:
        return pp_interp[0]
    return pp_interp

def getfqFd(bumap, powmap, coremap, fabulist, xsdict, faaxial, total_height):
    '''
    Get Fq, Fdel at one burn up step 
    bumap: input, core burn up distribution
    powmap: input, core relative power factor distribution
    coremap: input, core loading pattern 
    fabulist: input, list of FA depletion step for interpolation
    xsdict: input, dictionary input 
    '''
    _,nfa,nz = bumap.shape ## _, 81, 16
    powmap = powmap.reshape(bumap.shape)
    fqcheck =0
    fdelcheck =0
    for loc,fa in enumerate(coremap):
        temp2=0
        if fa not in ['10','00']:
            for node in range(nz):
                buval = bumap[0,loc,node]
                nodepower = powmap[0,loc,node]
                pp = getXSlist(xsdict,fa,'ppowers',str(node))
                pindata = interp_nd(fabulist, pp, buval).flatten()*nodepower
                temp2 += pindata*faaxial[node]
                fqcheck = max(fqcheck,np.max(pindata))
            fdelcheck = max(fdelcheck,np.max(temp2/total_height))

    return fdelcheck, fqcheck

def getfqFd_pinrecontruct(bumap, powmap, coremap, fabulist, xsdict, faaxial, total_height):
    '''
    Get Fq, Fdel at one burn up step 
    bumap: input, core burn up distribution
    powmap: input, core relative power factor distribution
    coremap: input, core loading pattern 
    fabulist: input, list of FA depletion step for interpolation
    xsdict: input, dictionary input 
    '''
    _,nfa,nz = bumap.shape ## _, 81, 16
    powmap = powmap.reshape(bumap.shape)
    fqcheck =0
    fdelcheck =0
    #pinpower = np.zeros((nz,153,153))
    ttt = np.zeros((nz,nfa,289)) 
    for loc,fa in enumerate(coremap):
        if fa not in ['10','00']:
            for node in range(nz):
                buval = bumap[0,loc,node]
                nodepower = powmap[0,loc,node]
                pp = getXSlist(xsdict,fa,'ppowers',str(node))
                pindata = interp_nd(fabulist, pp, buval).flatten()*nodepower # 17*17
                ttt[node,loc,:]=pindata
    pinpower = (
    ttt.reshape(16, 9, 9, 17, 17)      # break into blocks
       .transpose(0, 1, 3, 2, 4)       # reorder so dims match (9*17, 9*17)
       .reshape(16, 9*17, 9*17)        # final shape
       ).reshape(-1,153,153,1)
    
    return pinpower


def getfqFd_vectorized(pimodel,bumap, powmap, coremap, fabulist, xsdict, faaxial, total_height):
    """
    Fully vectorized version to calculate Fq and Fdel, eliminating all Python loops.
    """
    _, nfa, nz = bumap.shape
    powmap = powmap.reshape(bumap.shape)
    coremap_np = np.array(coremap)
    fabulist_np = np.sort(np.array(fabulist, dtype=np.float64))
    n_burnup_points = len(fabulist_np)

    # --- Pre-computation Step ---
    # 1. Get all unique fuel assembly types that need processing
    unique_fas = list(set(fa for fa in coremap if fa not in ['10', '00']))
    
    # 2. Create a large array to hold all pin power data
    # Dimensions: (num_unique_fas, num_burnup_points, nz, num_pins_per_fa)
    # Assuming a fixed number of pins, e.g., 289 (17x17)
    num_pins = 289 # Example for a 17x17 assembly
    all_pp_data = np.zeros((len(unique_fas), n_burnup_points, nz, num_pins))
    
    # 3. Create a mapping from FA type to an index in the large array
    fa_to_idx = {fa: i for i, fa in enumerate(unique_fas)}

    # 4. Populate the large array with data from xsdict
    for fa_idx, fa in enumerate(unique_fas):
        for bu_idx, bu in enumerate(fabulist_np):
            for node in range(nz):
                # Assuming getXSlist returns a dictionary-like object that needs to be converted
                # to an array. This part might need adjustment based on the exact structure.
                pp_data = xsdict[fa]['BU='+str(bu)]['ppowers'][str(node)]
                # Ensure pp_data is a flat numpy array of the correct size
                all_pp_data[fa_idx, bu_idx, node, :] = pp_data.flatten()

    # --- Vectorized Calculation Step ---
    # 1. Create masks for valid locations and map them to their pre-computed data index
    valid_mask = np.isin(coremap_np, unique_fas)
    fa_indices_for_core = np.array([fa_to_idx.get(fa, -1) for fa in coremap_np])
    
    # 2. Get all burnup and power values for valid locations
    buvals = bumap[0, valid_mask, :]
    powvals = powmap[0, valid_mask, :]
    
    # 3. Perform vectorized interpolation
    buvals_clipped = np.clip(buvals, fabulist_np[0], fabulist_np[-1])
    indices = np.searchsorted(fabulist_np, buvals_clipped)
    np.clip(indices, 1, n_burnup_points - 1, out=indices)
    
    x1 = fabulist_np[indices - 1]
    x2 = fabulist_np[indices]
    frac = (buvals_clipped - x1) / (x2 - x1 + 1e-9)
    
    # 4. Use advanced indexing to get y1 and y2 pin power values
    # Get the correct FA data index for each valid location
    fa_indices_at_valid_locs = fa_indices_for_core[valid_mask]
    
    # Prepare indices for broadcasting and advanced indexing
    loc_indices = np.arange(len(fa_indices_at_valid_locs))
    node_indices = np.arange(nz)

    # Index into the large pre-computed array
    y1 = all_pp_data[fa_indices_at_valid_locs[:, np.newaxis], indices - 1, node_indices, :]
    y2 = all_pp_data[fa_indices_at_valid_locs[:, np.newaxis], indices, node_indices, :]
    
    # 5. Calculate interpolated pin data for all valid locations, nodes, and pins at once
    interp_pin_data = y1 + frac[..., np.newaxis] * (y2 - y1)
    
    # 6. Apply the node power factor
    pindata = interp_pin_data * powvals[..., np.newaxis]
    
    # 7. Calculate Fq (max of all pin data)
    fqcheck = np.max(pindata) if pindata.size > 0 else 0
    
    # 8. Calculate Fdel
    # Sum pin data over the pin dimension (axis=2) and multiply by axial height factor
    faaxial_bcast = np.array(faaxial).reshape(1, nz)
    temp2 = np.sum(pindata * faaxial_bcast[:, :, np.newaxis], axis=1)

    fdelcheck = np.max(temp2 / total_height) if temp2.size > 0 else 0
    
    return fdelcheck, fqcheck


# 
def interpolatecycle(x, X, Y):
    """
    Parameters:
        x (float): ending boron ~10 ppm.
        X (list or array): Known boron concentration over depletion.
        Y (list or array): Coresponding cycle.

    Returns:
        float: Interpolated cycle length.
    """
    X=list(X)
    Y=list(Y)
    idx = -1
    ## test 
    # idx = next((i for i in range(len(X)-1) if abs(X[i+1][0] - X[i][0]) < 10 and X[i][0]<500) , -1)
    # idx=idx-4 ## take 3 step back 


    # idx = -5
    # for j in range(len(X)):
    #     if X[j][0]<x:
    #         idx = j-5 ## take one idx back
    #         break 
    ## only take 2 last element for cycle length calculation 
    a = (X[idx-1][0]-X[idx][0])/(Y[idx-1][0]-Y[idx][0])
    b = X[idx-1][0] - a*Y[idx-1][0]
    if Y[0][0]>0:
        y = (x-b)/a - Y[0][0] ## normalize to make the first value is 0 
    else:
        y = (x-b)/a
    return y



def precompute_xs_data(xsdict, unique_fas, fabulist_np):
    """
    Pre-computes and stacks all cross-section data for all unique fuel types upfront.
    This moves the expensive data preparation out of the main processing loop.
    """
    # print("Pre-computing and stacking all cross-section data...")
    precomputed_cache = {}
    bukeys = ['BU='+str(bu) for bu in fabulist_np]
    
    for fa in unique_fas:
        fa_cache = {}
        
        # --- FIX for 'removal' data structure ---
        # Get the raw removal data which may have multiple columns
        removal_data_raw = getXSlist_vectorized(xsdict, fa, 'removal')
        
        # --- FIX for TypeError: unhashable type: 'slice' ---
        # Use a nested comprehension to iterate into the inner dictionary of nodes
        # before applying the index to the innermost numpy array.
        removal_data_corrected = {
            bu_key: {
                node_key: node_array[0] for node_key, node_array in node_dict.items()
            }
            for bu_key, node_dict in removal_data_raw.items()
        }

        xs_data_for_fa = {
            'difc': getXSlist_vectorized(xsdict, fa, 'difc'),
            'nufission': getXSlist_vectorized(xsdict, fa, 'nufission'),
            'absorption': getXSlist_vectorized(xsdict, fa, 'absorption'),
            'removal': removal_data_corrected
        }
        
        for xs_type in ['difc', 'nufission', 'absorption', 'removal']:
            # xs_data_dict is a dictionary of {burnup_key: {node_key: xs_array}}
            xs_data_dict = xs_data_for_fa[xs_type]
            
            # --- FIX for Data Structure Mismatch ---
            # Reconstruct the expected array structure from the dictionary of nodes.
            # The goal is to create a list of arrays, where each array represents
            # a burnup step and has the shape (nz, n_groups).
            
            reconstructed_xs_list = []
            for bu_key in bukeys:
                # node_dict is the dictionary for a single burnup step, e.g., {'0': array(2,), '1': array(2,), ...}
                node_dict = xs_data_dict[bu_key]
                
                # Sort keys numerically to ensure correct axial node order ('0', '1', ..., '15')
                #sorted_node_keys = sorted(node_dict.keys(), key=int)
                
                # Create a list of arrays for each node
                # node_arrays = [node_dict[key] for key in sorted_node_keys]
                node_arrays = [node_dict[key] for key in node_dict.keys()]
                
                # Stack the node arrays to form a single (nz, n_groups) array for this burnup step
                stacked_array_for_bu_step = np.stack(node_arrays)
                reconstructed_xs_list.append(stacked_array_for_bu_step)

            # Now `reconstructed_xs_list` is the list of arrays in the correct format.
            sorted_xs_list = reconstructed_xs_list

            # The original logic for handling 1D vs 2D arrays can now proceed correctly.
            processed_list = []
            for arr in sorted_xs_list:
                if arr.ndim == 1:
                    # Convert 1D array of shape (nz,) to a 2D column vector of shape (nz, 1)
                    processed_list.append(arr[:, np.newaxis])
                else:
                    processed_list.append(arr)
            
            # The expensive stacking operation is now done only once per FA type here.
            fa_cache[xs_type] = np.stack(processed_list, axis=1)
            
        precomputed_cache[fa] = fa_cache
        
    return precomputed_cache

def getdata_fully_vectorized(bumap, xsdict, coremap, fabulist, scalerparam):
    """
    Fully vectorized version that pre-computes data structures to minimize
    work inside the main loop.
    Includes an optional validation step against np.interp.
    """
    _, nfa, nz = bumap.shape  # _, 81, 16
    tempbranchXS = np.zeros((7, nfa, nz))
    # fabulist_np = np.sort(np.array(fabulist, dtype=np.float64))
    fabulist_np = fabulist
    n_burnup_points = len(fabulist_np)
    
    coremap_np = np.array(coremap)
    invalid_mask = (coremap_np == '10')
    tempbranchXS[:, invalid_mask, :] = 0.0
    
    unique_fas = list(set([fa for fa in coremap if fa != '10' and fa != '00']))
    
    # --- MAJOR OPTIMIZATION: Pre-compute all data structures ---
    xs_cache = precompute_xs_data(xsdict, unique_fas, fabulist_np)
    
    # Process each unique FA type
    for fa in unique_fas:
        fa_locations = np.where(coremap_np == fa)[0]
        if len(fa_locations) == 0:
            continue
            
        buvals = bumap[0, fa_locations, :]
        
        # --- OPTIMIZED: Simple dictionary lookup for pre-stacked data ---
        xs_arrays = xs_cache[fa]

        # --- FULLY VECTORIZED INTERPOLATION (CORRECTED) ---
        # Clip buvals to the range of fabulist_np to mimic np.interp's behavior 
        # for out-of-bounds values (clamping instead of extrapolating).
        buvals_clipped = np.clip(buvals, fabulist_np[0], fabulist_np[-1])
        
        indices = np.searchsorted(fabulist_np, buvals_clipped)
        np.clip(indices, 1, n_burnup_points - 1, out=indices)
        
        x1 = fabulist_np[indices - 1]
        x2 = fabulist_np[indices]
        
        # Use the clipped values here as well to ensure frac is always between 0 and 1.
        frac = (buvals_clipped - x1) / (x2 - x1 + 1e-9)
        frac_bcast = frac[..., np.newaxis]

        node_indices = np.arange(nz).reshape(1, nz)
        
        y1_difc = xs_arrays['difc'][node_indices, indices - 1]
        y2_difc = xs_arrays['difc'][node_indices, indices]
        
        y1_nufis = xs_arrays['nufission'][node_indices, indices - 1]
        y2_nufis = xs_arrays['nufission'][node_indices, indices]
        
        y1_abs = xs_arrays['absorption'][node_indices, indices - 1]
        y2_abs = xs_arrays['absorption'][node_indices, indices]
        
        y1_rem = xs_arrays['removal'][node_indices, indices - 1]
        y2_rem = xs_arrays['removal'][node_indices, indices]

        interp_difc = y1_difc + frac_bcast * (y2_difc - y1_difc)
        interp_nufis = y1_nufis + frac_bcast * (y2_nufis - y1_nufis)
        interp_abs = y1_abs + frac_bcast * (y2_abs - y1_abs)
        interp_rem = y1_rem + frac_bcast * (y2_rem - y1_rem)

        tempbranchXS[0:2, fa_locations, :] = interp_difc.transpose(2, 0, 1)
        tempbranchXS[2:4, fa_locations, :] = interp_nufis.transpose(2, 0, 1)
        tempbranchXS[4:6, fa_locations, :] = interp_abs.transpose(2, 0, 1)
        # --- FIX for ValueError: More robustly select the data instead of squeezing ---
        tempbranchXS[6, fa_locations, :] = interp_rem[..., 0]

        # # --- VALIDATION STEP (runs only if validate_results is True) ---
        # if validate_results:
        #     print(f"\n--- Validating results for FA type: {fa} ---")
        #     check_results = np.zeros_like(tempbranchXS[:, fa_locations, :])
        #     for i, loc in enumerate(fa_locations):
        #         for node in range(nz):
        #             buval = bumap[0, loc, node]
                    
        #             check_results[0, i, node] = np.interp(buval, fabulist_np, xs_arrays['difc'][node, :, 0])
        #             check_results[1, i, node] = np.interp(buval, fabulist_np, xs_arrays['difc'][node, :, 1])
        #             check_results[2, i, node] = np.interp(buval, fabulist_np, xs_arrays['nufission'][node, :, 0])
        #             check_results[3, i, node] = np.interp(buval, fabulist_np, xs_arrays['nufission'][node, :, 1])
        #             check_results[4, i, node] = np.interp(buval, fabulist_np, xs_arrays['absorption'][node, :, 0])
        #             check_results[5, i, node] = np.interp(buval, fabulist_np, xs_arrays['absorption'][node, :, 1])
        #             check_results[6, i, node] = np.interp(buval, fabulist_np, xs_arrays['removal'][node, :, 0])
            
        #     are_close = np.allclose(tempbranchXS[:, fa_locations, :], check_results)
        #     print(f"Vectorized results match np.interp results: {are_close}")
        #     if not are_close:
        #          print("Mismatch found!")
    
    # --- Loop-based Scaling using the optimized minmaxscale function ---
    tempbranchXS = tempbranchXS.reshape(7, 1, nfa, nz)
    scaled_results = []
    scale_keys = ['b1max', 'b2max', 'b3max', 'b4max', 'b5max', 'b6max', 'b7max']
    scale_mins = ['b1min', 'b2min', 'b3min', 'b4min', 'b5min', 'b6min', 'b7min']
    
    for i in range(7):
        data_slice = tempbranchXS[i]
        scaled = minmaxscale(data_slice, scalerparam[scale_keys[i]], scalerparam[scale_mins[i]])
        # Reshape the result to the desired (1, nfa, nz) shape
        scaled_results.append(scaled)
    
    return scaled_results


## core map
xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/updateXS.pkl','rb'))
## type --> BU value --> xs type --> axial node --> value

fabulist = [0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 
            17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,\
            65.0, 70.0, 75.0, 80.0]
fabulist = np.array(fabulist)
# y1 = np.interp(x1, x, y)
faaxial = np.array([15.24, 10.16, 5.08,
                  30.48, 30.48, 30.48, 30.48, 30.48,
                  30.48, 30.48, 30.48, 30.48, 30.48,
                  5.08, 10.16, 15.24])
total_height = np.sum(faaxial)

template_lines = [
    '      {fa5}  {fa1}  {fa1}  {fa2}  {fa4}  {fa5}  {fa1}  {fa3}  10 \n',
    '      {fa1}  {fa2}  {fa3}  {fa4}  {fa5}  {fa2}  {fa4}  {fa3}  10 \n',
    '      {fa3}  {fa5}  {fa2}  {fa4}  {fa1}  {fa4}  {fa1}  {fa4}  10 \n',
    '      {fa2}  {fa3}  {fa1}  {fa2}  {fa4}  {fa5}  {fa3}  {fa2}  10 \n',
    '      {fa1}  {fa5}  {fa2}  {fa1}  {fa4}  {fa5}  {fa2}  10   10 \n',
    '      {fa3}  {fa2}  {fa5}  {fa4}  {fa3}  {fa1}  {fa4}  10   00 \n',
    '      {fa1}  {fa3}  {fa4}  {fa2}  {fa1}  {fa5}  10   10   00 \n',
    '      {fa3}  {fa3}  {fa5}  {fa1}  10   10   10   00   00 \n',
    '      10   10   10   10   10   00   00   00   00 \n'
]
# Index weight mapping
idx11 = np.array([0])
idx22 = np.array([1, 2, 3, 4, 5, 6, 7, 9, 18, 27, 36, 45, 54, 63])

def compute_tempnorm(output,idx2, idx1,total_height, faaxial, count_per_sample):
    N, M, _ = output.shape

    # Default weight is 4
    weights = np.full(M, 4, dtype=np.float32)
    #counts = np.full(M, 0, dtype=np.float32)

    weights[idx1] = 1
    #counts[idx1] = 1

    weights[idx2] = 2
    #counts[idx2] = 2

    # Normalize axial heights
    axial_weights = faaxial / total_height  # shape (16,)

    # Compute per-sample weighted sum
    # shape: (N, M, 16) * (16,) → broadcast → (N, M, 16)
    weighted_output = output * axial_weights

    # Apply radial weights: multiply along axis=1 (M)
    # shape: (M,) → reshape to (1, M, 1)
    weighted_output *= weights.reshape(1, -1, 1)

    # Sum over M and 16 axes
    tolRFP = np.sum(weighted_output, axis=(1, 2))  # shape (N,)

    # # Compute total weight count for each sample
    # # If any value in output > 1e-5 and index not in idx1 or idx2, count += 4
    # mask_large = np.max(output, axis=2) > 1e-5  # shape (N, M)

    # is_other_idx = ~np.isin(np.arange(M), np.concatenate([idx1, idx2]))  # shape (M,)
    # extra_counts = mask_large[:, is_other_idx] * 4

    # count_per_sample = np.sum(counts) + np.sum(extra_counts, axis=1)  # shape (N,)
    # count_per_sample = 193
    # count_per_sample = 157
    # Final normalization
    tempnorm = (tolRFP / count_per_sample).reshape(-1, 1)  # shape (N, 1)
    return tempnorm

def normalize(output,idx2, idx1,total_height, faaxial,count):
    # tempnorm = compute_tempnorm(output,idx2, idx1,total_height, faaxial, 193)
    tempnorm = compute_tempnorm(output,idx2, idx1,total_height, faaxial,count)
    outnew = [output[i]*1.0/tempnorm[i] for i in range(len(output))]
    outnew = np.array(outnew).reshape(output.shape)
    return outnew

listFA = [461,462,501,502,526,566,586,250,280,320,400,567,587]
# combo_generator = itertools.product(listFA, repeat=5)
corebulist = [0.1, 0.4, 0.5,1.0,1.0]
for i in range(28):
    corebulist.append(1*1.0)
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
index=-1
initbumap = np.zeros((1,81,16))
def get_result(LPs, bustep=bustep, count=count, 
               cycle_length_pred=cycle_length_pred, 
               boronmax_pred = boronmax_pred,
               Fqmax=Fqmax, Fdmax= Fdmax, corebulist = corebulist, listFA=listFA, index=index, initbumap=initbumap) :
    test = " ".join(LPs).strip().split()
    ## initializing 
    # inputdata = getdata(initbumap, xsdict, trunks, scalerparam, test,fabulist)
    Fq_all = 0
    Fd_all = 0
    pinpower =[]
    t0 = TT.time()
    print('before predccccc')
    # inputdataold = getdata_optimized(initbumap, xsdict, trunks, scalerparam, test,fabulist)
    inputdata = getdata_fully_vectorized(initbumap, xsdict, test, fabulist, scalerparam)
    # tdata = TT.time()
    # print('initial data', tdata-st)
    print('after get data')
    # pow_0 = depletion_power(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22)
    pow_0 = depletion_power_157(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22)
    print('predictin')
    # Fd, Fq = getfqFd(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    # Fd, Fq = getfqFd_vectorized(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    pinpowerbu = getfqFd_pinrecontruct(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    pinpower.append(pinpowerbu)
    
    # tt = TT.time()
    # Fd, Fq = getfqFd(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    # print('Fqfd old time', TT.time()-tt)
    # stop
    # Fd_all =max(Fd_all,Fd)
    # Fq_all =max(Fq_all,Fq)
    bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, tuple(inputdata+[trunk_core]), scaler_core)
    # print('initial prediction', TT.time()-tdata)
    # bor_true0, cyc_true0 =bor_truth[bustep], cyc_truth [bustep]
    bustep+=1
    # t0 = TT.time()
    # print(f'Initial step take {(t0-st):.5f}')
    # pow_0 = model.predict(inputdata)
    # Initialize
    power_history = [pow_0]
    boron_his = [bor_pred0]
    cycle_his = [cyc_pred0]
    bu_his =[initbumap]
    complete_flag = False
    # boron_his_true = [bor_true0]
    # cycle_his_true = [cyc_true0]
    for bu in corebulist:
            converged = False
            # tobu = TT.time()
            for iteration in range(100):
                if iteration == 0:
                    delE = bu * pow_0
                    bumap = initbumap + delE.reshape((1, 81, 16))
                else:
                    delE = bu * (0.5 * pow_0 + 0.5 * pow_prev)
                    bumap = initbumap + delE.reshape((1, 81, 16))
                tdata1 = TT.time()
                # input_data_old = getdata_optimized(bumap, xsdict, trunks, scalerparam, test, fabulist)
                input_data = getdata_fully_vectorized(bumap, xsdict, test, fabulist, scalerparam)
                # tdata2 = TT.time()
                # print('dada cal', tdata2-tdata1)
                # pow_new = depletion_power(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22)
                pow_new = depletion_power_157(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22)
                # print('model cal', TT.time() -tdata2)
                # Check convergence
                if iteration > 0:
                    max_delta = np.max(np.abs(pow_new - pow_prev))
                    if max_delta < 1e-4:
                        # print('Time get new core pow and bumap ', TT.time()-tobu)
                        # tfq = TT.time()
                        # Fd, Fq = getfqFd(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        # Fd, Fq = getfqFd_vectorized(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        pinpowerbu= getfqFd_pinrecontruct(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        pinpower.append(pinpowerbu)
                        # print('time to get Fq fd', TT.time()-tfq)
                        # tbor = TT.time()
                        bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, 
                                                                    tuple(input_data+[trunk_core]), 
                                                                    scaler_core)
                        # print('Time to get boron ', TT.time() - tbor)
                        # print('Time to get all params ', TT.time() - tobu)
                        # bor_true0, cyc_true0 =bor_truth[bustep], cyc_truth [bustep]
                        # Fd_all =max(Fd_all,Fd)
                        # Fq_all =max(Fq_all,Fq)
                        # print(Fq, Fd)
                        # print(f"BU={bu}, iter={iteration} → Converged (tol={max_delta:.2e})")
                        # print(f"step takes {TT.time() - tobu:.5f}")
                        # Update for next burnup
                        initbumap = bumap
                        bu_his.append(initbumap)
                        power_history.append(pow_new)
                        boron_his.append(bor_pred0)
                        cycle_his.append(cyc_pred0)
                        # boron_his_true.append(bor_true0)
                        # cycle_his_true.append(cyc_true0)
                        pow_0 = pow_new
                        converged = True
                        bustep+=1
                        # if bor_pred0[0] < 20:
                        if bor_pred0[0] < 10:
                            complete_flag=True
                        break
        
                pow_prev = pow_new
        
            if not converged:
                print(f"BU={bu} → ❌ Unconverged after 100 iterations")
                raise ValueError("Solution did not converge.")
            if complete_flag:
                break
    print('Deeponet time', TT.time() -t0)
    ## get the Fdmax and Fqmax here 
    pinpower = np.array(pinpower).reshape(-1,153,153,1)
    pinpower_reconstruct = pinmodel.predict(pinpower, verbose=0) # verbose=0 suppresses Keras logs for speed
    print('pred pin time', TT.time()-t0)
    # 6. Apply Mask and Calculate Statistics
    # Create mask from the first node (assuming geometry is axially constant)
    # Using np.not_equal is slightly faster/cleaner than != 0
    keep_mask = np.not_equal(pinpower[0], 0) 
    
    # Broadcast multiply: (16,153,153,1) * (153,153,1)
    pinpower_reconstruct *= keep_mask
    # ## test plot
    # print("Generating and saving prediction plots...")
    # IMG_DIM =153
    # # Reshape for plotting
    # y_pred_corrected_img = pinpower_reconstruct.reshape(-1, IMG_DIM, IMG_DIM)
    # X_test_img = pinpower.reshape(-1, IMG_DIM, IMG_DIM)
    
    
    # # Plot 3 random examples with 4 columns
    # num_examples = 3
    # indices = np.random.choice(range(len(pinpower)), num_examples, replace=False)
    
    # fig, axes = plt.subplots(num_examples, 4, figsize=(20, 5 * num_examples))
    # fig.suptitle("DEEP U-Net Prediction (No Scaling) vs. Ground Truth", fontsize=16)
    
    # for i, idx in enumerate(indices):
    #     input_img = X_test_img[idx]
    #     true_img = y_pred_corrected_img[idx]
    #     pred_img = y_pred_corrected_img[idx]
        
    #     # Calculate Difference (Prediction - Truth)
    #     diff_img = pred_img - true_img
    
    #     # Determine consistent min/max for the main plots
    #     vmin = min(true_img.min(), pred_img.min())
    #     vmax = max(true_img.max(), pred_img.max())
        
    #     print(f'Sample {idx} - Max True: {true_img.max():.4f}, Max Pred: {pred_img.max():.4f}, Max old {input_img.max():.4f}')
    
    #     # 1. Plot Input
    #     ax = axes[i, 0]
    #     im0 = ax.imshow(input_img, cmap='viridis', origin='lower')
    #     ax.set_title(f"Input (Sample {idx})")
    #     fig.colorbar(im0, ax=ax, label='Pin Power')
    
    #     # 2. Plot Ground Truth
    #     ax = axes[i, 1]
    #     im1 = ax.imshow(true_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
    #     ax.set_title(f"True")
    #     fig.colorbar(im1, ax=ax, label='Pin Power')
    
    #     # 3. Plot Predicted
    #     ax = axes[i, 2]
    #     im2 = ax.imshow(pred_img, cmap='viridis', origin='lower')#, vmin=vmin, vmax=vmax)
    #     ax.set_title(f"Predicted (Corrected)")
    #     fig.colorbar(im2, ax=ax, label='Pin Power')
        
    #     # 4. Plot Difference (Pred - True)
    #     ax = axes[i, 3]
    #     # Center the colormap at 0 to show + (red) and - (blue) errors
    #     diff_limit = max(abs(diff_img.min()), abs(diff_img.max()))
    #     im3 = ax.imshow(diff_img, cmap='bwr', origin='lower', vmin=-diff_limit, vmax=diff_limit)
    #     ax.set_title(f"Difference (Pred - True)")
    #     fig.colorbar(im3, ax=ax, label='Error')
    
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.savefig('checkpre.png')
    # stop
    # Calculate Max Fq
    Fq_all = np.max(pinpower_reconstruct)
    print(pinpower_reconstruct.shape)
    print(Fq_all)
    tetetetetet = pinpower_reconstruct.reshape(37,16,153,153,1)
    print(np.max(tetetetetet, axis=(1, 2, 3, 4)))   # shape: (M,))
    print(np.unravel_index(np.argmax(tetetetetet), tetetetetet.shape))
    stop
    # Calculate Max Fdel (Fdelta-H equivalent?)
    # Pre-calculate axial shape factor
    ax_factor = (np.array(faaxial) / total_height).reshape(1, 16, 1, 1)
    
    # Multiply and find max
    Fd_all = np.max(np.sum(pinpower_reconstruct.reshape(-1,16,153,153) * ax_factor,axis = 1))
    Fdmax.append(Fd_all)  
    Fqmax.append(Fq_all) 
    boronmax_pred.append(np.max(boron_his))
    for kk in range(len(boron_his)):
        print(cycle_his[kk], boron_his[kk])
    
    
    # boronmax_true.append(np.max(boron_his_true))
    cycle_length_pred.append(interpolatecycle(10,boron_his,cycle_his))
    # cycle_length_true.append(interpolatecycle(10,boron_his_true,cycle_his_true))
    count+=1
    t2 = TT.time()
    # np.save('tempBUhist.npy',np.array(bu_his))
    print(f'one cycle take {(t2-t0):.5f}')
    cycle_his_adj = cycle_his-cycle_his[0]
    index = next(i for i, x in enumerate(boron_his) if x[0] < 50)
    print(Fd_all, Fq_all, np.max(boron_his), interpolatecycle(10,boron_his[:index],cycle_his_adj[:index]))
    #stop
    return Fd_all, Fq_all, np.max(boron_his), interpolatecycle(10,boron_his[:index],cycle_his_adj[:index])

# LPs = ['566  526  526  502  461  462  586  461  10 \n', 
#        '526  566  461  526  586  502  526  526  10 \n', 
#        '526  461  526  586  566  502  461  566  10 \n', 
#        '502  526  586  502  526  586  566  461  10 \n', 
#        '461  586  566  526  586  502  566  10   10 \n',
#        '462  502  502  586  502  586  586  10   00 \n',
#        '586  526  461  566  566  586  10   10   00 \n',
#        '461  526  566  461  10   10   10   00   00 \n',
#        '10   10   10   10   10   00   00   00   00 \n']

# LPs = " ".join(LPs).strip().split()
# for kk in range(1):
#     print(get_result(LPs))
# print((TT.time()-st)/100)