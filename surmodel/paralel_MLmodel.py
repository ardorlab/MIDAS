'''
Updated depletion module to calculate burnup in parallel.
Each worker process maintains its own instances of the three pretrained models.
'''

# import os
# import shutil
# import random
# import subprocess
# import itertools
# from pathlib import Path
# import pickle 
# import numpy as np
# import copy
# import time as TT
# import matplotlib.pyplot as plt
# import tensorflow as tf

# # Standard DeepXDE/Model imports
# from loadmodel_Nov25 import minmaxReverse, minmaxscale, scalerparam, trunks
# from loadmodelcore import trunks as trunk_core, scalerparam as scaler_core
# Note: We avoid importing 'model' or 'pinmodel' globally to prevent process deadlocks

# --- GLOBAL WORKER STATE ---
# This dictionary will store the models locally within each process memory
# worker_models = {}

import os
import importlib
from pathlib import Path
import sys
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["AUTOGRAPH_VERBOSITY"] = "0"
os.environ["TF_AUTOGRAPH_DISABLE"] = "1"
os.environ["TMPDIR"] = "/tmp"
import socket
import numpy as np
import tensorflow as tf
import deepxde as dde


def init_worker():
    """ 
    Initializes each worker process. This function loads the models 
    only once when the process starts.
    """
    
    global worker_data
    
    import gc
    if 'worker_data' in globals():
        del worker_data
    gc.collect()
    worker_data = {}
    current_file_path = Path(__file__).resolve()
    source_path = str(current_file_path.parent) 

    # 2. Add it to sys.path if it's not already there
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    import loadmodel_Nov25
    import loadmodelcore
    import loadmodelpin
    importlib.reload(loadmodel_Nov25)
    importlib.reload(loadmodelcore)
    importlib.reload(loadmodelpin)
    worker_data['model'] = loadmodel_Nov25.model
    worker_data['scalerparam'] = loadmodel_Nov25.scalerparam
    worker_data['trunks'] = loadmodel_Nov25.trunks
    worker_data['minmaxscale'] = loadmodel_Nov25.minmaxscale
    worker_data['minmaxReverse'] = loadmodel_Nov25.minmaxReverse
    worker_data['modelcore'] = loadmodelcore.model
    worker_data['scaler_core'] = loadmodelcore.scalerparam
    worker_data['trunk_core'] = loadmodelcore.trunks
    worker_data['pinmodel'] = loadmodelpin.best_model
    try:
      cpu_id = os.sched_getcpu()
    except AttributeError:
       cpu_id = "N/A"
    pid = os.getpid()
    hostname = socket.gethostname()
    slurm_task = os.environ.get("SLURM_PROCID", "N/A")
    slurm_localid = os.environ.get("SLURM_LOCALID", "N/A")
    
    print(
        f"[Node: {hostname}] "
        f"[PID: {pid}] "
        f"[CPU: {cpu_id}] "
        f"[SLURM_PROCID: {slurm_task}] "
        f"[SLURM_LOCALID: {slurm_localid}] "
        f"successfully loaded all 3 models."
    )
    #print(f"Worker {os.getpid()} successfully loaded all 3 models.")

def init_worker_updated():
    """ 
    Initializes each worker process. This function loads the models 
    only once when the process starts.
    """

    # from .loadmodelrpf_retrained import model as m1, scalerparam as s1, trunks as t1, minmaxReverse as mm1re, minmaxscale as mm1
    # from .loadmodelcore_retrained import model as m2, scalerparam as s2, trunks as t2
    # from .loadmodelpin import best_model as pin_m
    
    global worker_data
    
    import gc
    if 'worker_data' in globals():
        del worker_data
    gc.collect()
    worker_data = {}
    current_file_path = Path(__file__).resolve()
    source_path = str(current_file_path.parent)
    if source_path not in sys.path:
        sys.path.insert(0, source_path)
    import loadmodelrpf_retrained
    import loadmodelcore_retrained
    import loadmodelpin
    importlib.reload(loadmodelrpf_retrained)
    importlib.reload(loadmodelcore_retrained)
    importlib.reload(loadmodelpin)

    worker_data['model'] = loadmodelrpf_retrained.model
    worker_data['scalerparam'] = loadmodelrpf_retrained.scalerparam
    worker_data['trunks'] = loadmodelrpf_retrained.trunks
    worker_data['minmaxscale']=loadmodelrpf_retrained.minmaxscale
    worker_data['minmaxReverse']=loadmodelrpf_retrained.minmaxReverse
    worker_data['modelcore'] = loadmodelcore_retrained.model
    worker_data['scaler_core'] = loadmodelcore_retrained.scalerparam
    worker_data['trunk_core'] = loadmodelcore_retrained.trunks
    worker_data['pinmodel'] = loadmodelpin.best_model
    
    try:
      cpu_id = os.sched_getcpu()
    except AttributeError:
       cpu_id = "N/A"
    pid = os.getpid()
    hostname = socket.gethostname()
    slurm_task = os.environ.get("SLURM_PROCID", "N/A")
    slurm_localid = os.environ.get("SLURM_LOCALID", "N/A")
    
    print(
        f"[Node: {hostname}] "
        f"[PID: {pid}] "
        f"[CPU: {cpu_id}] "
        f"[SLURM_PROCID: {slurm_task}] "
        f"[SLURM_LOCALID: {slurm_localid}] "
        f"successfully loaded all 3 models."
    )

def getXSlist(xsdict, fatype, xstype, axialnode):
    return [xsdict[fatype][key][xstype][axialnode] for key in xsdict[fatype]]

def depletion_power(model, modelinput, scaler,idx1, idx2, total_height, faaxial):
    '''
    Predict the new powdistribution based on new burnup distribution for 193 core
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
    idx = np.r_[8, 17, 25, 26, 35, 43:45, 52:54, 60:63, 67:tmp_pred.shape[1]]

    # tmp_pred[abs(tmp_pred)<1e-4] = 0
    tmp_pred = normalize(tmp_pred,idx2, idx1,total_height, faaxial)
    # y_predr = minmaxReverse(tmp_pred.reshape(tmp_pred.shape[0],newpow.shape[-1]),scaler['omax'],scaler['omin'])
    y_predr = tmp_pred.reshape(newpow.shape)
    return y_predr

def depletion_power_157(model, modelinput, scaler,idx1, idx2, total_height, faaxial):
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
    tmp_pred = normalize(tmp_pred,idx2, idx1,total_height, faaxial)
    # y_predr = minmaxReverse(tmp_pred.reshape(tmp_pred.shape[0],newpow.shape[-1]),scaler['omax'],scaler['omin'])
    y_predr = tmp_pred.reshape(newpow.shape)
    return y_predr

def depletion_boroncycle(model, modelinput, scaler, minmaxReverse):
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

    ## test 
    idx = next((i for i in range(len(X)-1) if abs(X[i+1][0] - X[i][0]) < 10 and X[i][0]<500) , -1)
    idx=idx-4 ## take 3 step back 


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

def getdata_fully_vectorized(bumap, xsdict, coremap, fabulist, scalerparam, minmaxscale):
    """
    Fully vectorized version that pre-computes data structures to minimize
    work inside the main loop.
    Includes an optional validation step against np.interp.
    Optimized version of the interpolation process
    bumap for interpolation
    xsdict for data of xs 
    scalerparam for scaling XS max, min value 
    coremap for LPs that maps FA/reflector/void
    fabulist for interpolation as the XS data is the function of these burnup value
    minmaxscale is the function for maxminscaling (+buffer values)
    '''
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

def compute_tempnorm(output,idx2, idx1,total_height, faaxial):
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
    count_per_sample = 193
    # Final normalization
    tempnorm = (tolRFP / count_per_sample).reshape(-1, 1)  # shape (N, 1)
    return tempnorm

def normalize(output,idx2, idx1,total_height, faaxial):
    tempnorm = compute_tempnorm(output,idx2, idx1,total_height, faaxial)
    outnew = [output[i]*1.0/tempnorm[i] for i in range(len(output))]
    outnew = np.array(outnew).reshape(output.shape)
    return outnew


# --- MAIN EXECUTION FUNCTION ---

def get_result(LPs, corebulist, xsdict, fabulist, faaxial, total_height, idx11, idx22) :
    """
    Function executed by worker processes.
    """
    # Access the process-specific models
    global worker_data
    model = worker_data['model']
    modelcore = worker_data['modelcore']
    pinmodel = worker_data['pinmodel']
    scalerparam = worker_data['scalerparam']
    trunks = worker_data['trunks']
    scaler_core = worker_data['scaler_core']
    trunk_core = worker_data['trunk_core']
    minmaxReverse = worker_data['minmaxReverse']
    minmaxscale = worker_data['minmaxscale']
    
    pinpower =[]
    test = " ".join(LPs).strip().split()
    ## initializing 
    initbumap = np.zeros((1,81,16))
    inputdata = getdata_fully_vectorized(initbumap, xsdict, test, fabulist, scalerparam, minmaxscale)
    pow_0 = depletion_power(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22, total_height, faaxial)
    pinpowerbu = getfqFd_pinrecontruct(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    pinpower.append(pinpowerbu)
    bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, tuple(inputdata+[trunk_core]), scaler_core,minmaxReverse)
    # bustep+=1
    # Initialize
    power_history = [pow_0]
    boron_his = [bor_pred0]
    cycle_his = [cyc_pred0]
    for bu in corebulist:
            converged = False
            for iteration in range(100):
                if iteration == 0:
                    delE = bu * pow_0
                    bumap = initbumap + delE.reshape((1, 81, 16))
                else:
                    delE = bu * (0.5 * pow_0 + 0.5 * pow_prev)
                    bumap = initbumap + delE.reshape((1, 81, 16))
                input_data = getdata_fully_vectorized(bumap, xsdict, test, fabulist, scalerparam, minmaxscale)
                pow_new = depletion_power(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22, total_height,faaxial)
                # Check convergence
                if iteration > 0:
                    max_delta = np.max(np.abs(pow_new - pow_prev))
                    if max_delta < 1e-4:
                        pinpowerbu= getfqFd_pinrecontruct(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        pinpower.append(pinpowerbu)
                        bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, 
                                                                    tuple(input_data+[trunk_core]), 
                                                                    scaler_core, minmaxReverse)
                        # print(f"BU={bu}, iter={iteration} → Converged (tol={max_delta:.2e})")
                        # Update for next burnup
                        initbumap = bumap
                        power_history.append(pow_new)
                        boron_his.append(bor_pred0)
                        cycle_his.append(cyc_pred0)
                        pow_0 = pow_new
                        converged = True
                        # bustep+=1
                        break
        
                pow_prev = pow_new
        
            if not converged:
                print(f"BU={bu} → ❌ Unconverged after 100 iterations")
                raise ValueError("Solution did not converge.")
                stop
    ## get the Fdmax and Fqmax here 
    pinpower = np.array(pinpower).reshape(-1,153,153,1)
    pinpower_reconstruct = pinmodel.predict(pinpower, verbose=0) # verbose=0 suppresses Keras logs for speed
    # 6. Apply Mask and Calculate Statistics
    # Create mask from the first node (assuming geometry is axially constant)
    # Using np.not_equal is slightly faster/cleaner than != 0
    keep_mask = np.not_equal(pinpower[0], 0) 
    
    # Broadcast multiply: (16,153,153,1) * (153,153,1)
    pinpower_reconstruct *= keep_mask
    # Calculate Max Fq
    Fq_all = np.max(pinpower_reconstruct)
    
    # Calculate Max Fdel (Fdelta-H equivalent?)
    # Pre-calculate axial shape factor
    ax_factor = (np.array(faaxial) / total_height).reshape(1, 16, 1, 1)
    
    # Multiply and find max
    Fd_all = np.max(np.sum(pinpower_reconstruct.reshape(34,16,153,153) * ax_factor,axis = 1))
    cycle_length = interpolatecycle(10,boron_his,cycle_his)
    if cycle_length>1000:
       print('boron list',boron_his)
       print('cycle list',cycle_his)
       print(LPs)
       

    return Fd_all, Fq_all, np.max(boron_his), cycle_length

def get_result_157(LPs, corebulist, xsdict, fabulist, faaxial, total_height, idx11, idx22) :
    """
    Function executed by worker processes for core 157 FAs.
    """
    # Access the process-specific models
    global worker_data
    model = worker_data['model']
    modelcore = worker_data['modelcore']
    pinmodel = worker_data['pinmodel']
    scalerparam = worker_data['scalerparam']
    trunks = worker_data['trunks']
    scaler_core = worker_data['scaler_core']
    trunk_core = worker_data['trunk_core']
    minmaxReverse = worker_data['minmaxReverse']
    minmaxscale = worker_data['minmaxscale']
    
    pinpower =[]
    test = " ".join(LPs).strip().split()
    ## initializing 
    initbumap = np.zeros((1,81,16))
    inputdata = getdata_fully_vectorized(initbumap, xsdict, test, fabulist, scalerparam, minmaxscale)
    pow_0 = depletion_power_157(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22, total_height, faaxial)
    pinpowerbu = getfqFd_pinrecontruct(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    pinpower.append(pinpowerbu)
    bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, tuple(inputdata+[trunk_core]), scaler_core,minmaxReverse)
    # bustep+=1
    # Initialize
    power_history = [pow_0]
    boron_his = [bor_pred0]
    cycle_his = [cyc_pred0]
    for bu in corebulist:
            converged = False
            for iteration in range(100):
                if iteration == 0:
                    delE = bu * pow_0
                    bumap = initbumap + delE.reshape((1, 81, 16))
                else:
                    delE = bu * (0.5 * pow_0 + 0.5 * pow_prev)
                    bumap = initbumap + delE.reshape((1, 81, 16))
                input_data = getdata_fully_vectorized(bumap, xsdict, test, fabulist, scalerparam, minmaxscale)
                pow_new = depletion_power_157(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22, total_height,faaxial)
                # Check convergence
                if iteration > 0:
                    max_delta = np.max(np.abs(pow_new - pow_prev))
                    if max_delta < 1e-4:
                        pinpowerbu= getfqFd_pinrecontruct(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        pinpower.append(pinpowerbu)
                        bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, 
                                                                    tuple(input_data+[trunk_core]), 
                                                                    scaler_core, minmaxReverse)
                        # print(f"BU={bu}, iter={iteration} → Converged (tol={max_delta:.2e})")
                        # Update for next burnup
                        initbumap = bumap
                        power_history.append(pow_new)
                        boron_his.append(bor_pred0)
                        cycle_his.append(cyc_pred0)
                        pow_0 = pow_new
                        converged = True
                        # bustep+=1
                        break
        
                pow_prev = pow_new
        
            if not converged:
                print(f"BU={bu} → ❌ Unconverged after 100 iterations")
                raise ValueError("Solution did not converge.")
                stop
    ## get the Fdmax and Fqmax here 
    pinpower = np.array(pinpower).reshape(-1,153,153,1)
    pinpower_reconstruct = pinmodel.predict(pinpower, verbose=0) # verbose=0 suppresses Keras logs for speed
    # 6. Apply Mask and Calculate Statistics
    # Create mask from the first node (assuming geometry is axially constant)
    # Using np.not_equal is slightly faster/cleaner than != 0
    keep_mask = np.not_equal(pinpower[0], 0) 
    
    # Broadcast multiply: (16,153,153,1) * (153,153,1)
    pinpower_reconstruct *= keep_mask
    # Calculate Max Fq
    Fq_all = np.max(pinpower_reconstruct)
    
    # Calculate Max Fdel (Fdelta-H equivalent?)
    # Pre-calculate axial shape factor
    ax_factor = (np.array(faaxial) / total_height).reshape(1, 16, 1, 1)
    
    # Multiply and find max
    Fd_all = np.max(np.sum(pinpower_reconstruct.reshape(34,16,153,153) * ax_factor,axis = 1))
    cycle_length = interpolatecycle(10,boron_his,cycle_his)
    if cycle_length>1000:
       print('boron list',boron_his)
       print('cycle list',cycle_his)
       print(LPs)
       

    return Fd_all, Fq_all, np.max(boron_his), cycle_length


# def retrained_models_1(datafile):
#     from .loadmodel_Nov25 import model as m1, scalerparam as s1, trunks as t1, minmaxReverse as mm1re, minmaxscale as mm1
#     return 

# def retrained_models_2():
#     return 

# def retrained_models_3():
#     return 


def get_result_serial(LPs, corebulist, xsdict, fabulist, faaxial, total_height, idx11, idx22) :
    """
    Function executed by worker processes.
    """
    
    model = m1
    modelcore = m2
    pinmodel = pin_m
    scalerparam = s1
    trunks = t1
    scaler_core = s2
    trunk_core = t2
    minmaxReverse = mm1re
    minmaxscale = mm1
    
    pinpower =[]
    test = " ".join(LPs).strip().split()
    ## initializing 
    initbumap = np.zeros((1,81,16))
    inputdata = getdata_fully_vectorized(initbumap, xsdict, test, fabulist, scalerparam, minmaxscale)
    pow_0 = depletion_power(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22, total_height, faaxial)
    pinpowerbu = getfqFd_pinrecontruct(initbumap, pow_0, test, fabulist, xsdict, faaxial, total_height) # first step
    pinpower.append(pinpowerbu)
    bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, tuple(inputdata+[trunk_core]), scaler_core,minmaxReverse)
    # bustep+=1
    # Initialize
    power_history = [pow_0]
    boron_his = [bor_pred0]
    cycle_his = [cyc_pred0]
    for bu in corebulist:
            converged = False
            for iteration in range(100):
                if iteration == 0:
                    delE = bu * pow_0
                    bumap = initbumap + delE.reshape((1, 81, 16))
                else:
                    delE = bu * (0.5 * pow_0 + 0.5 * pow_prev)
                    bumap = initbumap + delE.reshape((1, 81, 16))
                input_data = getdata_fully_vectorized(bumap, xsdict, test, fabulist, scalerparam, minmaxscale)
                pow_new = depletion_power(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22, total_height,faaxial)
                # Check convergence
                if iteration > 0:
                    max_delta = np.max(np.abs(pow_new - pow_prev))
                    if max_delta < 1e-4:
                        pinpowerbu= getfqFd_pinrecontruct(bumap, pow_new, test, fabulist, xsdict, faaxial, total_height)
                        pinpower.append(pinpowerbu)
                        bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, 
                                                                    tuple(input_data+[trunk_core]), 
                                                                    scaler_core, minmaxReverse)
                        # print(f"BU={bu}, iter={iteration} → Converged (tol={max_delta:.2e})")
                        # Update for next burnup
                        initbumap = bumap
                        power_history.append(pow_new)
                        boron_his.append(bor_pred0)
                        cycle_his.append(cyc_pred0)
                        pow_0 = pow_new
                        converged = True
                        # bustep+=1
                        break
        
                pow_prev = pow_new
        
            if not converged:
                print(f"BU={bu} → ❌ Unconverged after 100 iterations")
                raise ValueError("Solution did not converge.")
    ## get the Fdmax and Fqmax here 
    pinpower = np.array(pinpower).reshape(-1,153,153,1)
    pinpower_reconstruct = pinmodel.predict(pinpower, verbose=0) # verbose=0 suppresses Keras logs for speed
    # 6. Apply Mask and Calculate Statistics
    # Create mask from the first node (assuming geometry is axially constant)
    # Using np.not_equal is slightly faster/cleaner than != 0
    keep_mask = np.not_equal(pinpower[0], 0) 
    
    # Broadcast multiply: (16,153,153,1) * (153,153,1)
    pinpower_reconstruct *= keep_mask
    # Calculate Max Fq
    Fq_all = np.max(pinpower_reconstruct)
    
    # Calculate Max Fdel (Fdelta-H equivalent?)
    # Pre-calculate axial shape factor
    ax_factor = (np.array(faaxial) / total_height).reshape(1, 16, 1, 1)
    
    # Multiply and find max
    Fd_all = np.max(np.sum(pinpower_reconstruct.reshape(34,16,153,153) * ax_factor,axis = 1))
    cycle_length = interpolatecycle(10,boron_his,cycle_his)[0]

    return Fd_all, Fq_all, np.max(boron_his), cycle_length