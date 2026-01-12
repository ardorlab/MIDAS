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


### import surrogate model here 
import itertools
import pickle 
import copy
import time as TT
import matplotlib.pyplot as plt
from surmodel.loadmodel_Nov25 import minmaxReverse, minmaxscale
from surmodel.loadmodel_Nov25 import model
from surmodel.loadmodel_Nov25 import scalerparam
from surmodel.loadmodel_Nov25 import trunks
from surmodel.loadmodelcore import trunks as trunk_core
from surmodel.loadmodelcore import scalerparam as scaler_core
from surmodel.loadmodelcore import model as modelcore
from surmodel.loadmodelpin import best_model as pinmodel
import tensorflow as tf


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

## require function for depletion module 
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
    tmp_pred[:, 68:73, :] = 0
    tmp_pred[:, 73:, :] = 0
    # tmp_pred[abs(tmp_pred)<1e-4] = 0
    tmp_pred = normalize(tmp_pred,idx2, idx1,total_height, faaxial)
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
    ## only take 2 last element for cycle length calculation 
    a = (Y[-2]-Y[-1])/(X[-2]-X[-1])
    b = Y[-2] - a*X[-2]
    y = a*x+b 
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

## core map

xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/axialFAtypeXS.pkl','rb'))
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
idx11 = np.array([0])
idx22 = np.array([1, 2, 3, 4, 5, 6, 7, 9, 18, 27, 36, 40, 54, 63])

listFA = [461,462,501,502,526,566,586]
corebulist = [0.1, 0.4, 0.5,1.0,1.0]
for i in range(28):
    corebulist.append(1*1.0)
def depltion_cal():
    return 

def comparewithPARCS():

    return

def retrain():

    return



## Functions ##
def evaluate_old(solution, input):
    """
    Interface used to run PARCSv343 calculations.
    
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
    
        filename = solution.name + '.inp'
        # create input file based on if an input template is given
        if not input.input_template['apply']: 
            without_template(solution, input, cwd, filename)
        elif input.input_template['apply']:
            with_template(solution, input, cwd, filename)

    ## Run PARCS INPUT DECK #!TODO: separate the input writing and execution into two different functions that are called in sequence.
        parcscmd = __parcs343exe__
        print(solution.parameters)
        stop
        try:
            output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
        ## Get Results
            if 'Finished' in str(output): #job completed
                logger.debug(f"Job {solution.name} completed successfully in PARCSv343.")
                solution.parameters = get_results(solution.parameters, solution.name)
                print(solution.parameters)
                stop
            
            else: #job failed
                if input.calculation_type in ['eq_cycle']:
                    try:
                        solution.parameters = eq_cycle_convergence(input, solution, filename, parcscmd, input.code_walltime) #iteratively try to find an intial guess that will converge
                    except Exception as e:
                        logger.error(f"Job {solution.name} has failed to converge with the following exception: {e}")
                        solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
                        
                else: #standard execution pathway
                    logger.warning(f"Job {solution.name} has failed!")
                    solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        
        except subprocess.TimeoutExpired: #job timed out
            os.system('rm -f {}.parcs_pin*'.format(solution.name))
            logger.error(f"Job {solution.name} has timed out!")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        except subprocess.CalledProcessError as e: #PARCS returned an abort signal
            logger.error(f"Job {solution.name} has failed with the following exception: {e}")
            solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        
        logger.debug(f"Returning to original working directory: {cwd}")
        os.chdir(cwd)
        gc.collect()
        print(solution.__dir__())
        stop

        return solution

def evaluate(solution, input):
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
    
        filename = solution.name + '.inp'
        # create input file based on if an input template is given
        if not input.input_template['apply']: 
            LPs = without_template(solution, input, cwd, filename)
        elif input.input_template['apply']:
            with_template(solution, input, cwd, filename)


        ## create input then take the LPs here 

    ## Run PARCS INPUT DECK #!TODO: separate the input writing and execution into two different functions that are called in sequence.
        parcscmd = __parcs343exe__
        ### surrogate depletion here 
        # print(LPs)
        LPs = "".join(LPs).strip().split()
        ## initializing 
        initbumap = np.zeros((1,81,16))
        Fq_all = 0
        Fd_all = 0
        pinpower =[]
        inputdata = getdata_fully_vectorized(initbumap, xsdict, LPs, fabulist, scalerparam)
        pow_0 = depletion_power(model, tuple(inputdata+[trunks]), scalerparam,idx11, idx22)
        pinpowerbu = getfqFd_pinrecontruct(initbumap, pow_0, LPs, fabulist, xsdict, faaxial, total_height) # first step
        pinpower.append(pinpowerbu)
        bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, tuple(inputdata+[trunk_core]), scaler_core)
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
                    input_data = getdata_fully_vectorized(bumap, xsdict, LPs, fabulist, scalerparam)
                    pow_new = depletion_power(model, tuple(input_data+[trunks]), scalerparam,idx11, idx22)
                    # Check convergence
                    if iteration > 0:
                        max_delta = np.max(np.abs(pow_new - pow_prev))
                        if max_delta < 1e-4:
                            pinpowerbu= getfqFd_pinrecontruct(bumap, pow_new, LPs, fabulist, xsdict, faaxial, total_height)
                            pinpower.append(pinpowerbu)
                            bor_pred0, cyc_pred0 = depletion_boroncycle(modelcore, 
                                                                        tuple(input_data+[trunk_core]), 
                                                                        scaler_core)
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

        # print( Fd_all, Fq_all, np.max(boron_his), interpolatecycle(10,boron_his,cycle_his))
        # param = ["cycle_length", "pinpowerpeaking", "fdeltah", "max_boron"]
        solution.parameters["cycle_length"]['value'] = interpolatecycle(10,boron_his,cycle_his)[0]
        solution.parameters["pinpowerpeaking"]['value'] = Fq_all
        solution.parameters["fdeltah"]['value'] = Fd_all
        solution.parameters["max_boron"]['value'] = np.max(boron_his)

        # try:
        #     output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=input.code_walltime) #wait until calculation finishes
        # ## Get Results
        #     if 'Finished' in str(output): #job completed
        #         logger.debug(f"Job {solution.name} completed successfully in PARCSv343.")
        #         solution.parameters = get_results(solution.parameters, solution.name)
            
        #     else: #job failed
        #         if input.calculation_type in ['eq_cycle']:
        #             try:
        #                 solution.parameters = eq_cycle_convergence(input, solution, filename, parcscmd, input.code_walltime) #iteratively try to find an intial guess that will converge
        #             except Exception as e:
        #                 logger.error(f"Job {solution.name} has failed to converge with the following exception: {e}")
        #                 solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
                        
        #         else: #standard execution pathway
        #             logger.warning(f"Job {solution.name} has failed!")
        #             solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        
        # except subprocess.TimeoutExpired: #job timed out
        #     os.system('rm -f {}.parcs_pin*'.format(solution.name))
        #     logger.error(f"Job {solution.name} has timed out!")
        #     solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        # except subprocess.CalledProcessError as e: #PARCS returned an abort signal
        #     logger.error(f"Job {solution.name} has failed with the following exception: {e}")
        #     solution.parameters = get_results(solution.parameters, solution.name, job_failed=True)
        
        logger.debug(f"Returning to original working directory: {cwd}")
        os.chdir(cwd)
        gc.collect()

        return solution

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
        #!ofile.write("      PLOT_OPTS 0 0 0 0 0 2\n")
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