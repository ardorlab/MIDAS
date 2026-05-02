import os
import random
import itertools
import numpy as np
import pickle
import h5py  # <-- Added for saving data iteratively

# --- Helper Functions (Cleaned Up) ---

def extractpin(lines):
    """Extracts PARCS pin data for one FA across all BU/axial nodes."""
    # This array holds data for ONE FA: [34 BU steps, 17x, 17y, 16 axial]
    pin = np.zeros((34, 17, 17, 16))
    xloc, yloc = None, None
    
    Nline = len(lines)
    for idx in range(Nline):
        if lines[idx].find('Case') >= 0:
            temp = lines[idx].split()
            buidx_str = temp[1]
            axialidx_str = temp[10]

            if xloc is None:  # Get location from the first 'Case' line
                xloc = int(temp[8])
                yloc = int(temp[9])
            
            # get the block
            block_lines = lines[idx + 2 : idx + 2 + 17]
            block = []
            for line in block_lines:
                # Split line into floats, ignore the first column (row number)
                nums = [float(x) for x in line.split()[1:]]
                block.append(nums)
            
            arr = np.array(block)
            
            # Store in 0-based indices
            # File axial '2'-'17' -> array index 0-15
            # File BU '1'-'34' -> array index 0-33
            if axialidx_str not in ['18', '1', '0']:
                buidx_0 = int(buidx_str) - 1
                axialidx_0 = int(axialidx_str) - 2
                pin[buidx_0, :, :, axialidx_0] = arr
            else:
                pass
    return xloc, yloc, pin

def parse_dep_file(ofile, nfa, z_id, refl_id):
    """
    Parses a PARCS .dep file ONCE to get both RPF and Burnup 3D maps.
    """
    with open(ofile, 'r') as f:
        txt = f.readlines()

    nz = len(z_id)
    map_keys = {
        'rpf': (' RPF 3D MAP', ' EXP 2D MAP'),
        'bu': ('EXP 3D MAP', ' I_D 2D MAP'),
    }
    
    results = {}

    for map_type, (k_sta, k_end) in map_keys.items():
        alldeplines = []
        i_sta = -1
        for i, line in enumerate(txt):
            if line.find(k_sta) >= 0:
                i_sta = i + 1
            if line.find(k_end) >= 0 and i_sta != -1:
                i_end = i - 1
                depline = "".join(txt[i_sta : i_end + 1])
                alldeplines.append(depline)
                i_sta = -1 # Reset for next block

        nbu = len(alldeplines)
        if nbu == 0:
            # Return empty arrays if no data was found
            results[map_type] = np.array([])
            continue

        data_3d = np.zeros((nbu, nfa, nz))
        
        for bu, step in enumerate(alldeplines):
            txt_dep = step.split('\n')
            txt_dep = [line for line in txt_dep if line.strip()] # Cleaner filter
            
            asb_counter = 0
            fasb_counter = 0
            ifass = 0
            iass = 0
            
            for i in range(1, len(txt_dep)):
                line_dep = txt_dep[i].split()
                if line_dep[0] == 'k':
                    asb_counter = iass  # No copy.deepcopy needed
                    fasb_counter = ifass
                elif int(line_dep[0]) not in z_id:
                    pass
                else:
                    iz = z_id[-1] - int(line_dep[0])
                    ifass = fasb_counter
                    iass = asb_counter
                    for j in range(1, len(line_dep)):
                        iass += 1
                        if iass in refl_id:
                            pass
                        else:
                            ifass += 1
                            val = float(line_dep[j])
                            data_3d[bu, ifass - 1, iz] = val
        results[map_type] = data_3d

    return results.get('rpf', np.array([])), results.get('bu', np.array([]))


def getXSlist(xsdict, fatype, xstype, axialnode):
    bukyes = xsdict[fatype].keys()
    xsval = [xsdict[fatype][key][xstype][axialnode] for key in bukyes]
    return xsval

def interp_nd(fabulist, pp, buval):
    """
    Interpolates pp values at query points buval based on fabulist.
    (Original function was correct)
    """
    fabulist = np.asarray(fabulist)
    pp = np.asarray(pp)
    buval = np.atleast_1d(buval)  # ensure array

    # Find interval indices
    idx = np.searchsorted(fabulist, buval) - 1
    idx = np.clip(idx, 0, len(fabulist) - 2)

    x0, x1 = fabulist[idx], fabulist[idx + 1]  # (M,)
    y0, y1 = pp[idx, :], pp[idx + 1, :]  # (M, F)

    # Handle divide-by-zero if x1 == x0
    dx = x1 - x0
    t = np.zeros_like(buval, dtype=float)
    np.divide(buval - x0, dx, out=t, where=dx != 0) # (M,)

    pp_interp = y0 + (y1 - y0) * t[:, None]  # (M, F)

    # If input was scalar → return 1D (F,)
    if np.isscalar(buval) or np.ndim(buval) == 0:
        return pp_interp[0]
    return pp_interp

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

# --- Main Script ---

pathtoextract = "/home/khnguy22/Deeponet-midas/MIDAS/corereload197-surrogate-nofb_new/output_files"
output_folder = pathtoextract

listFA = [250, 280, 461,462,501,502,526,566,586]
fabulist = [
    0, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.5, 15.0, 17.5, 20.0,
    25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0
]
fabulist = np.array(fabulist)
xsdict = pickle.load(open('/home/khnguy22/Deeponet-midas/testparcs/updateXS.pkl', 'rb'))

# --- Geometry parameters for parser ---
nfa_parcs = 56 # Number of non-reflector FAs
# nfa_parcs = 47 # Number of non-reflector FAs for 157core 
z_id_parcs = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
nz_parcs = len(z_id_parcs)
refl_id_parcs = [9, 18, 27, 36, 44, 45, 53, 60, 61, 66, 67, 68, 69, 70, 71, 72, 73]
# refl_id_parcs = [9, 18, 26, 27, 35, 42, 43, 49, 50, 55, 56, 59, 60, 61, 62, 63, 64] # for 157 core
# (Assuming FULL_CORE=False based on original code)

# --- REFACTOR: Initialize HDF5 file ---
# This file will store the output datasets on disk
output_h5_file = 'pin_data_193r_v51.h5'
hf = h5py.File(output_h5_file, 'w')

# Create resizable datasets. We will append data in chunks.
# (10, 153, 153) is a good chunk size for performance.
dset_parcs = hf.create_dataset(
    'parcpin', 
    (0, 153, 153), 
    maxshape=(None, 153, 153), 
    chunks=(10, 153, 153), 
    dtype='f'
)
dset_calcs = hf.create_dataset(
    'calpin', 
    (0, 153, 153), 
    maxshape=(None, 153, 153), 
    chunks=(10, 153, 153), 
    dtype='f'
)
# ---

# N_SAMPLES = 544 # Number of random samples to take per file

for root, dirs, files in os.walk(pathtoextract):
    # print(dirs, root)
    for file in files:
        if file.endswith('.inp'):
            inputfile = os.path.join(root, file)
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
            with open(outputdpl_path) as f:
                lines = f.readlines()
            rpf3d, bu3d = parse_dep_file(
                outputdep_path, nfa_parcs, z_id_parcs, refl_id_parcs
            )
            nbu = len(bu3d)
            # test = getLP(output_path)  
            test = getLPs_cyc(outputdep_path)  # for multi_cycle
            nofluxpin = np.zeros((nfa_parcs, nbu, 17, 17, nz_parcs))
            for bu_step in range(nbu):
                idx = 0  # Fuel assembly index (0 to 55)
                for loc, fa in enumerate(test):
                    if fa not in ['10', '00']:
                        for node in range(nz_parcs):
                            buval = bu3d[bu_step, idx, node]
                            nodepower = rpf3d[bu_step, idx, node]
                            
                            pp = getXSlist(xsdict, fa, 'ppowers', str(node))
                            pindata = interp_nd(fabulist, pp, buval).flatten() * nodepower
                            
                            nofluxpin[idx, bu_step, :, :, node] = pindata.reshape(17, 17)
                        idx += 1
            
            # 3. Get random samples
            # extractidx_flat = random.sample(list(range(nbu * nz_parcs)), N_SAMPLES)
            extractidx_flat = list(range(nbu * nz_parcs))
            N_SAMPLES = len(extractidx_flat) # take all 
            # These are the 4 maps (153x153) we will build and save
            parcs_maps_to_save = [np.zeros((153, 153)) for _ in range(N_SAMPLES)]
            calcs_maps_to_save = [np.zeros((153, 153)) for _ in range(N_SAMPLES)]
            append_flag = True

            # 4. Process one FA pin file at a time (low memory)
            for fa_file_idx in range(nfa_parcs): # 0 to 55
                pin_file_num = fa_file_idx + 1 # 1 to 56
                pin_file_name = f"{pinfilename}{pin_file_num:03}"
                pin_file_path = os.path.join(folder_path, pin_file_name)
                if not os.path.exists(pin_file_path):
                    print(f"{pinfilename} not existed?")
                    append_flag = False
                    break
                with open(pin_file_path) as f:
                    pinreadlines = f.readlines()
                
                # Extract all data for this single FA
                xloc, yloc, pin_data_all_fa = extractpin(pinreadlines)
                # pin_data_all_fa shape is (34, 17, 17, 16)
                
                if xloc is None:
                    print(f"Warning: Could not parse xloc/yloc from {pin_file_name}")
                    continue

                x_start, y_start = (xloc - 1) * 17, (yloc - 1) * 17
                x_end, y_end = x_start + 17, y_start + 17

                for i in range(N_SAMPLES):
                    idx_flat = extractidx_flat[i]
                    bu_idx = idx_flat // nz_parcs # 0-based BU index
                    ax_idx = idx_flat % nz_parcs  # 0-based axial index

                    # Get data from PARCS file
                    parcs_pin_slice = pin_data_all_fa[bu_idx, :, :, ax_idx]
                    parcs_maps_to_save[i][y_start:y_end, x_start:x_end] = parcs_pin_slice
                    
                    # Get data from our calculated 'nofluxpin'
                    calcs_pin_slice = nofluxpin[fa_file_idx, bu_idx, :, :, ax_idx]
                    calcs_maps_to_save[i][y_start:y_end, x_start:x_end] = calcs_pin_slice
            
            if append_flag: ## only append if non zero
                parcs_data_batch = np.array(parcs_maps_to_save) 
                calcs_data_batch = np.array(calcs_maps_to_save) 
                ## revert axial index 
                N = calcs_data_batch.shape[0] 
                calcs_data_batch = calcs_data_batch.reshape(1,nbu,nz_parcs,153,153)
                print(calcs_data_batch.shape, nbu,nz_parcs)
                calcs_data_batch = calcs_data_batch[:,:,::-1,:,:]
                calcs_data_batch = calcs_data_batch.reshape(N,153,153)
                
                
                dset_parcs.resize(dset_parcs.shape[0] + N_SAMPLES, axis=0)
                dset_parcs[-N_SAMPLES:] = parcs_data_batch
                
                dset_calcs.resize(dset_calcs.shape[0] + N_SAMPLES, axis=0)
                dset_calcs[-N_SAMPLES:] = calcs_data_batch
    
                print(f"{inputfile} ... Done.")
            
# --- Finalize ---
hf.close()  # Close the HDF5 file
print("✅ All perturbations processed and saved to pin_data.h5.")
