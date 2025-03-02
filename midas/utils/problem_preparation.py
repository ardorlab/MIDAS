## Import Block ##
import shutil
import os
import numpy as np


## Classes ##
class Problem_Preparation_Tools():
    """
    Class of functions necessary for formatting the geometry and labels from the input data.
    
    Written by Nicholas Rollins. 10/04/2024
    """
    def parse_FAlabels(input_obj): #!TODO: where should this be called from?
        """
        
        Written by Nicholas Rollins. 10/04/2024
        """
        xs_list = {'fuel':[], 'reflectors':{'top':[], 'radial':[], 'bot':[]}, 'blankets':[]}
        tag_list = {'fuel':[], 'reflectors':[]}
        
        for key, param in input_obj.fa_options['fuel'].items():
            xs_list['fuel'].append(param['serial'])
            tag_list['fuel'].append((str(param['type']),key))
        for key, param in input_obj.fa_options['reflectors'].items():
            if param['refl_type'] == 'all':
                xs_list['reflectors']['top'].append(param['serial'])
                xs_list['reflectors']['radial'].append(param['serial'])
                xs_list['reflectors']['bot'].append(param['serial'])
                tag_list['reflectors'].append((str(10+len(tag_list['reflectors'])),key))
            elif param['refl_type'] == 'top':
                xs_list['reflectors']['top'].append(param['serial'])
            elif param['refl_type'] == 'radial':
                xs_list['reflectors']['radial'].append(param['serial'])
                tag_list['reflectors']['radial'].append((str(10+len(tag_list['reflectors'])),key))
            elif param['refl_type'] == 'bot':
                xs_list['reflectors']['bot'].append(param['serial'])
        if 'blankets' in input_obj.fa_options:
            for key, param in input_obj.fa_options['blankets'].items():
                xs_list['blankets'].append(param['serial'])
        
        return xs_list, tag_list
    
    def generate_core(input_obj, core_map):
        """
        Generates the e.g. 17x17 core map with consistent identifiers and treatment of symmetry.

        Parameters: None
        Additional comments:
          - The core_map is manually defined as a list line by line starting from the top of the core.
          It is advised to use the following naming conention for reflector assemblies 'ABCDE' 
          with A being R, BC indicating the row number (00-17) and DE the column number (00-17). For 
          the fuel assemblies it is advised to use the following naming convention 'ABC' with A being the 
          row letter indetifier (A-O) and BC being the column number (00-17). 

        Written by Gregory Delipei. 7/12/2022
        Updated by Nicholas Rollins. 10/04/2024
        """
        nrows = input_obj.nrow
        ncols = input_obj.ncol
        
        core_id = []
        for i in range(nrows-1,-1,-1): #iterate over row numbers in reverse
            for j in range(ncols):
                core_id.append((i-np.floor(nrows/2),j-np.floor(ncols/2)))

        core_id = np.array(core_id).reshape((nrows,ncols,2))
        
        core_sym_map = Problem_Preparation_Tools.symmetric_core(input_obj.symmetry,nrows,ncols,core_map,core_id)
        
        return core_sym_map

    def symmetric_core(symmetry,nrows,ncols,core_map,core_id):
        """
        Generates the symmetric core map.
        Handles even cores (BWR) and odd cores (PWR) separately

        Parameters: 
           - core_map: an e.g. 17x17 numpy array with the fuel assembly location names.
           - core_id: an e.g. 17x17x2 numpy array with coordinate indices for each fuel assembly location
           ranging from e.g. -8 to +8.

        Written by Gregory Delipei 7/12/2022
        Updated by Nicholas Rollins. 10/04/2024
        Updated by Jake Mikouchi 3/2/25
        """
        #!TODO: I think this assumes only one layer of radial reflectors.
        sym_center   = (np.floor(nrows/2),np.floor(ncols/2)) #x,y-coordinates for centerpoint of symmetry
        sym_vertical = (nrows-1,np.floor(ncols/2))
        if symmetry == 'quarter':
            sym_corner   = (np.floor(nrows/2),ncols-1)
        elif symmetry == 'octant':
            sym_corner   = (nrows-1,ncols-1)
        
        row_iter = np.arange(sym_center[0],sym_vertical[0]+1,1)
        col_iter = np.arange(sym_center[1],sym_corner[1]+1,1)
        
        if nrows % 2 == 0: # find symetric locations for cores with even number of rows and cols (BWRS)
            if symmetry == 'quarter':
                core_dict={}
                for irow in row_iter:
                    for icol in col_iter:
                        dict_value={'Symmetric_Assemblies':[], 'Value': None}

                        # finds assembly location relative to center assembly
                        irow = int(irow); icol = int(icol)
                        true_row_center = nrows / 2 
                        true_col_center = ncols / 2 
                        row_diff = irow - true_row_center 
                        col_diff = icol - true_col_center 
                        # 0s and 1s account for no true center in BWRs
                        idxy_1 = (np.array([int(true_col_center + col_diff + 0)]), np.array([int(true_row_center - row_diff - 1)]))  
                        idxy_2 = (np.array([int(true_col_center - row_diff - 1)]), np.array([int(true_row_center - col_diff - 1)])) 
                        idxy_3 = (np.array([int(true_col_center - col_diff - 1)]), np.array([int(true_row_center + row_diff + 0)]))
                                
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]

                        if core_map[irow,icol]: #skip empty locations
                            core_dict[core_map[irow,icol]] = dict_value

            if symmetry == 'octant':
                core_dict={}
                for irow in row_iter:
                    for icol in col_iter:
                        if irow == icol:
                            dict_value={'Symmetric_Assemblies':[], 'Value': None}
                            
                            # finds assembly location relative to center assembly
                            irow = int(irow); icol = int(icol)
                            true_row_center = nrows / 2 
                            true_col_center = ncols / 2 
                            row_diff = irow - true_row_center 
                            col_diff = icol - true_col_center 
                            # 0s and 1s account for no true center in BWRs
                            idxy_1 = (np.array([int(true_col_center + col_diff + 0)]), np.array([int(true_row_center - row_diff - 1)]))  
                            idxy_2 = (np.array([int(true_col_center - row_diff - 1)]), np.array([int(true_row_center - col_diff - 1)])) 
                            idxy_3 = (np.array([int(true_col_center - col_diff - 1)]), np.array([int(true_row_center + row_diff + 0)]))
                                    
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]

                            if core_map[irow,icol]: #skip empty locations
                                core_dict[core_map[irow,icol]] = dict_value
                        if irow > icol:
                            dict_value={'Symmetric_Assemblies':[], 'Value': None}

                            # finds assembly location relative to center assembly
                            irow = int(irow); icol = int(icol)
                            true_row_center = nrows / 2 
                            true_col_center = ncols / 2 
                            row_diff = irow - true_row_center 
                            col_diff = icol - true_col_center 
                            # 0s and 1s account for no true center in BWRs
                            idxy_1 = (np.array([int(true_col_center + col_diff + 0)]), np.array([int(true_row_center + row_diff + 0)]))
                            idxy_2 = (np.array([int(true_col_center - row_diff - 1)]), np.array([int(true_row_center + col_diff + 0)]))
                            idxy_3 = (np.array([int(true_col_center - col_diff - 1)]), np.array([int(true_row_center + row_diff + 0)]))
                            idxy_4 = (np.array([int(true_col_center - col_diff -1)]), np.array([int(true_row_center - row_diff -1 )]))
                            idxy_5 = (np.array([int(true_col_center - row_diff -1)]), np.array([int(true_row_center - col_diff -1 )]))
                            idxy_6 = (np.array([int(true_col_center + row_diff + 0)]), np.array([int(true_row_center - col_diff - 1)]))
                            idxy_7 = (np.array([int(true_col_center + col_diff - 0)]), np.array([int(true_row_center - row_diff - 1)]))
                                    
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0], core_map[idxy_4][0], core_map[idxy_5][0], core_map[idxy_6][0], core_map[idxy_7][0]]

                            if core_map[irow,icol]: #skip empty locations
                                core_dict[core_map[irow,icol]] = dict_value



        else: # find symetric locations for cores with odd number of rows and cols (PWRS)
            ## Link locations across symmetry zones (e.g. quadrants or octants)
            if symmetry == 'quarter':
                core_dict={}
                for irow in row_iter:
                    for icol in col_iter:
                        dict_value={'Symmetric_Assemblies':[],
                                    'Value': None}
                        if (irow,icol) == sym_center:
                            irow = int(irow); icol = int(icol)
                            pass
                        elif irow == sym_corner[0] and icol != sym_center[1]:  
                            continue
                        elif icol == sym_vertical[1] and irow != sym_center[0]:
                            irow = int(irow); icol = int(icol)
                            idy = core_id[irow,icol][0]
                            idx = core_id[irow,icol][1]
                            idxy_1 = np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                            idxy_2 = np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                            idxy_3 = np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                        else:
                            irow = int(irow); icol = int(icol)
                            idy = core_id[irow,icol][0]
                            idx = core_id[irow,icol][1]
                            idxy_1 = np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                            idxy_2 = np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                            idxy_3 = np.where((core_id[:,:,0] == -idx) & (core_id[:,:,1] == idy))

                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                        if core_map[irow,icol]: #skip empty locations
                            core_dict[core_map[irow,icol]] = dict_value
            elif symmetry == 'octant':
                core_dict={}
                for irow in row_iter:
                    for icol in col_iter:
                        if icol>irow:
                            continue
                        dict_value={'Symmetric_Assemblies':[],
                                    'Value': None}
                        if (irow,icol) == sym_center:
                            irow = int(irow); icol = int(icol)
                            pass
                        elif icol == sym_vertical[1] and irow != sym_center[0]:
                            irow = int(irow); icol = int(icol)
                            idy = core_id[irow,icol][0]
                            idx = core_id[irow,icol][1]
                            idxy_1= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                            idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                            idxy_3= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                        elif icol == irow and irow != sym_center[0]:
                            irow = int(irow); icol = int(icol)
                            idy = core_id[irow,icol][0]
                            idx = core_id[irow,icol][1]
                            idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                            idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                            idxy_3= np.where((core_id[:,:,0] == idy)  & (core_id[:,:,1] == -idx))
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                        else:
                            irow = int(irow); icol = int(icol)
                            idy = core_id[irow,icol][0]
                            idx = core_id[irow,icol][1]
                            idxy_1= np.where((core_id[:,:,0] == -idx) & (core_id[:,:,1] == -idy))
                            idxy_2= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                            idxy_3= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                            idxy_4= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                            idxy_5= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                            idxy_6= np.where((core_id[:,:,0] == -idx) & (core_id[:,:,1] == idy))
                            idxy_7= np.where((core_id[:,:,0] == idy) & (core_id[:,:,1] == -idx))
                            dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0], core_map[idxy_4][0],
                                                                core_map[idxy_5][0], core_map[idxy_6][0], core_map[idxy_7][0]]
                        if core_map[irow,icol]: #skip empty locations
                            core_dict[core_map[irow,icol]] = dict_value

        return core_dict


class Prepare_Problem_Values():
    """
    Class for preparing the necessary input values while setting up 
    each type of calculation.
    
    Written by Nicholas Rollins. 10/04/2024
    """
    def prepare_cycle(input_obj):
        """
        Prepares the provided input parameters in the necessary format for writing to
        the external model's input file.
        
        Updated by Nicholas Rollins. 10/04/2024
        """
    ## Parse cross sections and FA labels
        xs_list, tag_list = Problem_Preparation_Tools.parse_FAlabels(input_obj)
        
    ## Generate core map matching the length and symmetry of the chromosome
        core_shape = LWR_Core_Shapes.get_core_shape(input_obj.nrow, input_obj.ncol, input_obj.num_assemblies)
        core_dict = Problem_Preparation_Tools.generate_core(input_obj, core_shape) #create LP maps that match the length and symmetry of the chromosome.
 
    ## Generate map of core for printing to file
        if input_obj.map_size == "quarter":
            size_x = int(np.ceil(input_obj.nrow/2))
            size_y = int(np.ceil(input_obj.ncol/2))
            core_lattice = np.zeros((size_y,size_x), dtype='<U8')
            for y in range(size_y):
                for x in range(size_x):
                    val = core_shape[size_y-1+y,size_x-1+x]
                    if val is None:
                        val = "00"
                    core_lattice[y,x] = val
        else: #assume full map
            size_x = int(input_obj.nrow)
            size_y = int(input_obj.ncol)
            core_lattice = np.zeros((size_y,size_x), dtype='<U8')
            for y in range(size_y):
                for x in range(size_x):
                    val = core_shape[y,x]
                    if val is None:
                        val = "00"
                    core_lattice[y,x] = val

        pincal_loc = np.zeros((core_lattice.shape[0],core_lattice.shape[1]), dtype='<U20')
        for x in range(core_lattice.shape[0]):
            for y in range(core_lattice.shape[1]):
                loc = core_lattice[x,y]
                if loc[0] == "R":
                    pincal_loc[x,y]=0
                elif loc == "00":
                    if input_obj.map_size == 'quarter':
                        pincal_loc[x,y]=float('NaN')
                    else: #assume full geometry for printing
                        pincal_loc[x,y] = ' '
                else:
                    pincal_loc[x,y]=1
        
        # full-core map of locations; used for multi-cycle calcs.
        if input_obj.map_size == "quarter":
            size_x = int(np.ceil(input_obj.nrow/2))
            size_y = int(np.ceil(input_obj.ncol/2))
            full_core_locs = np.zeros((size_y,size_x), dtype='<U8')
            for y in range(size_y):
                for x in range(size_x):
                    val = core_shape[size_y-1+y,size_x-1+x]
                    if val is None or val[0] == "R":
                        full_core_locs[y,x] = "    "
                    else:
                        full_core_locs[y,x] = val[0] + '-' + val[1:]
        else: #assume full map
            size_x = int(input_obj.nrow)
            size_y = int(input_obj.ncol)
            full_core_locs = np.zeros((size_y,size_x), dtype='<U8')
            for y in range(size_y):
                for x in range(size_x):
                    val = core_shape[y,x]
                    if val is None:
                        full_core_locs[y,x] = "    "
                    elif val[0] == "R":
                        if input_obj.code_interface in ["parcs342","parcs343"]:
                            full_core_locs[y,x] = "   0"
                        else:
                            full_core_locs[y,x] = "    "
                    else:
                        full_core_locs[y,x] = val[0] + '-' + val[1:]
        
    ## Store results alongside input data
        input_obj.xs_list = xs_list
        input_obj.tag_list = tag_list
        input_obj.core_dict = core_dict
        input_obj.core_lattice = core_lattice
        input_obj.pincal_loc = pincal_loc
        input_obj.full_core_locs = full_core_locs
        
        return input_obj


class LWR_Core_Shapes():
    """
    All supported core shapes must be defined here to be properly understood within MIDAS.
    
    Written by Nicholas Rollins. 10/04/2024
    """
    def get_core_shape(num_rows, num_cols, num_FA):
        core_shape = {}
        core_shape[(17,17)] = {}
        core_shape[(17,17)][193] = [None , None , None  , None ,"R0004","R0005","R0006","R0007","R0008","R0009","R0010","R0011","R0012",  None , None ,  None ,  None ,
                                    None , None ,"R0102","R0103","R0104", "E01" , "F01" , "G01" , "H01" , "I01" , "J01" , "K01" ,"R0112","R0113","R0114",  None ,  None ,
                                    None ,"R0201","R0202", "C02" , "D02" , "E02" , "F02" , "G02" , "H02" , "I02" , "J02" , "K02" , "L02" , "M02" ,"R0214","R0215",  None ,
                                    None ,"R0301", "B03" , "C03" , "D03" , "E03" , "F03" , "G03" , "H03" , "I03" , "J03" , "K03" , "L03" , "M03" , "N03" ,"R0315",  None ,
                                   "R0400","R0401", "B04" , "C04" , "D04" , "E04" , "F04" , "G04" , "H04" , "I04" , "J04" , "K04" , "L04" , "M04" , "N04" ,"R0415","R0416",
                                   "R0500", "A05" , "B05" , "C05" , "D05" , "E05" , "F05" , "G05" , "H05" , "I05" , "J05" , "K05" , "L05" , "M05" , "N05" , "O05" ,"R0516",
                                   "R0600", "A06" , "B06" , "C06" , "D06" , "E06" , "F06" , "G06" , "H06" , "I06" , "J06" , "K06" , "L06" , "M06" , "N06" , "O06" ,"R0616",
                                   "R0700", "A07" , "B07" , "C07" , "D07" , "E07" , "F07" , "G07" , "H07" , "I07" , "J07" , "K07" , "L07" , "M07" , "N07" , "O07" ,"R0716",
                                   "R0800", "A08" , "B08" , "C08" , "D08" , "E08" , "F08" , "G08" , "H08" , "I08" , "J08" , "K08" , "L08" , "M08" , "N08" , "O08" ,"R0816",
                                   "R0900", "A09" , "B09" , "C09" , "D09" , "E09" , "F09" , "G09" , "H09" , "I09" , "J09" , "K09" , "L09" , "M09" , "N09" , "O09" ,"R0916",
                                   "R1000", "A10" , "B10" , "C10" , "D10" , "E10" , "F10" , "G10" , "H10" , "I10" , "J10" , "K10" , "L10" , "M10" , "N10" , "O10" ,"R1016",
                                   "R1100", "A11" , "B11" , "C11" , "D11" , "E11" , "F11" , "G11" , "H11" , "I11" , "J11" , "K11" , "L11" , "M11" , "N11" , "O11" ,"R1116",
                                   "R1200","R1201", "B12" , "C12" , "D12" , "E12" , "F12" , "G12" , "H12" , "I12" , "J12" , "K12" , "L12" , "M12" , "N12" ,"R1215","R1216",
                                    None ,"R1301", "B13" , "C13" , "D13" , "E13" , "F13" , "G13" , "H13" , "I13" , "J13" , "K13" , "L13" , "M13" , "N13" ,"R1315",  None ,
                                    None ,"R1401","R1402", "C14" , "D14" , "E14" , "F14" , "G14" , "H14" , "I14" , "J14" , "K14" , "L14" , "M14" ,"R1414","R1415",  None ,
                                    None , None ,"R1502","R1503","R1504", "E15" , "F15" , "G15" , "H15" , "I15" , "J15" , "K15" ,"R1512","R1513","R1514",  None ,  None ,
                                    None , None , None  , None ,"R1604","R1605","R1606","R1607","R1608","R1609","R1610","R1611","R1612",  None , None ,  None ,  None ]
        core_shape[(17,17)][157] = [None , None , None  , None , None , None ,"R0006","R0007","R0008","R0009","R0010", None  , None ,  None , None ,  None ,   None ,
                                    None , None , None  , None ,"R0104","R0105","R0106", "G01" , "H01" , "I01" ,"R0110","R0111","R0112", None , None ,  None ,  None ,
                                    None , None , None  ,"R0203","R0204", "E02" , "F02" , "G02" , "H02" , "I02" , "J02" , "K02" ,"R0212","R0213", None , None ,  None ,
                                    None , None ,"R0302","R0303", "D03" , "E03" , "F03" , "G03" , "H03" , "I03" , "J03" , "K03" , "L03" ,"R0313","R0314", None ,  None ,
                                    None ,"R0401","R0402", "C04" , "D04" , "E04" , "F04" , "G04" , "H04" , "I04" , "J04" , "K04" , "L04" , "M04" ,"R0414","R0415",  None ,
                                    None ,"R0501", "B05" , "C05" , "D05" , "E05" , "F05" , "G05" , "H05" , "I05" , "J05" , "K05" , "L05" , "M05" , "N05" ,"R0515",  None ,
                                   "R0600","R0601", "B06" , "C06" , "D06" , "E06" , "F06" , "G06" , "H06" , "I06" , "J06" , "K06" , "L06" , "M06" , "N06" ,"R0615","R0616",
                                   "R0700", "A07" , "B07" , "C07" , "D07" , "E07" , "F07" , "G07" , "H07" , "I07" , "J07" , "K07" , "L07" , "M07" , "N07" , "O07" ,"R0716",
                                   "R0800", "A08" , "B08" , "C08" , "D08" , "E08" , "F08" , "G08" , "H08" , "I08" , "J08" , "K08" , "L08" , "M08" , "N08" , "O08" ,"R0816",
                                   "R0900", "A09" , "B09" , "C09" , "D09" , "E09" , "F09" , "G09" , "H09" , "I09" , "J09" , "K09" , "L09" , "M09" , "N09" , "O09" ,"R0916",
                                   "R1000","R1001", "B10" , "C10" , "D10" , "E10" , "F10" , "G10" , "H10" , "I10" , "J10" , "K10" , "L10" , "M10" , "N10" ,"R1015","R1016",
                                    None ,"R1101", "B11" , "C11" , "D11" , "E11" , "F11" , "G11" , "H11" , "I11" , "J11" , "K11" , "L11" , "M11" , "N11" ,"R1115",  None ,
                                    None ,"R1201","R1202", "C12" , "D12" , "E12" , "F12" , "G12" , "H12" , "I12" , "J12" , "K12" , "L12" , "M12" ,"R1214","R1215",  None ,
                                    None , None ,"R1302","R1303", "D13" , "E13" , "F13" , "G13" , "H13" , "I13" , "J13" , "K13" , "L13" ,"R1313","R1314", None ,  None ,
                                    None , None , None  ,"R1403","R1404", "E14" , "F14" , "G14" , "H14" , "I14" , "J14" , "K14" ,"R1412","R1413", None , None ,  None ,
                                    None , None , None  , None ,"R1504","R1505","R1506", "G15" , "H15" , "I15" ,"R1510","R1511","R1512", None , None ,  None ,  None ,
                                    None , None , None  , None , None  , None ,"R1606","R1607","R1608","R1609","R1610", None  , None ,  None , None ,  None ,  None ]
        core_shape[(17,17)][157] = [None , None , None  , None , None , None ,"R0006","R0007","R0008","R0009","R0010", None  , None ,  None , None ,  None ,   None ,
                                    None , None , None  , None ,"R0104","R0105","R0106", "G01" , "H01" , "I01" ,"R0110","R0111","R0112", None , None ,  None ,  None ,
                                    None , None , None  ,"R0203","R0204", "E02" , "F02" , "G02" , "H02" , "I02" , "J02" , "K02" ,"R0212","R0213", None , None ,  None ,
                                    None , None ,"R0302","R0303", "D03" , "E03" , "F03" , "G03" , "H03" , "I03" , "J03" , "K03" , "L03" ,"R0313","R0314", None ,  None ,
                                    None ,"R0401","R0402", "C04" , "D04" , "E04" , "F04" , "G04" , "H04" , "I04" , "J04" , "K04" , "L04" , "M04" ,"R0414","R0415",  None ,
                                    None ,"R0501", "B05" , "C05" , "D05" , "E05" , "F05" , "G05" , "H05" , "I05" , "J05" , "K05" , "L05" , "M05" , "N05" ,"R0515",  None ,
                                   "R0600","R0601", "B06" , "C06" , "D06" , "E06" , "F06" , "G06" , "H06" , "I06" , "J06" , "K06" , "L06" , "M06" , "N06" ,"R0615","R0616",
                                   "R0700", "A07" , "B07" , "C07" , "D07" , "E07" , "F07" , "G07" , "H07" , "I07" , "J07" , "K07" , "L07" , "M07" , "N07" , "O07" ,"R0716",
                                   "R0800", "A08" , "B08" , "C08" , "D08" , "E08" , "F08" , "G08" , "H08" , "I08" , "J08" , "K08" , "L08" , "M08" , "N08" , "O08" ,"R0816",
                                   "R0900", "A09" , "B09" , "C09" , "D09" , "E09" , "F09" , "G09" , "H09" , "I09" , "J09" , "K09" , "L09" , "M09" , "N09" , "O09" ,"R0916",
                                   "R1000","R1001", "B10" , "C10" , "D10" , "E10" , "F10" , "G10" , "H10" , "I10" , "J10" , "K10" , "L10" , "M10" , "N10" ,"R1015","R1016",
                                    None ,"R1101", "B11" , "C11" , "D11" , "E11" , "F11" , "G11" , "H11" , "I11" , "J11" , "K11" , "L11" , "M11" , "N11" ,"R1115",  None ,
                                    None ,"R1201","R1202", "C12" , "D12" , "E12" , "F12" , "G12" , "H12" , "I12" , "J12" , "K12" , "L12" , "M12" ,"R1214","R1215",  None ,
                                    None , None ,"R1302","R1303", "D13" , "E13" , "F13" , "G13" , "H13" , "I13" , "J13" , "K13" , "L13" ,"R1313","R1314", None ,  None ,
                                    None , None , None  ,"R1403","R1404", "E14" , "F14" , "G14" , "H14" , "I14" , "J14" , "K14" ,"R1412","R1413", None , None ,  None ,
                                    None , None , None  , None ,"R1504","R1505","R1506", "G15" , "H15" , "I15" ,"R1510","R1511","R1512", None , None ,  None ,  None ,
                                    None , None , None  , None , None  , None ,"R1606","R1607","R1608","R1609","R1610", None  , None ,  None , None ,  None ,  None ]
    
        core_shape[(32,32)] = {}
        core_shape[(32,32)][764] = [ None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  , 'R0008', 'R0009', 'R0010', 'R0011', 'R0012', 'R0013', 'R0014', 'R0015', 'R0016', 'R0017', 'R0018', 'R0019', 'R0020', 'R0021', 'R0022', 'R0023',  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  , 'R0107', 'R0108',   'I01',   'J01',   'K01',   'L01',   'M01',   'N01',   'O01',   'P01',   'Q01',   'R01',   'S01',   'T01',   'U01',   'V01', 'R0123', 'R0124',  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  , 'R0205', 'R0206', 'R0207',   'H02',   'I02',   'J02',   'K02',   'L02',   'M02',   'N02',   'O02',   'P02',   'Q02',   'R02',   'S02',   'T02',   'U02',   'V02',   'W02', 'R0224', 'R0225', 'R0226',  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  , 'R0305',   'F03',   'G03',   'H03',   'I03',   'J03',   'K03',   'L03',   'M03',   'N03',   'O03',   'P03',   'Q03',   'R03',   'S03',   'T03',   'U03',   'V03',   'W03',   'X03',   'Y03', 'R0326',  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  , 'R0404', 'R0405',   'F04',   'G04',   'H04',   'I04',   'J04',   'K04',   'L04',   'M04',   'N04',   'O04',   'P04',   'Q04',   'R04',   'S04',   'T04',   'U04',   'V04',   'W04',   'X04',   'Y04', 'R0426', 'R0427',  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  , 'R0502', 'R0503', 'R0504',   'E05',   'F05',   'G05',   'H05',   'I05',   'J05',   'K05',   'L05',   'M05',   'N05',   'O05',   'P05',   'Q05',   'R05',   'S05',   'T05',   'U05',   'V05',   'W05',   'X05',   'Y05',   'Z05', 'R0527', 'R0528', 'R0529',  None  ,  None  ,
                                     None  ,  None  , 'R0602',   'C06',   'D06',   'E06',   'F06',   'G06',   'H06',   'I06',   'J06',   'K06',   'L06',   'M06',   'N06',   'O06',   'P06',   'Q06',   'R06',   'S06',   'T06',   'U06',   'V06',   'W06',   'X06',   'Y06',   'Z06',   'AA06', 'AB06', 'R0629',  None  ,  None  , 
                                     None  , 'R0701', 'R0702',   'C07',   'D07',   'E07',   'F07',   'G07',   'H07',   'I07',   'J07',   'K07',   'L07',   'M07',   'N07',   'O07',   'P07',   'Q07',   'R07',   'S07',   'T07',   'U07',   'V07',   'W07',   'X07',   'Y07',   'Z07',   'AA07', 'AB07', 'R0729', 'R0730',  None  ,
                                    'R0800', 'R0801',   'B08',   'C08',   'D08',   'E08',   'F08',   'G08',   'H08',   'I08',   'J08',   'K08',   'L08',   'M08',   'N08',   'O08',   'P08',   'Q08',   'R08',   'S08',   'T08',   'U08',   'V08',   'W08',   'X08',   'Y08',   'Z08',   'AA08', 'AB08',  'AC08', 'R0830', 'R0831',
                                    'R0900',   'A09',   'B09',   'C09',   'D09',   'E09',   'F09',   'G09',   'H09',   'I09',   'J09',   'K09',   'L09',   'M09',   'N09',   'O09',   'P09',   'Q09',   'R09',   'S09',   'T09',   'U09',   'V09',   'W09',   'X09',   'Y09',   'Z09',   'AA09', 'AB09',  'AC09',  'AD09', 'R0931',
                                    'R1000',   'A10',   'B10',   'C10',   'D10',   'E10',   'F10',   'G10',   'H10',   'I10',   'J10',   'K10',   'L10',   'M10',   'N10',   'O10',   'P10',   'Q10',   'R10',   'S10',   'T10',   'U10',   'V10',   'W10',   'X10',   'Y10',   'Z10',   'AA10', 'AB10',  'AC10',  'AD10', 'R1031',
                                    'R1100',   'A11',   'B11',   'C11',   'D11',   'E11',   'F11',   'G11',   'H11',   'I11',   'J11',   'K11',   'L11',   'M11',   'N11',   'O11',   'P11',   'Q11',   'R11',   'S11',   'T11',   'U11',   'V11',   'W11',   'X11',   'Y11',   'Z11',   'AA11', 'AB11',  'AC11',  'AD11', 'R1131',
                                    'R1200',   'A12',   'B12',   'C12',   'D12',   'E12',   'F12',   'G12',   'H12',   'I12',   'J12',   'K12',   'L12',   'M12',   'N12',   'O12',   'P12',   'Q12',   'R12',   'S12',   'T12',   'U12',   'V12',   'W12',   'X12',   'Y12',   'Z12',   'AA12', 'AB12',  'AC12',  'AD12', 'R1231',
                                    'R1300',   'A13',   'B13',   'C13',   'D13',   'E13',   'F13',   'G13',   'H13',   'I13',   'J13',   'K13',   'L13',   'M13',   'N13',   'O13',   'P13',   'Q13',   'R13',   'S13',   'T13',   'U13',   'V13',   'W13',   'X13',   'Y13',   'Z13',   'AA13', 'AB13',  'AC13',  'AD13', 'R1331',
                                    'R1400',   'A14',   'B14',   'C14',   'D14',   'E14',   'F14',   'G14',   'H14',   'I14',   'J14',   'K14',   'L14',   'M14',   'N14',   'O14',   'P14',   'Q14',   'R14',   'S14',   'T14',   'U14',   'V14',   'W14',   'X14',   'Y14',   'Z14',   'AA14', 'AB14',  'AC14',  'AD14', 'R1431',
                                    'R1500',   'A15',   'B15',   'C15',   'D15',   'E15',   'F15',   'G15',   'H15',   'I15',   'J15',   'K15',   'L15',   'M15',   'N15',   'O15',   'P15',   'Q15',   'R15',   'S15',   'T15',   'U15',   'V15',   'W15',   'X15',   'Y15',   'Z15',   'AA15', 'AB15',  'AC15',  'AD15', 'R1531',
                                    'R1600',   'A16',   'B16',   'C16',   'D16',   'E16',   'F16',   'G16',   'H16',   'I16',   'J16',   'K16',   'L16',   'M16',   'N16',   'O16',   'P16',   'Q16',   'R16',   'S16',   'T16',   'U16',   'V16',   'W16',   'X16',   'Y16',   'Z16',   'AA16', 'AB16',  'AC16',  'AD16', 'R1631',
                                    'R1700',   'A17',   'B17',   'C17',   'D17',   'E17',   'F17',   'G17',   'H17',   'I17',   'J17',   'K17',   'L17',   'M17',   'N17',   'O17',   'P17',   'Q17',   'R17',   'S17',   'T17',   'U17',   'V17',   'W17',   'X17',   'Y17',   'Z17',   'AA17', 'AB17',  'AC17',  'AD17', 'R1731',
                                    'R1800',   'A18',   'B18',   'C18',   'D18',   'E18',   'F18',   'G18',   'H18',   'I18',   'J18',   'K18',   'L18',   'M18',   'N18',   'O18',   'P18',   'Q18',   'R18',   'S18',   'T18',   'U18',   'V18',   'W18',   'X18',   'Y18',   'Z18',   'AA18', 'AB18',  'AC18',  'AD18', 'R1831',
                                    'R1900',   'A19',   'B19',   'C19',   'D19',   'E19',   'F19',   'G19',   'H19',   'I19',   'J19',   'K19',   'L19',   'M19',   'N19',   'O19',   'P19',   'Q19',   'R19',   'S19',   'T19',   'U19',   'V19',   'W19',   'X19',   'Y19',   'Z19',   'AA19', 'AB19',  'AC19',  'AD19', 'R1931',
                                    'R2000',   'A20',   'B20',   'C20',   'D20',   'E20',   'F20',   'G20',   'H20',   'I20',   'J20',   'K20',   'L20',   'M20',   'N20',   'O20',   'P20',   'Q20',   'R20',   'S20',   'T20',   'U20',   'V20',   'W20',   'X20',   'Y20',   'Z20',   'AA20', 'AB20',  'AC20',  'AD20', 'R2031',
                                    'R2100',   'A21',   'B21',   'C21',   'D21',   'E21',   'F21',   'G21',   'H21',   'I21',   'J21',   'K21',   'L21',   'M21',   'N21',   'O21',   'P21',   'Q21',   'R21',   'S21',   'T21',   'U21',   'V21',   'W21',   'X21',   'Y21',   'Z21',   'AA21', 'AB21',  'AC21',  'AD21', 'R2131',
                                    'R2200',   'A22',   'B22',   'C22',   'D22',   'E22',   'F22',   'G22',   'H22',   'I22',   'J22',   'K22',   'L22',   'M22',   'N22',   'O22',   'P22',   'Q22',   'R22',   'S22',   'T22',   'U22',   'V22',   'W22',   'X22',   'Y22',   'Z22',   'AA22', 'AB22',  'AC22',  'AD22', 'R2231',
                                    'R2300', 'R2301',   'B23',   'C23',   'D23',   'E23',   'F23',   'G23',   'H23',   'I23',   'J23',   'K23',   'L23',   'M23',   'N23',   'O23',   'P23',   'Q23',   'R23',   'S23',   'T23',   'U23',   'V23',   'W23',   'X23',   'Y23',   'Z23',   'AA23', 'AB23',  'AC23', 'R2330', 'R2331',
                                     None  , 'R2401', 'R2402',   'C24',   'D24',   'E24',   'F24',   'G24',   'H24',   'I24',   'J24',   'K24',   'L24',   'M24',   'N24',   'O24',   'P24',   'Q24',   'R24',   'S24',   'T24',   'U24',   'V24',   'W24',   'X24',   'Y24',   'Z24',   'AA24', 'AB24', 'R2429', 'R2430',  None  ,
                                     None  ,  None  , 'R2502',   'C25',   'D25',   'E25',   'F25',   'G25',   'H25',   'I25',   'J25',   'K25',   'L25',   'M25',   'N25',   'O25',   'P25',   'Q25',   'R25',   'S25',   'T25',   'U25',   'V25',   'W25',   'X25',   'Y25',   'Z25',   'AA25', 'AB25', 'R2529',  None  ,  None  ,
                                     None  ,  None  , 'R2602', 'R2603', 'R2604',   'E26',   'F26',   'G26',   'H26',   'I26',   'J26',   'K26',   'L26',   'M26',   'N26',   'O26',   'P26',   'Q26',   'R26',   'S26',   'T26',   'U26',   'V26',   'W26',   'X26',   'Y26',   'Z26', 'R2627', 'R2628', 'R2629',  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  , 'R2704', 'R2705',   'F27',   'G27',   'H27',   'I27',   'J27',   'K27',   'L27',   'M27',   'N27',   'O27',   'P27',   'Q27',   'R27',   'S27',   'T27',   'U27',   'V27',   'W27',   'X27',   'Y27', 'R2726', 'R2727',  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  , 'R2805',   'F28',   'G28',   'H28',   'I28',   'J28',   'K28',   'L28',   'M28',   'N28',   'O28',   'P28',   'Q28',   'R28',   'S28',   'T28',   'U28',   'V28',   'W28',   'X28',   'Y28', 'R2826',  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  , 'R2905', 'R2906', 'R2907',   'H29',   'I29',   'J29',   'K29',   'L29',   'M29',   'N29',   'O29',   'P29',   'Q29',   'R29',   'S29',   'T29',   'U29',   'V29',   'W29', 'R2924', 'R2925', 'R2926',  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  , 'R3007', 'R3008',   'I30',   'J30',   'K30',   'L30',   'M30',   'N30',   'O30',   'P30',   'Q30',   'R30',   'S30',   'T30',   'U30',   'V30', 'R3023', 'R3024',  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,
                                     None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  , 'R3108', 'R3109', 'R3110', 'R3111', 'R3112', 'R3113', 'R3114', 'R3115', 'R3116', 'R3117', 'R3118', 'R3119', 'R3120', 'R3121', 'R3122', 'R3123',  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None  ,  None ]


    
        return np.array(core_shape[(num_rows,num_cols)][num_FA]).reshape((num_rows,num_cols))
    
    def get_symmetry_multiplicity(num_rows, num_cols, num_FA, symmetry):
        multdict = {}
        multdict[(17,17)] = {}
        multdict[(17,17)][193] = {}
        multdict[(17,17)][193]['octant'] = { 0:1,\
                                             1:4,  2:4,\
                                             3:4,  4:8,  5:4,\
                                             6:4,  7:8,  8:8,  9:4,\
                                            10:4, 11:8, 12:8, 13:8, 14:4,\
                                            15:4, 16:8, 17:8, 18:8, 19:8, 20:4,\
                                            21:4, 22:8, 23:8, 24:8, 25:8, 26:8,\
                                            27:4, 28:8, 29:8, 30:8}
        multdict[(17,17)][193]['quarter'] = { 0:1,\
                                              1:4,  2:4,  3:4,  4:4,  5:4,  6:4,  7:4,  8:4,\
                                              9:4, 10:4, 11:4, 12:4, 13:4, 14:4, 15:4, 16:4,\
                                             17:4, 18:4, 19:4, 20:4, 21:4, 22:4, 23:4, 24:4,\
                                             25:4, 26:4, 27:4, 28:4, 29:4, 30:4, 31:4,\
                                             32:4, 33:4, 34:4, 35:4, 36:4, 37:4, 38:4,\
                                             39:4, 40:4, 41:4, 42:4, 43:4, 44:4,\
                                             45:4, 46:4, 47:4, 48:4}
        multdict[(17,17)][157] = {}
        multdict[(17,17)][157]['quarter'] = { 0:1,\
                                              1:4,  2:4,  3:4,  4:4,  5:4,  6:4,  7:4,  8:4,\
                                              9:4, 10:4, 11:4, 12:4, 13:4, 14:4, 15:4,\
                                             16:4, 17:4, 18:4, 19:4, 20:4, 21:4, 22:4,\
                                             23:4, 24:4, 25:4, 26:4, 27:4, 28:4,\
                                             29:4, 30:4, 31:4, 32:4, 33:4,\
                                             34:4, 35:4, 36:4, 37:4,\
                                             38:4, 39:4}
        multdict[(32,32)] = {}
        multdict[(32,32)][764] = {}
        multdict[(32,32)][764]['octant'] = { 0:4,\
                                             1:8, 2:4, \
                                             3:8, 4:8, 5:4, \
                                             6:8, 7:8, 8:8, 9:4, \
                                             10:8, 11:8, 12:8, 13:8, 14:4, \
                                             15:8, 16:8, 17:8, 18:8, 19:8, 20:4, \
                                             21:8, 22:8, 23:8, 24:8, 25:8, 26:8, 27:4, \
                                             28:8, 29:8, 30:8, 31:8, 32:8, 33:8, 34:8, 35:4,  \
                                             36:8, 37:8, 38:8, 39:8, 40:8, 41:8, 42:8, 43:8, 44:4, \
                                             45:8, 46:8, 47:8, 48:8, 49:8, 50:8, 51:8, 52:8, 53:8, 54:4, \
                                             55:8, 56:8, 57:8, 58:8, 59:8, 60:8, 61:8, 62:8, 63:8, 64:8, 65:4, \
                                             66:8, 67:8, 68:8, 69:8, 70:8, 71:8, 72:8, 73:8, 74:8, 75:8, \
                                             76:8, 77:8, 78:8, 79:8, 80:8, 81:8, 82:8, 83:8, 84:8, 85:8, \
                                             86:8, 87:8, 88:8, 89:8, 90:8, 91:8, 92:8, 93:8, \
                                             94:8, 95:8, 96:8, 97:8, 98:8, 99:8, 100:8} 

        multdict[(32,32)][764]['quarter'] = { 0:4,  1:4,  2:4,  3:4,  4:4,  5:4,  6:4,  7:4,  8:4,  9:4, 10:4, 11:4, 12:4, 13:4, 14:4,\
                                             15:4, 16:4, 17:4, 18:4, 19:4, 20:4, 21:4, 22:4, 23:4, 24:4, 25:4, 26:4, 27:4, 28:4, 29:4,\
                                             30:4, 31:4, 32:4, 33:4, 34:4, 35:4, 36:4, 37:4, 38:4, 39:4, 40:4, 41:4, 42:4, 43:4, 44:4,\
                                             45:4, 46:4, 47:4, 48:4, 49:4, 50:4, 51:4, 52:4, 53:4, 54:4, 55:4, 56:4, 57:4, 58:4, 59:4,\
                                             60:4, 61:4, 62:4, 63:4, 64:4, 65:4, 66:4, 67:4, 68:4, 69:4, 70:4, 71:4, 72:4, 73:4, 74:4,\
                                             75:4, 76:4, 77:4, 78:4, 79:4, 80:4, 81:4, 82:4, 83:4, 84:4, 85:4, 86:4, 87:4, 88:4, 89:4,\
                                             90:4, 91:4, 92:4, 93:4, 94:4, 95:4, 96:4, 97:4, 98:4, 99:4, 100:4, 101:4, 102:4, 103:4, 104:4,\
                                             105:4, 106:4, 107:4, 108:4, 109:4, 110:4, 111:4, 112:4, 113:4, 114:4, 115:4, 116:4, 117:4, 118:4, \
                                             119:4, 120:4, 121:4, 122:4, 123:4, 124:4, 125:4, 126:4, 127:4, 128:4, 129:4, 130:4, 131:4, \
                                             132:4, 133:4, 134:4, 135:4, 136:4, 137:4, 138:4, 139:4, 140:4, 141:4, 142:4, 143:4, 144:4, \
                                             145:4, 146:4, 147:4, 148:4, 149:4, 150:4, 151:4, 152:4, 153:4, 154:4, 155:4, \
                                             156:4, 157:4, 158:4, 159:4, 160:4, 161:4, 162:4, 163:4, 164:4, 165:4, \
                                             166:4, 167:4, 168:4, 169:4, 170:4, 171:4, 172:4, 173:4, 174:4, 175:4, \
                                             176:4, 177:4, 178:4, 179:4, 180:4, 181:4, 182:4, 183:4, \
                                             184:4, 185:4, 186:4, 187:4, 188:4, 189:4, 190:4} 

        return multdict[(num_rows,num_cols)][num_FA][symmetry]
    
    def count_in_LP(multdict, chromosome):
        gene_type_count = {'total':0}
        for i in range(len(chromosome)):
            multiplicity = multdict[i]
            gene = chromosome[i]
            if gene in gene_type_count:
                gene_type_count[gene] += multiplicity
            else:
                gene_type_count[gene] = multiplicity
            gene_type_count['total'] += multiplicity
        return gene_type_count
        
        