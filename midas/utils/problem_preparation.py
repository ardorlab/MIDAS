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
                core_id.append((i-8,j-8))
        core_id = np.array(core_id).reshape((nrows,ncols,2))
        
        core_sym_map = Problem_Preparation_Tools.symmetric_core(input_obj.symmetry,nrows,ncols,core_map,core_id)
        
        return core_sym_map

    def symmetric_core(symmetry,nrows,ncols,core_map,core_id):
        """
        Generates the symmetric core map.

        Parameters: 
           - core_map: an e.g. 17x17 numpy array with the fuel assembly location names.
           - core_id: an e.g. 17x17x2 numpy array with coordinate indices for each fuel assembly location
           ranging from e.g. -8 to +8.

        Written by Gregory Delipei 7/12/2022
        Updated by Nicholas Rollins. 10/04/2024
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
        
        