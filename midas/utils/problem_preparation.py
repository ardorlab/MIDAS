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
        xs_list = {'fuel':[], 'top':[], 'radial':[], 'bot':[]}
        tag_list = {'fuel':[], 'reflectors':[], 'blankets':[]}
        
        for param in input_obj['assembly_options']['fuel']:
            xs_list['fuel'].append(param['serial'])
            tag_list['fuel'].append(param['type'])
        for param in input_obj['assembly_options']['reflectors']:
            if param['refl_type'] == 'all':
                xs_list['reflectors']['top'].append(param['serial'])
                xs_list['reflectors']['radial'].append(param['serial'])
                xs_list['reflectors']['bot'].append(param['serial'])
                tag_list['reflectors']['top'].append(10)
                tag_list['reflectors']['radial'].append(10)
                tag_list['reflectors']['bot'].append(10)
            elif param['refl_type'] == 'top':
                xs_list['reflectors']['top'].append(param['serial'])
                tag_list['reflectors']['top'].append(11)
            elif param['refl_type'] == 'radial':
                xs_list['reflectors']['radial'].append(param['serial'])
                tag_list['reflectors']['radial'].append(10)
            elif param['refl_type'] == 'bot':
                xs_list['reflectors']['bot'].append(param['serial'])
                tag_list['reflectors']['bot'].append(12)
        
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
        
        core_sym_map, fuel_sym_map = Problem_Preparation_Tools.symmetric_core(inp_obj.symmetry,nrows,ncols,core_map,core_id)
        
        return core_sym_map, fuel_sym_map

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
        sym_center   = (np.floor(nrow/2)+1,np.floor(ncol/2)+1) #x,y-coordinates for centerpoint of symmetry
        sym_vertical = (nrow-1,np.floor(ncol/2)+1)
        if symmetry == 'quarter':
            sym_corner   = (np.floor(nrow/2)+1,ncol-1)
        elif symmetry == 'octant':
            sym_corner   = (nrow-1,ncol-1)
        
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
                        pass
                    elif irow == sym_horizontal[0] and icol != sym_center[1]:                 
                        idy = core_id[irow,icol][0]
                        idx = core_id[irow,icol][1]
                        idxy_1= np.where((core_id[:,:,0] == idy) & (core_id[:,:,1] == -idx))
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0]]
                    elif icol == sym_vertical[1] and irow != sym_center[0]:
                        idy = core_id[irow,icol][0]
                        idx = core_id[irow,icol][1]
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0]]
                    else:
                        idy = core_id[irow,icol][0]
                        idx = core_id[irow,icol][1]
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                        idxy_3= np.where((core_id[:,:,0] == idy) & (core_id[:,:,1] == -idx))
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
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
                        pass
                    elif icol == sym_vertical[1] and irow != sym_center[0]:                 
                        idy = core_id[irow,icol][0]
                        idx = core_id[irow,icol][1]
                        idxy_1= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                        idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_3= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                    elif icol == irow and irow != sym_center[0]:
                        idy = core_id[irow,icol][0]
                        idx = core_id[irow,icol][1]
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                        idxy_3= np.where((core_id[:,:,0] == idy)  & (core_id[:,:,1] == -idx))
                        dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                    else:
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
                    core_dict[core_map[irow,icol]] = dict_value
        
        fuel_dict={} # extract fuel assembly locations
        for key, value in core_dict.items():
            if key is None:
                continue
            elif key[0]=="R":
                continue
            else:
                fuel_dict[key]=value
        
        return(core_dict,fuel_dict)


class Prepare_Problem_Values():
    """
    Class for preparing the necessary input values while setting up 
    each type of calculation.
    
    Written by Nicholas Rollins. 10/04/2024
    """
    def single_cycle_preparation(input_obj):
        """
        #!TODO: write docstring.
        
        Updated by Nicholas Rollins. 10/04/2024
        """
    ## Generate core map matching the length and symmetry of the chromosome
        core_shape = LWR_Core_Shapes.get_core_shape(input_obj.nrow, input_obj.ncol)
        core_dict = {}
        core_dict['core'], core_dict['fuel'] = Problem_Preparation_Tools.generate_core(input_obj, core_shape) #create LP maps that match the length and symmetry of the chromosome.
        
    ## Expand to full core map while maintaining symmetry
        full_core  = {}
        for key, value in core_dict['fuel'].items():
            full_core[key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core[skey]=value['Value'] 
        
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
                    #!core_lattice[x,y] = core_dict['Reflectors']['Radial']['Tag']
                    pincal_loc[x,y]=0
                elif loc == "00":
                    if input_obj.map_size == 'quarter':
                        pincal_loc[x,y]=float('NaN')
                    else: #assume full geometry for printing
                        pincal_loc[x,y] = ' '
                else:
                    pincal_loc[x,y]=1
        
    ## Store results alongside input data
        inp_obj.core_dict = core_dict
        inp_obj.full_core = full_core
        inp_obj.core_lattice = core_lattice
        inp_obj.pincal_loc = pincal_loc
        
        return inp_obj


class LWR_Core_Shapes():
    """
    All supported core shapes must be defined here to be properly understood within MIDAS.
    
    Written by Nicholas Rollins. 10/04/2024
    """
    def get_core_shape(num_rows, num_cols):
        core_shape = {}
        core_shape[(17,17)] = [None , None , None  , None ,"R0004","R0005","R0006","R0007","R0008","R0009","R0010","R0011","R0012",  None , None ,  None ,   None  ,
                               None , None ,"R0102","R0103","R0104", "E01" , "F01" , "G01" , "H01" , "I01" , "J01" , "K01" ,"R0112","R0113","R0114",  None ,   None  ,
                               None ,"R0201","R0202", "C02" , "D02" , "E02" , "F02" , "G02" , "H02" , "I02" , "J02" , "K02" , "L02" , "M02" ,"R0214","R0215",   None  ,
                               None ,"R0301", "B03" , "C03" , "D03" , "E03" , "F03" , "G03" , "H03" , "I03" , "J03" , "K03" , "L03" , "M03" , "N03" ,"R0315",   None  ,
                              "R0400","R0401", "B04" , "C04" , "D04" , "E04" , "F04" , "G04" , "H04" , "I04" , "J04" , "K04" , "L04" , "M04" , "N04" ,"R0415", "R0416" ,
                              "R0500", "A05" , "B05" , "C05" , "D05" , "E05" , "F05" , "G05" , "H05" , "I05" , "J05" , "K05" , "L05" , "M05" , "N05" , "O05" , "R0516" ,
                              "R0600", "A06" , "B06" , "C06" , "D06" , "E06" , "F06" , "G06" , "H06" , "I06" , "J06" , "K06" , "L06" , "M06" , "N06" , "O06" , "R0616" ,
                              "R0700", "A07" , "B07" , "C07" , "D07" , "E07" , "F07" , "G07" , "H07" , "I07" , "J07" , "K07" , "L07" , "M07" , "N07" , "O07" , "R0716" ,
                              "R0800", "A08" , "B08" , "C08" , "D08" , "E08" , "F08" , "G08" , "H08" , "I08" , "J08" , "K08" , "L08" , "M08" , "N08" , "O08" , "R0816" ,
                              "R0900", "A09" , "B09" , "C09" , "D09" , "E09" , "F09" , "G09" , "H09" , "I09" , "J09" , "K09" , "L09" , "M09" , "N09" , "O09" , "R0916" ,
                              "R1000", "A10" , "B10" , "C10" , "D10" , "E10" , "F10" , "G10" , "H10" , "I10" , "J10" , "K10" , "L10" , "M10" , "N10" , "O10" , "R1016" ,
                              "R1100", "A11" , "B11" , "C11" , "D11" , "E11" , "F11" , "G11" , "H11" , "I11" , "J11" , "K11" , "L11" , "M11" , "N11" , "O11" , "R1116" ,
                              "R1200","R1201", "B12" , "C12" , "D12" , "E12" , "F12" , "G12" , "H12" , "I12" , "J12" , "K12" , "L12" , "M12" , "N12" ,"R1215", "R1216" ,
                               None ,"R1301", "B13" , "C13" , "D13" , "E13" , "F13" , "G13" , "H13" , "I13" , "J13" , "K13" , "L13" , "M13" , "N13" ,"R1315",   None  ,
                               None ,"R1401","R1402", "C14" , "D14" , "E14" , "F14" , "G14" , "H14" , "I14" , "J14" , "K14" , "L14" , "M14" ,"R1414","R1415",   None  ,
                               None , None ,"R1502","R1503","R1504", "E15" , "F15" , "G15" , "H15" , "I15" , "J15" , "K15" ,"R1512","R1513","R1514",  None ,   None  ,
                               None , None , None  , None ,"R1604","R1605","R1606","R1607","R1608","R1609","R1610","R1611","R1612",  None , None ,  None ,   None  ]
        
        return np.array(core_shape[(num_rows,num_cols)]).reshape((num_rows,num_cols))