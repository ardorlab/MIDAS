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
from pathlib import Path
import shutil 
from matplotlib import pyplot as plt
import subprocess
from subprocess import PIPE, STDOUT
from midas.utils.solution_types import Solution
from midas.utils.lcoe import LCOE, LCOE_MCYC
import xgboost as xgb

## Trial decorator ##
def trials(func):
    def wrapper(*args, **kwargs):
        MAX_TRIALS = 1000
        for i in range(MAX_TRIALS):
            try:
                return(func(*args, **kwargs))
            except Exception as e:
                if i == MAX_TRIALS - 1:
                    return(None)
                print(f"Trial {i+1} failed. Retrying...")
    return(wrapper)

class MCycle_Inventory_Loading_Pattern_Solution(Solution):
    """
    Solution class for designing multi-cycle loading patterns using PARCS code.
    Each cycle is run as a single-cycle.
    The reloaded fuel assemblies are grouped in pre-defined groups by reactivity.
    The inventory of reloaded assemblies includes assemblies from all cycles.
    Design limits are imposed on the maximum reload burnup and fresh feed.

    Parameters: None

    Written by Gregory Delipei. 02/12/2024
    """
    def __init__(self):
        Solution.__init__(self)
        self.type        = None
        self.number_pins = None
        self.model       = None
        self.symmetry    = None
        self.nrow=17
        self.ncol=17
        self.action_lower=-1
        self.action_upper=1
        self.core_dict={}
       
    def add_additional_information(self,settings):
        """
        Adds information on reactor parameters.

        Parameters
            settings: The settings dictionary for the parameters.
        
        Written by Gregory Delipei. 01/08/2022
        """
        self.symmetry=settings["genome"]["parcs_data"]["symmetry"]
        self.core_dict['core_map'], self.core_dict['core_id'] = self.generate_core()
        if 'inventory' in list(settings["genome"].keys()):
            self.core_dict['Inventory'] = settings["genome"]['inventory']
        else:
            inventory = {}
            for key,value in settings["genome"]["chromosomes"].items():
                inv_it = {}
                inv_it['Max_Limit']=np.inf
                inv_it['In_Design']=0
                inv_it['Cost']=250
                inv_it['Tag']=str(value['type'])
                inv_it['Cross_Section']=value['serial']
                inventory[key]=inv_it
            self.core_dict['Inventory'] = inventory
       
        info = settings['genome']['parcs_data']
        if 'xs_library' in info:
            self.library = info['xs_library']
        if 'power' in info:
            self.power = float(info['power'])
        if 'flow' in info:
            self.flow = float(info['flow'])
        if 'inlet_temperature' in info:
            self.inlet_temperature = float(info['inlet_temperature'])
        if 'map_size' in info:
            self.map_size= info['map_size']
        if 'number_assemblies' in info:
            self.number_assemblies = int(info['number_assemblies'])
        if 'ncycles' in info:
            self.ncycles = int(info['ncycles'])
            if 'design_limits' in settings['genome']:
                if 'optimize' in settings['genome']['design_limits']:
                    if settings['genome']['design_limits']['optimize']=='single_cycle':
                        self.ncycles=1
        if 'fixed_problem' == settings['optimization']['reproducer']:
            self.fixed_genome = True
        elif 'unique_genes' == settings['optimization']['reproducer']:
            self.fixed_genome = True
        self.settings = settings

    def generate_core(self):
        """
        Generates the 17x17 core map with consistent identifiers and treatment of symmetry.

        Parameters: None
        Additional comments:
          - The core_map is manually defined as a list line by line starting from the top of the core.
          It is advised to use the following naming conention for reflector assemblies 'ABCDE' 
          with A being R, BC indicating the row number (00-17) and DE the column number (00-17). For 
          the fuel assemblies it is advised to use the following naming convention 'ABC' with A being the 
          row letter indetifier (A-O) and BC being the column number (00-17). 

        Written by Gregory Delipei 7/12/2022
        """
        core_map = [  None ,  None ,  None ,  None ,"R0004","R0005","R0006","R0007","R0008","R0009","R0010","R0011","R0012",  None ,  None ,  None ,  None ,
                      None ,  None ,"R0102","R0103","R0104", "A05" , "A06" , "A07" , "A08" , "A09" , "A10" , "A11" ,"R0112","R0113","R0114",  None ,  None ,
                      None ,"R0201","R0202", "B03" , "B04" , "B05" , "B06" , "B07" , "B08" , "B09" , "B10" , "B11" , "B12" , "B13" ,"R0214","R0215",  None ,
                      None ,"R0301", "C02" , "C03" , "C04" , "C05" , "C06" , "C07" , "C08" , "C09" , "C10" , "C11" , "C12" , "C13" , "C14" ,"R0315",  None ,
                    "R0400","R0401", "D02" , "D03" , "D04" , "D05" , "D06" , "D07" , "D08" , "D09" , "D10" , "D11" , "D12" , "D13" , "D14" ,"R0415","R0416",
                    "R0500", "E01" , "E02" , "E03" , "E04" , "E05" , "E06" , "E07" , "E08" , "E09" , "E10" , "E11" , "E12" , "E13" , "E14" , "E15" ,"R0516",
                    "R0600", "F01" , "F02" , "F03" , "F04" , "F05" , "F06" , "F07" , "F08" , "F09" , "F10" , "F11" , "F12" , "F13" , "F14" , "F15" ,"R0616",
                    "R0700", "G01" , "G02" , "G03" , "G04" , "G05" , "G06" , "G07" , "G08" , "G09" , "G10" , "G11" , "G12" , "G13" , "G14" , "G15" ,"R0716",
                    "R0800", "H01" , "H02" , "H03" , "H04" , "H05" , "H06" , "H07" , "H08" , "H09" , "H10" , "H11" , "H12" , "H13" , "H14" , "H15" ,"R0816",
                    "R0900", "I01" , "I02" , "I03" , "I04" , "I05" , "I06" , "I07" , "I08" , "I09" , "I10" , "I11" , "I12" , "I13" , "I14" , "I15" ,"R0916",
                    "R1000", "J01" , "J02" , "J03" , "J04" , "J05" , "J06" , "J07" , "J08" , "J09" , "J10" , "J11" , "J12" , "J13" , "J14" , "J15" ,"R1016",
                    "R1100", "K01" , "K02" , "K03" , "K04" , "K05" , "K06" , "K07" , "K08" , "K09" , "K10" , "K11" , "K12" , "K13" , "K14" , "K15" ,"R1116",
                    "R1200","R1201", "L02" , "L03" , "L04" , "L05" , "L06" , "L07" , "L08" , "L09" , "L10" , "L11" , "L12" , "L13" , "L14" ,"R1215","R1216",
                      None ,"R1301", "M02" , "M03" , "M04" , "M05" , "M06" , "M07" , "M08" , "M09" , "M10" , "M11" , "M12" , "M13" , "M14" ,"R1315",  None ,
                      None ,"R1401","R1402", "N03" , "N04" , "N05" , "N06" , "N07" , "N08" , "N09" , "N10" , "N11" , "N12" , "N13" ,"R1414","R1415",  None ,
                      None ,  None ,"R1502","R1503","R1504", "O05" , "O06" , "O07" , "O08" , "O09" , "O10" , "O11" ,"R1512","R1513","R1514",  None ,  None ,
                      None ,  None ,  None ,  None ,"R1604","R1605","R1606","R1607","R1608","R1609","R1610","R1611","R1612",  None ,  None ,  None ,  None ]
        core_map = np.array(core_map).reshape((self.nrow,self.ncol))
        core_id = []
        for i in range(self.nrow-1,-1,-1):
            for j in range(self.ncol):
                core_id.append((i-8,j-8))
        core_id=np.array(core_id).reshape((self.nrow,self.ncol,2))
        if self.symmetry[0] == 'quarter':
            self.symmetry_axes = ((8,8),(8,16),(16,8))
            core, fuel = self.quarter_core(core_map,core_id,sym_type=self.symmetry[1])
            self.core_dict['core']= {
                'C1': copy.deepcopy(core),
                'C2':copy.deepcopy(core),
                'C3':copy.deepcopy(core)
            }
            self.core_dict['fuel']= {
                'C1': copy.deepcopy(fuel),
                'C2':copy.deepcopy(fuel),
                'C3':copy.deepcopy(fuel)
            }
        elif self.symmetry[0] == 'octant':
            self.symmetry_axes = ((8,8),(16,16),(16,8))
            core, fuel = self.octant_core(core_map,core_id)
            self.core_dict['core']= {
                'C1': copy.deepcopy(core),
                'C2':copy.deepcopy(core),
                'C3':copy.deepcopy(core)
            }
            self.core_dict['fuel']= {
                'C1': copy.deepcopy(fuel),
                'C2':copy.deepcopy(fuel),
                'C3':copy.deepcopy(fuel)
            }
        else:
            raise ValueError(
                f"The selected symmetry ({self.symmetry[0]}) is not valid."
            )
        return(core_map, core_id)

    def get_full_core(self):
        """
        Generates the 17x17 full fuel core from symmetry.

        Parameters: None
    
        Written by Gregory Delipei 7/24/2022
        """
        full_core  = {}
        full_core['C1']  = {}
        full_core['C2']  = {}
        full_core['C3']  = {}
        for key, value in self.core_dict['fuel']['C1'].items():
            full_core['C1'][key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core['C1'][skey]=value['Value'] 
        for key, value in self.core_dict['fuel']['C2'].items():
            full_core['C2'][key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core['C2'][skey]=value['Value'] 
        for key, value in self.core_dict['fuel']['C3'].items():
            full_core['C3'][key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core['C3'][skey]=value['Value'] 
        return(full_core)

    def quarter_core(self,core_map,core_id, sym_type='mirror'):
        """
        Generates the quarter core symmetry map.

        Parameters: 
           - core_map: a 17x17 numpy array with the fuel assembly location names.
           - core_id: a 17x17x2 numpy array with coordinate indices for each fuel assembly location
           ranging from -8 to +8.

        Written by Gregory Delipei 7/12/2022
        """
        sym_center = self.symmetry_axes[0]
        sym_horizontal = self.symmetry_axes[1]
        sym_vertical = self.symmetry_axes[2]
        if sym_vertical[0] > sym_center[0]:
            row_iter = np.arange(sym_center[0],sym_vertical[0]+1,1)
        else:
            row_iter = np.arange(sym_center[0],sym_vertical[0]-1,-1)
        if sym_horizontal[1] > sym_center[1]:
            col_iter = np.arange(sym_center[1],sym_horizontal[1]+1,1)
        else:
            col_iter = np.arange(sym_center[1],sym_horizontal[1]-1,-1)

        core_dict={}
        
        # Center assembly
        dict_value={'Symmetric_Assemblies':[],
                            'Value': None}
        irow = row_iter[0]
        icol = col_iter[0]
        core_dict[core_map[irow,icol]] = dict_value
        
        for irow in row_iter[1:]:
            for icol in col_iter:
                dict_value={'Symmetric_Assemblies':[],
                            'Value': None}
                if icol == sym_vertical[1]:                 
                    idy = core_id[irow,icol][0]
                    idx = core_id[irow,icol][1]
                    if sym_type == 'mirror':
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_2= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                        idxy_3= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                    elif sym_type == 'rotational':
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_2= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == idy))
                        idxy_3= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                    else:
                        raise ValueError(f"The selected symmetry ({sym_type}) is not valid.")
                    dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                else:
                    idy = core_id[irow,icol][0]
                    idx = core_id[irow,icol][1]
                    if sym_type == 'mirror':
                        idxy_1= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == idx))
                        idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                        idxy_3= np.where((core_id[:,:,0] == idy) & (core_id[:,:,1] == -idx))
                    elif sym_type == 'rotational':
                        idxy_1= np.where((core_id[:,:,0] == idx) & (core_id[:,:,1] == -idy))
                        idxy_2= np.where((core_id[:,:,0] == -idy) & (core_id[:,:,1] == -idx))
                        idxy_3= np.where((core_id[:,:,0] == -idx) & (core_id[:,:,1] == idy))
                    else:
                        raise ValueError(f"The selected symmetry ({sym_type}) is not valid.")
                    dict_value['Symmetric_Assemblies'] = [core_map[idxy_1][0], core_map[idxy_2][0], core_map[idxy_3][0]]
                core_dict[core_map[irow,icol]] = dict_value
        fuel_dict = self.extract_fuel(core_dict)
        return(core_dict,fuel_dict)
        
    def octant_core(self, core_map, core_id):
        """
        Generates the octant core symmetry map.

        Parameters: 
           - core_map: a 17x17 numpy array with the fuel assembly location names.
           - core_id: a 17x17x2 numpy array with coordinate indices for each fuel assembly location
           ranging from -8 to +8.

        Written by Gregory Delipei 7/12/2022
        """
        sym_center = self.symmetry_axes[0]
        sym_corner = self.symmetry_axes[1]
        sym_vertical = self.symmetry_axes[2]
        if sym_corner[0] > sym_center[0]:
            row_iter = np.arange(sym_center[0],sym_corner[0]+1,1)
        else:
            row_iter = np.arange(sym_center[0],sym_corner[0]-1,-1)
        if sym_corner[1] > sym_center[1]:
            col_iter = np.arange(sym_center[1],sym_corner[1]+1,1)
        else:
            col_iter = np.arange(sym_center[1],sym_corner[1]-1,-1)

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
                    idxy_3= np.where((core_id[:,:,0] == idy) & (core_id[:,:,1] == -idx))
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
                fuel_dict = self.extract_fuel(core_dict)
        return(core_dict,fuel_dict)

    def extract_fuel(self, core_dict):
        """
        Extracts the unique fuel assemblies from the core in a dictionary.

        Parameters: 
           - core_dict: a core dictionary including all the geometrical parameters.

        Written by Gregory Delipei 7/12/2022
        """
        fuel_dict={}
        for key, value in core_dict.items():
            if key is None:
                continue
            elif key[0]=="R":
                continue
            else:
                fuel_dict[key]=value
        return(fuel_dict)
 
    def plot_design(self,filepath):
        """
        Plot current loading pattern design.

        Parameters: None

        Written by Gregory Delipei 7/13/2022
        """
        color_fill = ['maroon','red','darkorange','limegreen','turquoise','pink','purple','plum','teal']
        nfa = 1
        tags=[]
        for key in self.core_dict['Inventory'].keys():
            if self.core_dict['Inventory'][key]['In_Design'] > 0:
                nfa+=1
                tags.append(key)
        ass_pitch = 21.21      

        plt.rcParams.update({'font.size': 4})
        fig=plt.figure()
        plt.axes()

        value_map = copy.deepcopy(self.core_dict['core_map'])
        value_map = value_map.astype('<U12')
        for key, value in self.core_dict['core'].items():
            value_map[np.where(self.core_dict['core_map']==key)] = value['Value']
            for isym in value['Symmetric_Assemblies']:
                value_map[np.where(self.core_dict['core_map']==isym)] = value['Value']

        for i in range(17):
            for j in range(17):
                ass_type=value_map[i,j]
                yloc=self.core_dict['core_id'][i,j][0]
                xloc=self.core_dict['core_id'][i,j][1]
                yid = yloc*ass_pitch
                xid = xloc*ass_pitch
                if ass_type == 'None':
                    continue
                else:
                    ass_color = color_fill[tags.index(ass_type)]
                    rectangle = plt.Rectangle((xid,yid),ass_pitch, ass_pitch,fc=ass_color,ec='grey')
                    plt.gca().add_patch(rectangle)
                    plt.text(xid+5,yid+9,self.core_dict['core_map'][i,j])


        plt.rcParams.update({'font.size': 8})
        rectangle = plt.Rectangle((12*ass_pitch,8*ass_pitch),0.3*ass_pitch,0.2*ass_pitch,1,fc="White")
        plt.gca().add_patch(rectangle)
        for i in range(len(tags)):    
            rectangle = plt.Rectangle((9.0*ass_pitch,8*ass_pitch-i*ass_pitch),0.4*ass_pitch,0.3*ass_pitch,1,fc=color_fill[i])
            plt.gca().add_patch(rectangle)
            plt.text(9.5*ass_pitch,8*ass_pitch-i*ass_pitch-0.03*ass_pitch,tags[i])

        plt.axis('scaled')
        plt.tick_params(left = False, right = False , labelleft = False ,
                        labelbottom = False, bottom = False)

        plt.savefig(filepath,bbox_inches='tight',dpi=300) 
        plt.close(fig)

    def genes_in_group(self,chromosome_map,group_name):
        """
        Returns a list of the genes in the chosen group
        """
        gene_list = []
        for gene in chromosome_map:
            if gene == 'symmetry_list':
                pass
            else:
                if group_name == chromosome_map[gene]['gene_group']:
                    gene_list.append(gene)

        return gene_list

    def is_gene_ok(self,chromosome_map,gene,space):
        """
        Checks if the gene is allowed in the desired location
        """
        gene_is_ok = True
        if not chromosome_map[gene]['map'][space]:
            gene_is_ok = False
        if space in chromosome_map['symmetry_list']:
            if self.my_group[chromosome_map[gene]['gene_group']] <= 1:
                gene_is_ok = False
        else:
            if not self.my_group[chromosome_map[gene]['gene_group']]:
                gene_is_ok = False
        if 'unique' in chromosome_map[gene]:
            if gene in self.genome:
                gene_is_ok = False

        return gene_is_ok

    def feed_counter(self,genome):
        """
        Computes the fresh feed for a specific genome

        Parameters:
            genome: a loading pattern corresponding to a single cycle
        """
        burnt_fuel_list = []
        fresh_fuel_list = []
        chromosomes = self.settings['genome']['chromosomes']
        for key in chromosomes.keys():
            if chromosomes[key]['type'] =='reload':
                burnt_fuel_list.append(key) 
            else:
                fresh_fuel_list.append(key) 
        feed_count = 0
        loc_tag = list(self.core_dict['fuel']['C1'].keys())
        loc_dict = self.core_dict['fuel']['C1']
        for i in range(len(genome)):
            iloc_tag = loc_tag[i]
            isym = 1 + len(loc_dict[iloc_tag]['Symmetric_Assemblies'])  
            igene=genome[i]
            if igene in fresh_fuel_list:
                feed_count+=isym
        return(feed_count)

    def generate_initial(self,chromosome_map):
        """
        Generates the initial solutions to the optimization problem.

        Parameters: 
            chromosome_map: Dictionary
                The genome portion of the dictionary settings file. 

        Written by Brian Andersen. 1/9/2020
        """
        chromosome_length = None
        chromosome_list = list(chromosome_map.keys())
        if 'symmetry_list' in chromosome_list:
            chromosome_list.remove('symmetry_list')

        for chromosome in chromosome_list:
            if chromosome_length is None:
                chromosome_length = len(chromosome_map[chromosome]['map'])
            elif len(chromosome_map[chromosome]['map']) == chromosome_length:
                pass
            else:
                raise ValueError("Chromosome Maps are of unequal length")

        fresh_fuel_list = []
        burnt_fuel_list = []
        self.genome = []
        # Separate fresh from burnt assemblies
        for i in range(len(chromosome_list)):
            if chromosome_map[chromosome_list[i]]['type'] =='reload':
                burnt_fuel_list.append(chromosome_list[i]) 
            else:
                fresh_fuel_list.append(chromosome_list[i]) 
        
        map_id=np.array(chromosome_map[chromosome_list[0]]['map'] )
        nfa = map_id[map_id > 0].shape[0]
        nfa_reload = int(49/len(burnt_fuel_list))
        fuel_types =  ['Fresh', 'Burnt']
        loc_tag = list(self.core_dict['fuel']['C1'].keys())
        loc_dict = self.core_dict['fuel']['C1']
        feed_count = np.zeros(self.ncycles, dtype=np.int32)
        fresh_feed_bool = False
        max_fresh_feed = np.inf
        if 'design_limits' in self.settings['genome']:
            if 'fresh_feed' in self.settings['genome']['design_limits']:
                max_fresh_feed = self.settings['genome']['design_limits']['fresh_feed']
                fresh_feed_bool = True
        
        for ci in range(self.ncycles):   
            burnt_fuel_cycle_list = copy.deepcopy(burnt_fuel_list)    
            avail_rfa=np.ones(len(burnt_fuel_list), dtype=np.int32)*nfa_reload       
            avail_rfa=avail_rfa.tolist()                              
            for i in range(chromosome_length):              
                no_gene_found = True 
                iloc_tag = loc_tag[i]
                isym = 1 + len(loc_dict[iloc_tag]['Symmetric_Assemblies'])                       
                while no_gene_found:
                    gene_type = random.choice(fuel_types)
                    if gene_type == 'Fresh':
                        gene = random.choice(fresh_fuel_list)
                    else:
                        gene = random.choice(burnt_fuel_cycle_list)
                    if chromosome_map[gene]['map'][i]:
                        if gene in burnt_fuel_cycle_list:
                            gene_id = burnt_fuel_cycle_list.index(gene)
                            avail_rfa[gene_id]-=1
                            if avail_rfa[gene_id]==0:
                                burnt_fuel_cycle_list.remove(gene)
                                avail_rfa.remove(0)
                            self.genome.append(gene)
                            no_gene_found = False
                        else:
                            ifeed = feed_count[ci] + isym
                            if fresh_feed_bool and ifeed>max_fresh_feed:
                                pass
                            else:
                                self.genome.append(gene)
                                feed_count[ci] = ifeed
                                no_gene_found = False

    def generate_initial_fixed(self,chromosome_map,gene_groups):
        """
        Generates initial solution when only specific number of assemblies
        may be used.

        Written by Brian Andersen 3/15/2020
        """
        chromosome_length = None
        chromosome_list = list(chromosome_map.keys())
        if 'symmetry_list' in chromosome_list:
            chromosome_list.remove('symmetry_list')

        for chromosome in chromosome_list:
            if chromosome_length is None:
                chromosome_length = len(chromosome_map[chromosome]['map'])
            elif len(chromosome_map[chromosome]['map']) == chromosome_length:
                pass
            else:
                raise ValueError("Chromosome Maps are of unequal length")

        no_valid_solution = True
        while no_valid_solution:
            no_valid_solution = False
            my_group = copy.deepcopy(gene_groups)
            self.genome = [None]*chromosome_length
            for i in range(chromosome_length):
                no_gene_found = True
                attempt_counter = 0
                while no_gene_found:
                    gene = random.choice(chromosome_list)
                    if 'unique' in chromosome_map[gene]:
                        if chromosome_map[gene]['unique']:
                            if gene in self.genome:
                                pass
                            else:
                                #This else loop activates if the gene is labeled unique but is not used. 
                                if chromosome_map[gene]['map'][i] == 1:
                                    if i in chromosome_map['symmetry_list']:
                                        if my_group[chromosome_map[gene]['gene_group']] > 1:
                                            self.genome[i] = gene
                                            no_gene_found = False
                                            my_group[chromosome_map[gene]['gene_group']] -= 2
                                    else:
                                        if my_group[chromosome_map[gene]['gene_group']] > 0:
                                            self.genome[i] = gene
                                            no_gene_found = False
                                            my_group[chromosome_map[gene]['gene_group']] -= 1            
                        else:
                            #adding unique loop above this code
                            if chromosome_map[gene]['map'][i] == 1:
                                if i in chromosome_map['symmetry_list']:
                                    if my_group[chromosome_map[gene]['gene_group']] > 1:
                                        self.genome[i] = gene
                                        no_gene_found = False
                                        my_group[chromosome_map[gene]['gene_group']] -= 2
                                else:
                                    if my_group[chromosome_map[gene]['gene_group']] > 0:
                                        self.genome[i] = gene
                                        no_gene_found = False
                                        my_group[chromosome_map[gene]['gene_group']] -= 1
                    else:
                        #adding unique loop above this code
                        if chromosome_map[gene]['map'][i] == 1:
                            if i in chromosome_map['symmetry_list']:
                                if my_group[chromosome_map[gene]['gene_group']] > 1:
                                    self.genome[i] = gene
                                    no_gene_found = False
                                    my_group[chromosome_map[gene]['gene_group']] -= 2
                            else:
                                if my_group[chromosome_map[gene]['gene_group']] > 0:
                                    self.genome[i] = gene
                                    no_gene_found = False
                                    my_group[chromosome_map[gene]['gene_group']] -= 1
                    attempt_counter += 1
                    if attempt_counter == 100:
                        no_gene_found = False
                        no_valid_solution = True

    def get_clength(self,efpd,boron,keff):
        if 0.1 in boron:
            eoc1_ind = 0
            eco2_ind = 0
            first_appear = False # bolean to identify boron=0.1 first appearance
            for i in range(len(efpd)-1):
                if boron[i] > 0.1 and boron[i+1] == 0.1 and first_appear == False:
                    eoc1_ind = i
                    eco2_ind = i+1
                    first_appear = True
            dbor = abs(boron[eoc1_ind-1]-boron[eoc1_ind])
            defpd = abs(efpd[eoc1_ind-1]-efpd[eoc1_ind])
            if dbor == 0 or eoc1_ind==0:
                def_dbor = 0.0
            else:
                def_dbor = defpd/dbor
            eoc = efpd[eoc1_ind] + def_dbor*(boron[eoc1_ind]-0.1)
            eoc = min(eoc, efpd[eoc1_ind+1])
        elif boron[-1]==boron[0]==1800.0:
            drho_dcb=10 
            drho1 = (keff[-2]-1.0)*10**5
            dcb1 = drho1/drho_dcb
            cb1= boron[-2] + dcb1
            drho2 = (keff[-1]-1.0)*10**5
            dcb2 = drho2/drho_dcb
            cb2= boron[-1] + dcb2
            dbor = abs(cb1-cb2)
            defpd = abs(efpd[-2]-efpd[-1])
            def_dbor = defpd/dbor
            eoc = efpd[-1] + def_dbor*(cb2-0.1)
        else:
            dbor = abs(boron[-2]-boron[-1])
            defpd = abs(efpd[-2]-efpd[-1])
            if dbor == 0:
                def_dbor = 0
            else:
                def_dbor = defpd/dbor
            eoc = efpd[-1] + def_dbor*(boron[-1]-0.1)
        if eoc == 0.0:
            eoc += 0.1
        return(eoc)

    def get_max_boron(self,boron,keff):
        res = max(boron)
        if res == 1800.0:
            max_boron =0
            for i in range(len(boron)):
                if boron[i]== 1800.0:
                    drho_dcb=10 
                    drho = (keff[i]-1.0)*10**5
                    dcb = drho/drho_dcb
                    mboron = 1800.0+dcb
                    if mboron > max_boron:
                        max_boron = mboron
            res = max_boron
        return(res)

    def get_lcoe(self):
        
        enr_map = {}
        enr_map['FE400']=4.00
        enr_map['FE401']=4.00
        enr_map['FE402']=4.00
        enr_map['FE450']=4.50
        enr_map['FE460']=4.60
        enr_map['FE461']=4.60
        enr_map['FE462']=4.60
        enr_map['FE501']=4.95
        enr_map['FE502']=4.95
        enr_map['FE526']=5.25
        enr_map['FE550']=5.50
        enr_map['FE566']=5.65
        enr_map['FE586']=5.85
        enr_map['FE600']=6.00
        enr_map['FE601']=6.00

        cycle1_param={'EFPD': self.parameters['cycle1_length']['target'],
                    'Batches': 3,
                    'Thermal_Power': self.power,
                    'Efficiency': 0.33,
                    'Fuel_Assemblies': self.number_assemblies}

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
                    'Storage_Time': 5.0+cycle1_param['EFPD']*cycle1_param['Batches']/365.25,
                    'Disposal_Time': cycle1_param['EFPD']*cycle1_param['Batches']/365.25}

        as_values=list(self.full_core['C1'].values())
        for i in range(len(as_values)):
            elt = as_values[i]
            if elt in self.core_dict['fuel']['C1']:
                as_values[i]='BRN'

        unique_fa =  np.unique(as_values)
        asb_param = {}
        for i in range(len(unique_fa)):
            nfa = as_values.count(unique_fa[i])
            if 'ASB' in unique_fa[i]:
                enr = 0.01
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 0.0
                            }
            else:
                asb=unique_fa[i]
                enr = enr_map[asb]/100
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 250
                            }
            asb_param[unique_fa[i]]=asb_dict

        lcoe_c1, bu, asb_cost = LCOE_MCYC(cycle1_param,lcoe_param, asb_param)

        cycle2_param={'EFPD': self.parameters['cycle2_length']['target'],
                    'Batches': 3,
                    'Thermal_Power': self.power,
                    'Efficiency': 0.33,
                    'Fuel_Assemblies': self.number_assemblies}

        lcoe2_param={'Discount_Rate': 0.07,
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
                    'Uranium_Ore_Time': -2.0 + cycle1_param['EFPD']/365.25,
                    'Conversion_Time': -1.5 + cycle1_param['EFPD']/365.25,
                    'Enrichment_Time': -1.0 + cycle1_param['EFPD']/365.25,
                    'Fabrication_Time': -0.5 + cycle1_param['EFPD']/365.25,
                    'Storage_Time': 5.0+ cycle1_param['EFPD']/365.25 +cycle2_param['EFPD']*cycle2_param['Batches']/365.25,
                    'Disposal_Time': cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']*cycle2_param['Batches']/365.25}

        as_values=list(self.full_core['C2'].values())
        for i in range(len(as_values)):
            elt = as_values[i]
            if elt in self.core_dict['fuel']['C2']:
                as_values[i]='BRN'

        unique_fa =  np.unique(as_values)
        asb2_param = {}
        for i in range(len(unique_fa)):
            nfa = as_values.count(unique_fa[i])
            if 'ASB' in unique_fa[i]:
                enr = 0.01
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 0.0
                            }
            else:
                asb=unique_fa[i]
                enr = enr_map[asb]/100
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 250
                            }
            asb2_param[unique_fa[i]]=asb_dict

        lcoe_c2, bu, asb_cost = LCOE_MCYC(cycle2_param,lcoe2_param, asb2_param)

        cycle3_param={'EFPD': self.parameters['cycle3_length']['target'],
                    'Batches': 3,
                    'Thermal_Power': self.power,
                    'Efficiency': 0.33,
                    'Fuel_Assemblies': self.number_assemblies}

        lcoe3_param={'Discount_Rate': 0.07,
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
                    'Uranium_Ore_Time': -2.0 + cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25,
                    'Conversion_Time': -1.5 + cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25,
                    'Enrichment_Time': -1.0 + cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25,
                    'Fabrication_Time': -0.5 + cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25,
                    'Storage_Time': 5.0 + cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25 +cycle3_param['EFPD']*cycle3_param['Batches']/365.25,
                    'Disposal_Time': cycle1_param['EFPD']/365.25 + cycle2_param['EFPD']/365.25 + cycle3_param['EFPD']*cycle3_param['Batches']/365.25}

        as_values=list(self.full_core['C3'].values())
        for i in range(len(as_values)):
            elt = as_values[i]
            if elt in self.core_dict['fuel']['C3']:
                as_values[i]='BRN'

        unique_fa =  np.unique(as_values)
        asb3_param = {}
        for i in range(len(unique_fa)):
            nfa = as_values.count(unique_fa[i])
            if 'ASB' in unique_fa[i]:
                enr = 0.01
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 0.0
                            }
            else:
                asb=unique_fa[i]
                enr = enr_map[asb]/100
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 250
                            }
            asb3_param[unique_fa[i]]=asb_dict

        lcoe_c3, bu, asb_cost = LCOE(cycle3_param,lcoe3_param, asb3_param)
        lcoe = (lcoe_c1 + lcoe_c2 + lcoe_c3)/3
        self.parameters["lcoe"]['value'] = lcoe
        return()

    def get_lcoe_cycle(self,ncycle):
        
        enr_map = {}
        enr_map['FE400']=4.00
        enr_map['FE401']=4.00
        enr_map['FE402']=4.00
        enr_map['FE450']=4.50
        enr_map['FE460']=4.60
        enr_map['FE461']=4.60
        enr_map['FE462']=4.60
        enr_map['FE501']=4.95
        enr_map['FE502']=4.95
        enr_map['FE526']=5.25
        enr_map['FE550']=5.50
        enr_map['FE566']=5.65
        enr_map['FE586']=5.85
        enr_map['FE600']=6.00
        enr_map['FE601']=6.00

        cycle_param={'EFPD': self.parameters['cycle_length']['target'],
                    'Batches': 3,
                    'Thermal_Power': self.power,
                    'Efficiency': 0.33,
                    'Fuel_Assemblies': self.number_assemblies}

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

        as_values=list(self.full_core['C'+str(ncycle)].values())
        as_values = [x[0] for x in as_values]
        for i in range(len(as_values)):
            elt = as_values[i]
            if elt in self.core_dict['fuel']['C'+str(ncycle)]:
                as_values[i]='BRN'

        unique_fa =  np.unique(as_values)
        asb_param = {}
        for i in range(len(unique_fa)):
            nfa = as_values.count(unique_fa[i])
            if 'ASB' in unique_fa[i]:
                enr = 0.01
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 0.0
                            }
            else:
                asb=unique_fa[i]
                enr = enr_map[asb]/100
                asb_dict = {'Number': nfa,
                            'Fuel_Rods': 264,
                            'Fuel_Radius': 0.41,
                            'Fuel_Height': 365.76,
                            'Enrichment': enr,
                            'Fuel_Density': 10.23,
                            'Fabrication_Price': 250
                            }
            asb_param[unique_fa[i]]=asb_dict

        lcoe, bu, asb_cost = LCOE_MCYC(cycle_param,lcoe_param, asb_param)
        self.parameters["lcoe"]['value'] = lcoe
        return()

    def get_cycle_results(self,filepath,ncyc=1):
        efpd=[]
        boron =[]
        fq=[]
        fdh=[]
        fxy=[]
        fxyz=[]
        keff = []
        read_bool  = False
        ofile = open(filepath + ".parcs_dpl", "r")
        filestr = ofile.read()
        ofile.close()
        res_str = filestr.split('===============================================================================')
        res_str = res_str[1].split('_______________________________________________________________________________')
        res_str = res_str[0].split('\n')
        for i in range(2, len(res_str)-1):
            res_val=res_str[i].split()
            fxy.append(float(res_val[6]))
            fxyz.append(float(res_val[7]))
            efpd.append(float(res_val[9]))
            boron.append(float(res_val[14]))
            keff.append(float(res_val[2]))
            fq.append(float(res_val[22]))
            fdh.append(float(res_val[21]))
        res = {}
        if ncyc==1:
            self.cycle_parameters = {}
        self.cycle_parameters['C'+str(ncyc)] = {}
        self.cycle_parameters['C'+str(ncyc)]["cycle_length"]= self.get_clength(efpd,boron,keff)       
        self.cycle_parameters['C'+str(ncyc)]["PinPowerPeaking"] = max(fq)
        self.cycle_parameters['C'+str(ncyc)]["FDeltaH"] = max(fdh)
        self.cycle_parameters['C'+str(ncyc)]["FXY"] = max(fxy)
        self.cycle_parameters['C'+str(ncyc)]["FXYZ"] = max(fxyz)
        self.cycle_parameters['C'+str(ncyc)]["max_boron"]= max(boron)
        if self.cycle_parameters['C'+str(ncyc)]["max_boron"] == 1800.0:
            max_boron =0
            for i in range(len(boron)):
                if boron[i]== 1800.0:
                    drho_dcb=10 
                    drho = (keff[i]-1.0)*10**5
                    dcb = drho/drho_dcb
                    mboron = 1800.0+dcb
                    if mboron > max_boron:
                        max_boron = mboron
            self.cycle_parameters['C'+str(ncyc)]["max_boron"] = max_boron

    def get_quarter_lattice(self):
            """
            Generates the 17x17 quarter core lattice.

            Parameters: None
        
            Written by Gregory Delipei 7/24/2022
            """
            core_map = self.core_dict['core_map']
            core_id = self.core_dict['core_id']
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

    def get_full_lattice(self):
            """
            Generates the 17x17 quarter core lattice.

            Parameters: None
        
            Written by Gregory Delipei 7/24/2022
            """
            core_map = self.core_dict['core_map']
            core_id = self.core_dict['core_id']
            nx = int(len(core_map[0]))
            ny = int(len(core_map))
            quarter_core = np.zeros((ny,nx), dtype='<U8')
            for y in range(ny):
                for x in range(nx):
                    val =core_map[y,x]
                    if val is None:
                        val = "00"
                    quarter_core[y,x] = val
            return(quarter_core)

    def compute_nasb(self):
        nfa=0
        for x in range(self.core_lattice.shape[0]):
            for y in range(self.core_lattice.shape[1]):
                loc = self.core_lattice[x,y]
                if loc != "00 " and loc != "10 ":
                    nfa+=1
                else:
                    pass
        return(nfa)

    def reproduce(self):
        if self.ncycles==1:
            new_genome=self.reproduce_cycle()
        else:
            new_genome=self.reproduce_mcycle()
        return(new_genome)

    def reproduce_cycle(self):
        old_genome = copy.deepcopy(self.genome)
        nfuel = len(self.core_dict['fuel']['C1'].keys())
        no_gene_found = True
        ncount = 0
        while no_gene_found:  
            pos_id = random.choice(np.arange(nfuel))
            cyc_genome = old_genome[0:56]
            inv = list(self.core_dict['Inventory'].keys())
            old_gene = old_genome[pos_id]
            fresh_feed_bool = False
            max_fresh_feed = np.inf
            if 'design_limits' in self.settings['genome']:
                if 'fresh_feed' in self.settings['genome']['design_limits']:
                    max_fresh_feed = self.settings['genome']['design_limits']['fresh_feed']
                    fresh_feed_bool = True
            burnt_fuel_list = []
            fresh_fuel_list = []
            for i in range(len(inv)):
                if 'FE' in inv[i]:
                    fresh_fuel_list.append(inv[i])
                else:
                    burnt_fuel_list.append(inv[i]) 

            nfa_reload = int(56/len(burnt_fuel_list))        
            inv.remove(old_gene)
       
            new_genome = copy.deepcopy(old_genome)
            new_gene = random.choice(inv)
            nfa=cyc_genome.count(new_gene)
            if new_gene in burnt_fuel_list and nfa == nfa_reload:
                cgen_id = random.choice(np.where(np.array(cyc_genome)==new_gene)[0])
                new_genome[pos_id] = new_gene
                new_genome[cgen_id] = old_gene
                ifeed = self.feed_counter(new_genome)
                if fresh_feed_bool and ifeed > max_fresh_feed:
                    pass
                else:
                    no_gene_found = False
            else:
                new_genome[pos_id] = new_gene
                ifeed = self.feed_counter(new_genome)
                if fresh_feed_bool and ifeed>max_fresh_feed:
                    pass
                else:
                    no_gene_found = False
            ncount +=1
            if ncount > 10000:
                new_genome = copy.deepcopy(old_genome)
                no_gene_found = False
        feed = self.feed_counter(new_genome)
        print('Number of tries: {}'.format(ncount))
        print('Fresh feed: {}'.format(feed))
        return(new_genome)

    def reproduce_mcycle(self):
        old_genome = self.genome
        nfuel = len(self.core_dict['fuel']['C1'].keys())
        no_gene_found = True
        ncount = 0
        while no_gene_found:  
            cyc_id=random.choice([1, 2, 3])
            pos_id = random.choice(np.arange(nfuel)) + (cyc_id-1)*nfuel
            cyc_genome = old_genome[(cyc_id-1)*56:(cyc_id)*56]
            inv = list(self.core_dict['Inventory'].keys())
            old_gene = old_genome[pos_id]
            fresh_feed_bool = False
            max_fresh_feed = np.inf
            if 'design_limits' in self.settings['genome']:
                if 'fresh_feed' in self.settings['genome']['design_limits']:
                    max_fresh_feed = self.settings['genome']['design_limits']['fresh_feed']
                    fresh_feed_bool = True
            burnt_fuel_list = []
            fresh_fuel_list = []
            for i in range(len(inv)):
                if 'FE' in inv[i]:
                    fresh_fuel_list.append(inv[i])
                else:
                    burnt_fuel_list.append(inv[i]) 

            nfa_reload = int(56/len(burnt_fuel_list))        
            inv.remove(old_gene)
       
            new_genome = copy.deepcopy(old_genome)
            new_gene = random.choice(inv)
            nfa=cyc_genome.count(new_gene)
            if new_gene in burnt_fuel_list and nfa == nfa_reload:
                cgen_id = random.choice(np.where(np.array(cyc_genome)==new_gene)[0]) + (cyc_id-1)*nfuel
                new_genome[pos_id] = new_gene
                new_genome[cgen_id] = old_gene
                no_gene_found = False
            else:
                new_genome[pos_id] = new_gene
                icyc_start = (cyc_id-1)*nfuel
                icyc_end = (cyc_id-1)*nfuel + nfuel
                ifeed = self.feed_counter(new_genome[icyc_start:icyc_end])
                if fresh_feed_bool and ifeed>max_fresh_feed:
                    pass
                else:
                    no_gene_found = False
            ncount +=1
            if ncount > 10000:
                new_genome = copy.deepcopy(old_genome)
                no_gene_found = False
        feed1 = self.feed_counter(new_genome[0:56])
        feed2 = self.feed_counter(new_genome[56:112])
        feed3 = self.feed_counter(new_genome[112:])
        print('Number of tries: {}'.format(ncount))
        print('Fresh feed: {}, {}, {}'.format(feed1,feed2, feed3))
        return(new_genome)

    def get_burnup(self,ofile,FULL_CORE=False):
        '''
        Read 2D and 3D burnup from PARCS .dep output files
        Some geometry predefined parameters are required.
        '''
        nfa=56
        if FULL_CORE:
            nfa=193
        bu_2d=[]
        full2quarter_indices = [97,  98,  99, 100, 101, 102, 103, 104,
                                112, 113, 114, 115, 116, 117, 118, 119,
                                127, 128, 129, 130, 131, 132, 133, 134,
                                142, 143, 144, 145, 146, 147, 148, 149,
                                156, 157, 158, 159, 160, 161, 162,
                                169, 170, 171, 172, 173, 174, 175,
                                181, 182, 183, 184, 185, 186,
                                190, 191, 192, 193]
        txt = Path(ofile).read_text()
        txt_dep = txt.split('   PT')[-1].split('EXP 2D MAP')[1].split('EXP 1D MAP')[0]
        txt_dep=txt_dep.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        counter=0
        for i in range(1,len(txt_dep)-1):
            line_dep = txt_dep[i].split()
            for j in range(1,len(line_dep)):
                val = float(line_dep[j])
                if val > 0.0:
                    bu_2d.append(val)
                else:
                    pass
        
        bu_2d = np.array(bu_2d)
        if FULL_CORE:
            new_bu2d = []
            for i in range(len(full2quarter_indices)):
                new_bu2d.append(bu_2d[full2quarter_indices[i]-1])
            bu_2d = new_bu2d

        z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
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
        bu_3d=np.zeros((nfa,nz))
        txt = Path(ofile).read_text()
        txt_dep = txt.split('   PT')[-1].split('EXP 3D MAP')[1].split(' I_D 2D MAP')[0]
        txt_dep=txt_dep.split('\n')
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
                        bu_3d[ifass-1, iz]=val
        
        if FULL_CORE:
            new_bu3d = np.zeros((56,nz))
            for i in range(len(full2quarter_indices)):
                new_bu3d[i,:] = bu_3d[full2quarter_indices[i]-1,:]
            bu_3d = new_bu3d
        # filter symmetric quarter core assemblies
        bu_2d_filtered = self.get_quarter_symmetry_values(bu_2d)
        bu_3d_filtered = self.get_quarter_symmetry_values(bu_3d,axis=0)
        return((bu_2d_filtered, bu_3d_filtered))
    
    def get_alldep(self,ofile):
        '''
        Extract all results from a PARCS .parcs_dep file
        '''
        nfa=len(self.core_dict['fuel']['C1'].keys())
        txt = Path(ofile).read_text()
        txt_dep = txt.split('   PT')[1:]
        nsteps = len(txt_dep)
        res_dir = {}
        for i in range(nsteps):
            ires = {}
            itxt_dep = txt_dep[i]
            # Scalar Values
            itxt_scalar = itxt_dep.split('\n')[1]
            itxt_scalar = itxt_scalar.split()
            ires['KEFF'] = float(itxt_scalar[2])
            ires['AXOFF'] = float(itxt_scalar[4])
            ires['PZ'] = float(itxt_scalar[5])
            ires['PXY'] = float(itxt_scalar[6])
            ires['PXYZ'] = float(itxt_scalar[7])
            ires['PPIN'] = float(itxt_scalar[8])
            ires['EFPD'] = float(itxt_scalar[9])
            ires['BUAVE'] = float(itxt_scalar[10])
            ires['BUMAX'] = float(itxt_scalar[11])
            ires['BETA'] = float(itxt_scalar[12])
            ires['BORON'] = float(itxt_scalar[14])
            ires['TFUEL'] = float(itxt_scalar[15])
            ires['TMOD'] = float(itxt_scalar[16])
            ires['DMOD'] = float(itxt_scalar[17])
            # 2D power
            pow_2d=np.zeros(nfa)
            itxt_pow2d = itxt_dep.split(' RPF 2D MAP')[1].split(' RPF 1D MAP')[0]
            itxt_pow2d = itxt_pow2d.split('\n')
            itxt_pow2d=list(filter(lambda a: a != '', itxt_pow2d))
            itxt_pow2d=list(filter(lambda a: a != ' ', itxt_pow2d))
            counter=0
            for j in range(1,len(itxt_pow2d)-1):
                line_dep = itxt_pow2d[j].split()
                for k in range(1,len(line_dep)):
                    val = float(line_dep[k])
                    if val > 0.0:
                        pow_2d[counter]=val
                        counter+=1
                    else:
                        pass
            ires['POW2D']=np.array(pow_2d)
            # 2D burnup
            bu_2d=np.zeros(nfa)
            itxt_bu2d = itxt_dep.split(' EXP 2D MAP')[1].split(' EXP 1D MAP')[0]
            itxt_bu2d = itxt_bu2d.split('\n')
            itxt_bu2d=list(filter(lambda a: a != '', itxt_bu2d))
            itxt_bu2d=list(filter(lambda a: a != ' ', itxt_bu2d))
            counter=0
            for j in range(1,len(itxt_bu2d)-1):
                line_dep = itxt_bu2d[j].split()
                for k in range(1,len(line_dep)):
                    val = float(line_dep[k])
                    if val > 0.0:
                        bu_2d[counter]=val
                        counter+=1
                    else:
                        pass
            ires['BU2D']=np.array(bu_2d)

            # 3D 
            z_id=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
            nz=len(z_id)
            refl_id=[9,18,27,36,44,45,53,60,61,66,67,68,69,70,71,72,73]

            # 3D Power
            pow_3d=np.zeros((nfa,nz))
            itxt_pow3d = itxt_dep.split(' RPF 3D MAP 1.0E+00')[1].split(' EXP 2D MAP')[0]
            itxt_pow3d=itxt_pow3d.split('\n')
            itxt_pow3d=list(filter(lambda a: a != '', itxt_pow3d))
            itxt_pow3d=list(filter(lambda a: a != ' ', itxt_pow3d))
            asb_counter=0
            fasb_counter=0
            ifass = 0
            iass = 0
            for j in range(1,len(itxt_pow3d)):
                line_dep = itxt_pow3d[j].split()
                if line_dep[0]=='k':
                    asb_counter=copy.deepcopy(iass)
                    fasb_counter=copy.deepcopy(ifass)
                elif int(line_dep[0]) not in z_id: 
                    pass
                else:
                    iz = z_id[-1] - int(line_dep[0])
                    ifass = copy.deepcopy(fasb_counter)
                    iass = copy.deepcopy(asb_counter)
                    for k in range(1,len(line_dep)):
                        iass +=1
                        if iass in refl_id:
                            pass
                        else:
                            ifass +=1
                            val = float(line_dep[k])
                            pow_3d[ifass-1, iz]=val
            ires['POW3D'] =  np.flip(pow_3d, axis=1)

            # 3D Burnup
            bu_3d=np.zeros((nfa,nz))
            itxt_bu_3d = itxt_dep.split(' EXP 3D MAP 1.0E+00')[1].split(' END STEP')[0]
            itxt_bu_3d=itxt_bu_3d.split('\n')
            itxt_bu_3d=list(filter(lambda a: a != '', itxt_bu_3d))
            itxt_bu_3d=list(filter(lambda a: a != ' ', itxt_bu_3d))
            asb_counter=0
            fasb_counter=0
            ifass = 0
            iass = 0
            for j in range(1,len(itxt_bu_3d)):
                line_dep = itxt_bu_3d[j].split()
                if line_dep[0]=='k':
                    asb_counter=copy.deepcopy(iass)
                    fasb_counter=copy.deepcopy(ifass)
                elif int(line_dep[0]) not in z_id: 
                    pass
                else:
                    iz = z_id[-1] - int(line_dep[0])
                    ifass = copy.deepcopy(fasb_counter)
                    iass = copy.deepcopy(asb_counter)
                    for k in range(1,len(line_dep)):
                        iass +=1
                        if iass in refl_id:
                            pass
                        else:
                            ifass +=1
                            val = float(line_dep[k])
                            bu_3d[ifass-1, iz]=val
            ires['BU3D'] =  np.flip(bu_3d, axis=1)
            res_dir['STEP_'+str(i+1)] = ires
        return(res_dir)

    def get_reac(self,ass,bu):
        reac = []
        xsDict =  pickle.load(open( self.settings['genome']['parcs_data']['xs_library'] + '/xs_dts_mcyc.p', "rb" ) )
        for i in range(len(ass)):
            for key, value in self.settings['genome']['chromosomes'].items():
                if str(value['type']) == str(ass[i]):
                    xs_tag=value['serial'][4:]
            bulist= list(xsDict[xs_tag].keys())
            bu0=0
            bu1=0
            reac0=0
            reac1=0
            for j in range(len(bulist)-1):
                ibu0 = float(bulist[j][3:])
                ibu1 = float(bulist[j+1][3:])
                if bu[i]>ibu0 and bu[i]<ibu1:
                    bu0=ibu0
                    bu1=ibu1
                    reac0 = xsDict[xs_tag][bulist[j]]['keff']
                    reac1 = xsDict[xs_tag][bulist[j+1]]['keff']
            if bu0==0 and bu1==0:
                bu0=float(bulist[-2][3:])
                bu1=float(bulist[-1][3:])
                reac0 = xsDict[xs_tag][bulist[-2]]['keff']
                reac1 = xsDict[xs_tag][bulist[-1]]['keff']
                dreac_dbu = (reac1-reac0)/(bu1-bu0)
                ireac = reac1 + dreac_dbu*(bu[i]-bu1)
            else:
                dreac_dbu = (reac1-reac0)/(bu1-bu0)
                ireac = reac0 + dreac_dbu*(bu[i]-bu0)
            reac.append(ireac)
        return(reac)

    def rank_assb(self):
        reac = []
        asb = []
        for key, value in self.reload_inventory.items():
            asb.append(key)
            reac.append(value['REAC'])
    
        reac = np.array(reac)
        reac_sorted = np.argsort(reac)
        chromosomes=self.settings['genome']['chromosomes']
        burnt_fuel_list=[]
        for key,value in chromosomes.items():
            if value['type'] =='reload':
                burnt_fuel_list.append(key)

        nasb = len(self.reload_inventory)      
        nfa_reload = int(nasb/len(burnt_fuel_list))
        ranked_assb = {}
        asb_group_count=0
        assb_fuel_tags=[]
        assb_fuel_reac=[]
        for i in range(nasb):
            if i>0 and i%nfa_reload==0 and i<nfa_reload*len(burnt_fuel_list):
                assb_dict = {}
                assb_dict['fuel_tag']=assb_fuel_tags
                assb_dict['fuel_reactivity']=assb_fuel_reac
                ranked_assb[burnt_fuel_list[asb_group_count]]=assb_dict
                assb_fuel_tags=[]
                assb_fuel_reac=[]
                asb_group_count+=1
            reac_id = reac_sorted[i]
            assb_fuel_reac.append(reac[reac_id])
            assb_fuel_tags.append(asb[reac_id])
            self.reload_inventory[asb[reac_id]]['GROUP']=burnt_fuel_list[asb_group_count]
        assb_dict = {}
        assb_dict['fuel_tag']=assb_fuel_tags
        assb_dict['fuel_reactivity']=assb_fuel_reac
        ranked_assb[burnt_fuel_list[asb_group_count]]=assb_dict
        return(ranked_assb)

    @trials
    def getlp_from_genome(self,cycle_genome,cycle_assb,fuel_locations,cycle_tag):
        avail_assb = copy.deepcopy(cycle_assb)
        reload_inventory = copy.deepcopy(self.reload_inventory)
        cycle_lp=[]
        
        for i in range(len(cycle_genome)):
            gene = cycle_genome[i]
            loc = fuel_locations[i]
            loc_symmetry = len(self.core_dict['fuel'][cycle_tag][loc]['Symmetric_Assemblies']) + 1 
            if self.settings['genome']['chromosomes'][gene]['type']=='reload':
                avail_reload_genes = []
                reload_genes=avail_assb[gene]['fuel_tag']
                for rgene in reload_genes:
                    if reload_inventory[rgene]['QTY'] >= loc_symmetry:
                        avail_reload_genes.append(rgene)
                selected_rgene = random.choice(avail_reload_genes)
                cycle_lp.append(selected_rgene)
                reload_inventory[selected_rgene]['QTY'] -= loc_symmetry
                if reload_inventory[selected_rgene]['QTY'] == 0:
                    selected_rgene_id=avail_assb[gene]['fuel_tag'].index(selected_rgene)
                    selected_rgene_reac = avail_assb[gene]['fuel_reactivity'][selected_rgene_id]
                    avail_assb[gene]['fuel_tag'].remove(selected_rgene)
                    avail_assb[gene]['fuel_reactivity'].remove(selected_rgene_reac)
            else:
                cycle_lp.append(gene)
        self.reload_inventory = copy.deepcopy(reload_inventory)
        return(cycle_lp)

    def test_fail(self,fpath):
        val=False
        if os.path.exists(fpath):
            val = True
        return(val)

    def write_exp(self,exp_file,ncyc):
        full_core_lp = []
        for iy in range(self.core_lattice.shape[0]):
            for ix in range(self.core_lattice.shape[1]):
                loc = self.core_lattice[iy,ix]
                if loc != "00" and loc[0] != "R":
                    full_core_lp.append(loc)
                else:
                    pass
        with open(exp_file,"w") as ofile:             
            ofile.write("\n")
            ofile.write(" BEGIN STEP\n")
            ofile.write("\n")
            ofile.write(" EXP 3D MAP 1.0E+00\n")
            ofile.write("\n")
            asb_counter = 0
            ncol = 10
            nz=28
            nblocks = int(len(full_core_lp)/ncol)
            for nb in range(nblocks):
                ofile.write(" k lb")
                for nc in range(1,ncol+1):
                   ofile.write('{:7d}'.format(asb_counter+nc)) 
                   ofile.write(" ")
                ofile.write("\n")
                for iz in range(nz):
                    zid = nz-iz
                    ofile.write('{:3d}'.format(zid))
                    ofile.write("  ")
                    for nc in range(ncol):
                        iasb = asb_counter + nc
                        iloc = full_core_lp[iasb]
                        fasb = self.full_core['C'+str(ncyc)][iloc][0]
                        if fasb in self.reload_inventory:
                            val = self.reload_inventory[fasb]['BU3D'][iz]
                        else:
                            val = 0.0
                        ofile.write(r'{:7.3f}'.format(val))
                        ofile.write(' ')
                    ofile.write("\n")
                asb_counter+=nc+1
                ofile.write("\n")
            if asb_counter<len(full_core_lp):
                ofile.write(" k lb")
                for nc in range(asb_counter+1,len(full_core_lp)+1):
                   ofile.write('{:7d}'.format(nc)) 
                   ofile.write(" ")
                ofile.write("\n")
                for iz in range(nz):
                    zid = nz-iz
                    ofile.write('{:3d}'.format(zid))
                    ofile.write("  ")
                    for nc in range(asb_counter+1,len(full_core_lp)+1):
                        iasb = nc-1
                        iloc = full_core_lp[iasb]
                        fasb = self.full_core['C'+str(ncyc)][iloc][0]
                        if fasb in self.reload_inventory:
                            val = self.reload_inventory[fasb]['BU3D'][iz]
                        else:
                            val = 0.0
                        ofile.write(r'{:7.3f}'.format(val))
                        ofile.write(' ')
                    ofile.write("\n")
                asb_counter=nc
            ofile.write("\n")
            ofile.write(" END STEP")
            ofile.write("\n")
        return

    def clean_inventory(self):
        keys = list(self.reload_inventory.keys())
        for key in keys:
            idict = self.reload_inventory[key]
            if idict['QTY']==0:
                del self.reload_inventory[key]
        return

    def update_discharge_inventory(self,max_bu):
        keys = list(self.reload_inventory.keys())
        for key in keys:
            idict = self.reload_inventory[key]
            if idict['BU2D']>=max_bu:
                self.discharge_inventory[key]=idict
                del self.reload_inventory[key]
        return

    def store_cycle(self,opt,ftag,cycle_lp,bu2d_eoc,ncyc):
        '''
        Store cycle depletion results in a binary dictionary (ldts.p or hdts.p)
        Two options: 
            - light storage stores only cycle objectives and loading pattern
            - heavy storage stores more quantities
        '''
        dts_fpath = "ldts.p" if opt=='light' else "hdts.p"
        rdict = {}
        for key, value in self.cycle_parameters['C'+str(ncyc)].items():
            rdict[key]= value
        
        xs_type = []
        fuel_type= []
        burnup_2d = []
        lp = []
        for i in range(len(cycle_lp)):
            iasb = cycle_lp[i]
            if 'FE' in iasb:
                xs_val = self.core_dict['Inventory'][iasb]['Cross_Section']
                fuel_val = int(self.core_dict['Inventory'][iasb]['Tag'])
                bu_val = 0.0
                lp_val = iasb
            else:
                fresh_iasb = 'FE' + str(self.reload_inventory[iasb]['TYPE'])
                xs_val = self.core_dict['Inventory'][fresh_iasb]['Cross_Section']
                fuel_val = int(self.reload_inventory[iasb]['TYPE'])
                bu_val = self.reload_inventory[iasb]['BU2D']
                lp_val = self.reload_inventory[iasb]['GROUP']
            xs_type.append(xs_val)
            fuel_type.append(fuel_val)
            burnup_2d.append(bu_val)
            lp.append(lp_val)
        rdict['LP_XS']=np.array(xs_type)
        rdict['LP_TYPE']=np.array(fuel_type)
        rdict['LP_BU2D_BOC']=np.array(burnup_2d)
        rdict['LP_BU2D_EOC']=np.array(bu2d_eoc)
        rdict['LP_ASB']=np.array(cycle_lp)
        rdict['LP_FUEL']=np.array(lp)
        if opt == 'heavy':
            nz = 16
            burnup_3d = np.zeros((len(burnup_2d),nz))
            for i in range(len(cycle_lp)):
                iasb = cycle_lp[i]
                if 'FE' in iasb:
                    bu_val = np.zeros(nz)
                else:
                    bu_val = self.reload_inventory[iasb]['BU3D']
                burnup_3d[i,:] = bu_val
            rdict['LP_BU3D']=np.flip(burnup_3d, axis=1)
            ofile = ftag + '.parcs_dep'
            res_dir = self.get_alldep(ofile)
            rdict['RES'] = res_dir
        pickle.dump( rdict, open( dts_fpath, "wb" ) )
        return
    
    def set_poor_results(self, ncyc=None):
        '''
        A function that assigns bad results to the parameters in purpose
        '''
        if ncyc is None:
            self.parameters['max_boron']['value']=np.random.uniform(0,10)
            self.parameters['PinPowerPeaking']['value']=np.random.uniform(10,20)
            self.parameters['FDeltaH']['value']=np.random.uniform(10,20)
            self.parameters['cycle_length']['value'] =  np.random.uniform(5000,10000) 
            self.parameters["lcoe"]['value'] = np.random.uniform(100,200) 
        else:
            self.cycle_parameters['C'+str(ncyc)] = {}
            self.cycle_parameters['C'+str(ncyc)]["cycle_length"] = np.random.uniform(0,10)
            self.cycle_parameters['C'+str(ncyc)]["PinPowerPeaking"] = np.random.uniform(10,20)
            self.cycle_parameters['C'+str(ncyc)]["FDeltaH"]= np.random.uniform(10,20)
            self.cycle_parameters['C'+str(ncyc)]["FXY"] = np.random.uniform(10,20)
            self.cycle_parameters['C'+str(ncyc)]["FXYZ"]= np.random.uniform(10,20)
            self.cycle_parameters['C'+str(ncyc)]["max_boron"] = np.random.uniform(5000,10000) 

    def run_cycle(self,fuel_loc,c0_assb,cycle_lp,rdir,ncyc=1):
        pwd = Path(os.getcwd())
        run_directory = pwd / rdir
        if not os.path.exists(run_directory):
            os.makedirs(run_directory)
        else:
            shutil.rmtree(run_directory, ignore_errors=True)
            os.makedirs(run_directory)
        
        self.update_core_lp(cycle_lp,c0_assb,fuel_loc,ncyc)

        cdir = self.library
        exp_file = 'cyc_exp.dep' 
        dist_file = run_directory / exp_file
        self.write_exp(dist_file,ncyc)
        
        os.chdir(run_directory)
        print('{} cycle {} calculation starts at {}!'.format(self.name,ncyc,os.getcwd()))

        core_cycle_lattice = copy.deepcopy(self.core_lattice)
        xs_array_cycle = np.zeros((core_cycle_lattice.shape[0],core_cycle_lattice.shape[1]), dtype='<U20')
        pincal_loc = np.zeros((core_cycle_lattice.shape[0],core_cycle_lattice.shape[1]))
        for x in range(core_cycle_lattice.shape[0]):
            for y in range(core_cycle_lattice.shape[1]):
                loc = core_cycle_lattice[x,y]
                if loc != "00" and loc[0] != "R":
                    sel_asb=self.full_core['C'+str(ncyc)][loc][0]
                    if sel_asb in self.reload_inventory:
                        asb_type = self.reload_inventory[sel_asb]['TYPE']
                        fresh_ass = 'FE' + str(asb_type)
                        xs_val = self.core_dict['Inventory'][fresh_ass]['Cross_Section']
                    else:
                        fresh_ass = sel_asb
                        xs_val = self.core_dict['Inventory'][fresh_ass]['Cross_Section']
                    xs_array_cycle[x,y] = xs_val
                    if loc in self.core_dict['fuel']['C'+str(ncyc)] or loc in ['H09','H10','H11','H12','H13','H14','H15']:
                        pincal_loc[x,y]=1
                    else:
                        pincal_loc[x,y]=0
                elif loc[0] == "R":
                    xs_array_cycle[x,y] = None
                    pincal_loc[x,y]=0
                elif loc == "00":
                    xs_array_cycle[x,y] = None
                    pincal_loc[x,y]=np.nan
        
        xs_unique_cycle = np.unique(xs_array_cycle)
        xs_unique_cycle = np.delete(xs_unique_cycle, np.argwhere(xs_unique_cycle == 'None'))
        xs_forced_cycle = np.array([self.core_dict['Inventory']['FE461']['Cross_Section'],
                        self.core_dict['Inventory']['FE462']['Cross_Section'],
                        self.core_dict['Inventory']['FE501']['Cross_Section'],
                        self.core_dict['Inventory']['FE502']['Cross_Section']])
        xs_unique_cycle = np.append(xs_unique_cycle,xs_forced_cycle)
        xs_unique_cycle = np.unique(xs_unique_cycle)

        tag_unique_cycle = copy.deepcopy(xs_unique_cycle)
        xs_ref = np.arange(5,5+len(xs_unique_cycle)) # 1-3 for reflectors and 4 for blankets
        for key,value in self.core_dict["Inventory"].items():
            for i in range(xs_unique_cycle.shape[0]):
                if value['Cross_Section']==xs_unique_cycle[i]:
                    tag_unique_cycle[i]=value['Tag']
        fname = 'solution'
        filename = fname + '.inp'
        with open(filename,"w") as ofile:             
            ofile.write("!******************************************************************************\n")
            ofile.write('CASEID {}  \n'.format(fname))
            ofile.write("!******************************************************************************\n\n")

        with open(filename,"a") as ofile:             
            ofile.write("CNTL\n")
            ofile.write("     RUN_OPTS F T F F\n")
            ofile.write("     INT_TH     T -1\n")
            ofile.write("     CORE_POWER 100.0\n")
            ofile.write("     CORE_TYPE  PWR\n")
            ofile.write("     PPM        1000\n")
            ofile.write("     DEPLETION  T  1.0E-5 T\n")
            ofile.write("     TREE_XS    T  0 T  T  F  F  T  F  T  F  T  F  T  T  T  F \n")
            ofile.write("     BANK_POS   100 100 100 100 100 100\n")
            ofile.write("     XE_SM      1 1 1 1\n")
            ofile.write("     SEARCH     PPM 1.0 1800.0 10.0\n")
            ofile.write("     MULT_CYC   F\n")
            ofile.write("     PIN_POWER  T\n")
            ofile.write("     PLOT_OPTS 0 0 0 0 0 2\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
            
        with open(filename,"a") as ofile:             
            ofile.write("PARAM\n")
            ofile.write("     LSOLVER  1 1 20\n")
            ofile.write("     NODAL_KERN     HYBRID\n")
            ofile.write("     CMFD     2\n")
            ofile.write("     DECUSP   2\n")
            ofile.write("     INIT_GUESS 0\n")
            ofile.write("     CONV_SS   1.e-6 1.e-5 1.e-5\n")
            ofile.write("     EPS_ANM   0.000001\n")
            ofile.write("     NLUPD_SS  3 5 1\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
        

        with open(filename,"a") as ofile:             
            ofile.write("GEOM\n")
            ofile.write("     GEO_DIM 17 17 30 1 1\n")
            ofile.write("     RAD_CONF\n")
            for iy in range(self.core_lattice.shape[0]):
                ofile.write('     ')
                for ix in range(self.core_lattice.shape[1]):
                    iloc = self.core_lattice[iy,ix]
                    if iloc == '00':
                        value = '00 '
                    elif iloc[0] == 'R':
                        value = '10 '
                    else:
                        value = str(self.full_core['C'+str(ncyc)][iloc][1])
                    ofile.write(value + '  ')
                ofile.write('\n')
            ofile.write('\n')
            ofile.write("     GRID_X      17*21.50\n")
            ofile.write("     NEUTMESH_X  17*1\n")
            ofile.write("     GRID_Y      17*21.50\n")
            ofile.write("     NEUTMESH_Y  17*1\n")
            ofile.write("     GRID_Z      30.48 15.24 10.16 5.08 22*13.85455 5.08 10.16 15.24 30.48\n")            
            ofile.write("     ASSY_TYPE   10   1*2   28*2    1*2 REFL\n")
            for i in range(xs_unique_cycle.shape[0]):
                if 'gd_0' in xs_unique_cycle[i]:
                    ofile.write("     ASSY_TYPE   {}   1*1 1*4  26*{}  1*4  1*3 FUEL\n".format(tag_unique_cycle[i],xs_ref[i]))
                else:
                    ofile.write("     ASSY_TYPE   {}   1*1 1*4  1*4 24*{} 1*4 1*4  1*3 FUEL\n".format(tag_unique_cycle[i],xs_ref[i]))
            ofile.write("\n")

            ofile.write("     boun_cond   2 2 2 2 2 2\n")
            ofile.write("     PINCAL_LOC\n")
            for x in range(pincal_loc.shape[0]):
                ofile.write("      ")
                for y in range(pincal_loc.shape[1]):
                    val = pincal_loc[x,y]
                    if np.isnan(val):
                        val = ' '
                        ofile.write(val)
                        ofile.write("  ")
                    else:
                        ofile.write(str(int(pincal_loc[x,y])))
                        ofile.write("  ")
                ofile.write("\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
        
        with open(filename,"a") as ofile:             
            ofile.write("FDBK\n")
            ofile.write("     FA_POWPIT     {} 21.5\n".format(np.round(self.power/193,4)))
            ofile.write("     GAMMA_FRAC    0.0208    0.0    0.0\n")
            ofile.write("     CDC_DED   0 1\n")
            ofile.write("     EFF_DOPLT   T  0.5556\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")


        with open(filename,"a") as ofile:   
            ofile.write("TH\n")          
            ofile.write("     FLU_TYP       0\n")
            ofile.write("     N_PINGT    264 25\n")
            ofile.write("     PIN_DIM      4.1 4.75 0.58 6.13\n")
            ofile.write("     FLOW_COND    {}  {}\n".format(np.round(self.inlet_temperature-273.15,2),np.round(self.flow/193,4)))
            ofile.write("     HGAP     10000.0\n")
            ofile.write("     N_RING   6\n")
            ofile.write("     THMESH_X       17*1\n")
            ofile.write("     THMESH_Y       17*1\n")
            ofile.write("     THMESH_Z       1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30\n")
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")

        with open(filename,"a") as ofile:             
            ofile.write("DEPL\n")
            if ncyc==1:
                ofile.write("     TIME_STP  1 1 17*30 13\n")
            else:
                ofile.write("     TIME_STP  1 1 22*30 14 14\n")
            ofile.write("     INP_HST   './cyc_exp.dep' -2 1\n")
            ofile.write("     PMAXS_F   1 '{}' 1\n".format(cdir + '/' + 'xs_gbot'))
            ofile.write("     PMAXS_F   2 '{}' 2\n".format(cdir + '/' + 'xs_grad'))
            ofile.write("     PMAXS_F   3 '{}' 3\n".format(cdir + '/' + 'xs_gtop'))
            ofile.write("     PMAXS_F   4 '{}' 4\n".format(cdir + '/' + 'xs_g250_gd_0_wt_0'))
            for i in range(xs_unique_cycle.shape[0]):
                ofile.write("     PMAXS_F   {} '{}' {}\n".format(5+i,cdir + '/' + xs_unique_cycle[i],5+i))
            ofile.write("\n")
            ofile.write("!******************************************************************************\n\n")
            ofile.write(".")

        # Run PARCS INPUT DECK
        
        parcscmd = "/cm/shared/nuclearCodes/parcs-3.4.3/PARCS-v343_Exe/Executables/Linux/parcs-v343-linux2-intel-x64-release.x"
        print('Execute PARCS')
        print('Running in process')
        try:
            #
            output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=180)
            # Get Results
            if 'Finished' in str(output):
                ofile = fname + '.out'
                self.get_cycle_results(fname,ncyc)
                

                ## update reload inventory
                bu_2d, bu_3d=self.get_burnup(fname+'.parcs_dep',FULL_CORE=True)
                reac = self.get_reac(c0_assb,bu_2d)

                # Store Optionally 
                if 'options' in self.settings:
                    if 'store' in self.settings['options']:
                        store_op = self.settings['options']['store']
                        self.store_cycle(store_op,fname,cycle_lp,bu_2d,ncyc)

                for i in range(len(cycle_lp)):
                    iasb = cycle_lp[i]
                    if iasb in self.reload_inventory:
                        self.reload_inventory[iasb]['BU2D']=bu_2d[i]
                        self.reload_inventory[iasb]['BU3D']=bu_3d[i,:]
                        self.reload_inventory[iasb]['REAC']=reac[i]
                    else:
                        idict = {}
                        idict['BU2D']=bu_2d[i]
                        idict['BU3D']=bu_3d[i,:]
                        idict['REAC']=reac[i]
                        idict['TYPE']=c0_assb[i]
                        idict['QTY']= len(self.core_dict['fuel']['C'+str(ncyc)][fuel_loc[i]]['Symmetric_Assemblies']) + 1
                        idict['LOC'+str(ncyc)]=fuel_loc[i]
                        self.reload_inventory_counter +=1 
                        self.reload_inventory['ASB'+str(self.reload_inventory_counter)] = idict
                
                # Clean

                os.system('rm -f {}.parcs_pin*'.format(fname))
                # os.system('rm -f {}.inp'.format(fname))
                os.system('rm -f {}.inp_parcs_err'.format(fname))
                os.system('rm -f {}.inp_paths_err'.format(fname))
                os.system('rm -f {}.parcs_dep'.format(fname))
                os.system('rm -f {}.parcs_itr'.format(fname))
                os.system('rm -f {}.parcs_msg'.format(fname))
                os.system('rm -f {}.parcs_out'.format(fname))
                os.system('rm -f {}.parcs_sum'.format(fname))
                os.system('rm -f {}.parcs_xml'.format(fname))
                os.system('rm -f {}.parcs_rst'.format(fname))
                #os.system('rm -f cyc_exp.dep')
            else:
                if ncyc==1:
                    self.cycle_parameters = {}
                self.set_poor_results(ncyc=ncyc)
                os.system('rm -f ./*')

        except subprocess.TimeoutExpired:
            print('Timed out - killing')
            
            os.system('rm -f {}.parcs_pin*'.format(fname))
            if ncyc==1:
                self.cycle_parameters = {}
            self.set_poor_results(ncyc=ncyc)
            os.system('rm -f ./*')

        print('{} cycle {} calculation is done at {}!'.format(self.name,ncyc,os.getcwd()))
        os.chdir(pwd)
        gc.collect()  
        print('finished collecting garbage...')
        print('exiting cycle evaluate...')
        return()

    def run_ml_cycle(self,fuel_loc,c0_assb,cycle_lp,rdir,ncyc=1):
        pwd = Path(os.getcwd())
        run_directory = pwd / rdir
        if not os.path.exists(run_directory):
            os.makedirs(run_directory)
        else:
            shutil.rmtree(run_directory, ignore_errors=True)
            os.makedirs(run_directory)
        
        os.chdir(run_directory)
        print('{} cycle {} calculation starts at {} with ML model!'.format(self.name,ncyc,os.getcwd()))
        

        # Load ML model
        ml_model_path = self.settings['genome']['ml_model']['model']
        # Load the XGBoost model
        with open(ml_model_path, 'rb') as model_file:
            ml_model = pickle.load(model_file)
        import pdb; pdb.set_trace()

        
        
        # Run PARCS INPUT DECK
        
        parcscmd = "/cm/shared/nuclearCodes/parcs-3.4.3/PARCS-v343_Exe/Executables/Linux/parcs-v343-linux2-intel-x64-release.x"

        print('Execute PARCS')
        print('Running in process')
        #
        output = subprocess.check_output([parcscmd, filename], stderr=STDOUT, timeout=120)
        # Get Results
       
        ofile = fname + '.out'
        self.get_cycle_results(fname,ncyc)

        # Store Optionally 
        if 'options' in self.settings:
            if 'store' in self.settings['options']:
                store_op = self.settings['options']['store']
                self.store_cycle(store_op,fname,cycle_lp,ncyc)

        ## update reload inventory
        bu_2d, bu_3d=self.get_burnup(fname+'.parcs_dep')
        reac = self.get_reac(c0_assb,bu_2d)
        for i in range(len(cycle_lp)):
            iasb = cycle_lp[i]
            if iasb in self.reload_inventory:
                self.reload_inventory[iasb]['BU2D']=bu_2d[i]
                self.reload_inventory[iasb]['BU3D']=bu_3d[i,:]
                self.reload_inventory[iasb]['REAC']=reac[i]
            else:
                idict = {}
                idict['BU2D']=bu_2d[i]
                idict['BU3D']=bu_3d[i,:]
                idict['REAC']=reac[i]
                idict['TYPE']=c0_assb[i]
                idict['QTY']= len(self.core_dict['fuel']['C'+str(ncyc)][fuel_loc[i]]['Symmetric_Assemblies']) + 1
                idict['LOC'+str(ncyc)]=fuel_loc[i]
                self.reload_inventory_counter +=1 
                self.reload_inventory['ASB'+str(self.reload_inventory_counter)] = idict

        print('{} cycle {} calculation is done at {}!'.format(self.name,ncyc,os.getcwd()))
        os.chdir(pwd)
        gc.collect()  
        print('finished collecting garbage...')
        print('exiting cycle evaluate...')
        return()

    def get_quarter_symmetry_values(self,values,axis=None):
        exclude_indices = [1,2,3,4,5,6,7]
        values = np.array(values)
        values_quarter = np.delete(values, exclude_indices,axis=axis)
        return(values_quarter)

    def update_core_lp(self,asb_values,asb_types,fuel_locations,ncycle):
        '''
        Function that updates the whole core loading pattern information based on a list of values in symmetry
        '''
        
        cyc_tag = 'C' + str(ncycle) 
            
        for i in range(len(fuel_locations)):
            self.core_dict['fuel'][cyc_tag][fuel_locations[i]]['Value']=[asb_values[i], asb_types[i]]
            self.core_dict['core'][cyc_tag][fuel_locations[i]]['Value']=[asb_values[i], asb_types[i]]

        self.full_core = self.get_full_core()
        self.core_lattice = self.get_full_lattice()
        return()

    def evaluate(self):
        """
        Redirects to the appropriate evaluate function

        Written by Gregory Delipei 7/29/2023
        """

        if self.ncycles==1:
            self.evaluate_cycle()
        else:
            self.evaluate_mcycle()
        return()
       
    def evaluate_cycle(self):
        """
        Creates the input deck, runs the calculation and retrieves the results and the cost.

        Single-Cycle problem.

        Parameters: 
        loc: String - Directory of execution
        fname: String - File name

        Written by Gregory Delipei 7/29/2023
        """
        start = time.time()
        print('Start of Single-Cycle Run')

        # Pre-processing steps for all cycles

        pwd = Path(os.getcwd())

        if not os.path.exists(self.name):
            os.makedirs(self.name)
        else:
            shutil.rmtree(self.name, ignore_errors=True)
            os.makedirs(self.name)
        
        if 'ml_model' in self.settings['genome']:
            RUN_ML=True
        else:
            RUN_ML=False

        cdir = self.library
        ncycle = int(self.settings['genome']['design_limits']['ncycle'])
        depfile = self.settings['genome']['design_limits']['depfile']
        os.chdir(self.name)
        shutil.copyfile(depfile, 'mcyc_exp.dep')
        fuel_locations = list(self.core_dict['fuel']['C'+str(ncycle)].keys())

        if 'inventory_binary' in self.settings['genome']['design_limits']:
            self.reload_inventory = pickle.load(open( self.settings['genome']['design_limits']['inventory_binary'], "rb" ) )
            if 'dinventory_binary' in self.settings['genome']['design_limits']:
                self.discharge_inventory = pickle.load(open( self.settings['genome']['design_limits']['dinventory_binary'], "rb" ) )
            else:
                self.discharge_inventory = {}
            self.cycle_parameters = {}
            self.reload_inventory_counter = 0
            for key in self.reload_inventory.keys():
                key_id = int(key[3:])
                if key_id > self.reload_inventory_counter:
                    self.reload_inventory_counter = key_id
            
            c_ass_groups = self.rank_assb()
            
            # same_asb = []

            # inv_keys=list(self.reload_inventory.keys())
            # for i in range(len(inv_keys)):
            #     for key,value in self.reload_inventory.items():
            #         if key != inv_keys[i] and value['BU3D'][5] == self.reload_inventory[inv_keys[i]]['BU3D'][5]:
            #             same_asb.append([inv_keys[i], key])

            self.lp_dict = {}
            self.lp_dict['C'+str(ncycle)]={}
            c_genome = self.genome
            c_lp = self.getlp_from_genome(c_genome,c_ass_groups,fuel_locations,'C'+str(ncycle))
            if c_lp is None:
                if ncycle==1:
                    self.cycle_parameters = {}
                self.set_poor_results(ncyc=ncycle)
                self.set_poor_results()
                reload_inv_path = 'inv_dts.p'
                discharge_inv_path = 'dinv_dts.p'
                pickle.dump( self.reload_inventory, open( reload_inv_path, "wb" ) )
                pickle.dump( self.discharge_inventory, open( discharge_inv_path, "wb" ) )
                os.system('rm -f mcyc_exp.dep')
                os.chdir(pwd)
                print('{} calculation is skipped'.format(self.name))
                print('End of Single-Cycle Run ... Duration: {} s'.format(time.time()-start))
                gc.collect()  
                print('finished collecting garbage...')
                print('exiting evaluate...')
                return()
                
            ## update inventory
            for i in range(len(c_lp)):
                key = c_lp[i]
                if key in self.reload_inventory:
                    self.reload_inventory_counter +=1 
                    new_value = copy.deepcopy(self.reload_inventory[key])
                    new_key = 'ASB'+str(self.reload_inventory_counter)
                    iloc = fuel_locations[i]
                    isym = len(self.core_dict['fuel']['C'+str(ncycle)][iloc]['Symmetric_Assemblies']) + 1
                    new_value['LOC'+str(ncycle)] = iloc
                    new_value['QTY'] = isym
                    self.reload_inventory[new_key] = new_value
                    c_lp[i] = new_key
                else:
                    pass

            cycle_dir = 'c' + str(ncycle)
            cyc_assemblies = []
            for i in range(len(c_lp)):
                iasb = c_lp[i]
                if iasb in self.reload_inventory:
                    itag = self.reload_inventory[iasb]['TYPE']
                else:
                    itag = int(iasb[2:])
                cyc_assemblies.append(itag)
            
            if RUN_ML:
                self.run_ml_cycle(fuel_locations,cyc_assemblies,c_lp,cycle_dir,ncyc=ncycle)
            else:
                self.run_cycle(fuel_locations,cyc_assemblies,c_lp,cycle_dir,ncyc=ncycle)

            if 'initial' in self.name and self.cycle_parameters['C'+str(ncycle)]["max_boron"] > 5000:
                print('Re-run initial case due to non-convergence')
                self.generate_initial(self.settings['genome']['chromosomes'])
                os.chdir(pwd)
                self.evaluate()
                return()
             
            self.clean_inventory()
            if 'design_limits' in self.settings['genome']:
                if 'reload_burnup' in self.settings['genome']['design_limits']:
                    max_bu = self.settings['genome']['design_limits']['reload_burnup']
                    self.update_discharge_inventory(max_bu)
        else:
            ## Get initial fuel assembly types and burnup
            cyc0_assemblies = []
            txt = Path('mcyc_exp.dep').read_text()
            txt_dep = txt.split(' Assembly Type Layout')[1].split('User')[0]
            txt_dep=txt_dep.split('\n')
            txt_dep=list(filter(lambda a: a != '', txt_dep))
            txt_dep=list(filter(lambda a: a != ' ', txt_dep))
            for i in range(len(txt_dep)):
                line_txt = txt_dep[i].split()
                for j in range(len(line_txt)):
                    val = int(line_txt[j])
                    if val != 10:
                        cyc0_assemblies.append(val)
            
            cyc0_assemblies = self.get_quarter_symmetry_values(cyc0_assemblies)
            nfa=len(cyc0_assemblies)
            bu_file = pwd / self.name
            bu_file = bu_file / 'mcyc_exp.dep'
            bu_2d, bu_3d=self.get_burnup(bu_file)

            ## Compute the reactivity of each fuel assembly and create initial reload inventory
            reac = self.get_reac(cyc0_assemblies,bu_2d)
            self.reload_inventory ={}
            self.discharge_inventory ={}
            self.reload_inventory_counter = 0
            for i in range(len(bu_2d)):
                idict = {}
                idict['BU2D']=bu_2d[i]
                idict['BU3D']=bu_3d[i,:]
                idict['REAC']=reac[i]
                idict['TYPE']=cyc0_assemblies[i]
                idict['QTY']= len(self.core_dict['fuel']['C'+str(ncycle)][fuel_locations[i]]['Symmetric_Assemblies']) + 1
                idict['LOC1']='OUT'
                idict['LOC2']='OUT'
                idict['LOC3']='OUT'
                self.reload_inventory_counter +=1 
                self.reload_inventory['ASB'+str(self.reload_inventory_counter)] = idict
            
            ## Group assemblies by reactivity

            c1_ass_groups = self.rank_assb()

            # Cycle 1 Calculation

            self.lp_dict = {}
            self.lp_dict['C'+str(ncycle)]={}
            c1_genome = self.genome
            c1_lp = self.getlp_from_genome(c1_genome,c1_ass_groups,fuel_locations,'C'+str(ncycle))
            if c1_lp is None:
                if ncycle==1:
                    self.cycle_parameters = {}
                self.set_poor_results(ncyc=ncycle)
                self.set_poor_results()
                reload_inv_path = 'inv_dts.p'
                discharge_inv_path = 'dinv_dts.p'
                pickle.dump( self.reload_inventory, open( reload_inv_path, "wb" ) )
                pickle.dump( self.discharge_inventory, open( discharge_inv_path, "wb" ) )
                os.system('rm -f mcyc_exp.dep')
                os.chdir(pwd)
                print('{} calculation is skipped'.format(self.name))
                print('End of Single-Cycle Run ... Duration: {} s'.format(time.time()-start))
                gc.collect()  
                print('finished collecting garbage...')
                print('exiting evaluate...')
                return()
                
            for i in range(len(c1_lp)):
                key = c1_lp[i]
                if key in self.reload_inventory:
                    self.reload_inventory_counter +=1 
                    new_value = copy.deepcopy(self.reload_inventory[key])
                    new_key = 'ASB'+str(self.reload_inventory_counter)
                    iloc = fuel_locations[i]
                    isym = len(self.core_dict['fuel']['C'+str(ncycle)][iloc]['Symmetric_Assemblies']) + 1
                    new_value['LOC'+str(ncycle)] = iloc
                    new_value['QTY'] = isym
                    self.reload_inventory[new_key] = new_value
                    c1_lp[i] = new_key
                else:
                    pass

           
            cycle_dir = 'c' + str(ncycle)
            cyc1_assemblies = []
            for i in range(len(c1_lp)):
                iasb = c1_lp[i]
                if iasb in self.reload_inventory:
                    itag = self.reload_inventory[iasb]['TYPE']
                else:
                    itag = int(iasb[2:])
                cyc1_assemblies.append(itag)

            if RUN_ML:
                self.run_ml_cycle(fuel_locations,cyc1_assemblies,c1_lp,cycle_dir,ncyc=ncycle)
            else:
                self.run_cycle(fuel_locations,cyc1_assemblies,c1_lp,cycle_dir,ncyc=ncycle)

            if 'initial' in self.name and self.cycle_parameters['C'+str(ncycle)]["max_boron"] > 5000:
                print('Re-run initial case due to non-convergence')
                self.generate_initial(self.settings['genome']['chromosomes'])
                os.chdir(pwd)
                self.evaluate()
                return()
             
            self.clean_inventory()
            if 'design_limits' in self.settings['genome']:
                if 'reload_burnup' in self.settings['genome']['design_limits']:
                    max_bu = self.settings['genome']['design_limits']['reload_burnup']
                    self.update_discharge_inventory(max_bu)

        # Get Single-Cycle Results

        self.parameters['max_boron']['value']=self.cycle_parameters['C'+str(ncycle)]['max_boron']
        self.parameters['PinPowerPeaking']['value']=self.cycle_parameters['C'+str(ncycle)]['PinPowerPeaking']
        self.parameters['FDeltaH']['value']=self.cycle_parameters['C'+str(ncycle)]['FDeltaH']
        self.parameters['cycle_length']['value'] =  self.cycle_parameters['C'+str(ncycle)]['cycle_length']    
        self.get_lcoe_cycle(ncycle)

        # Store Results 
        if 'options' in self.settings:
            if 'store' in self.settings['options']:
                opt = self.settings['options']['store']
                dts_fpath = "ldts.p" if opt=='light' else "hdts.p" 
                if dts_fpath in os.listdir('./c' + str(ncycle) + '/'):
                    c_dts = pickle.load(open( './c' + str(ncycle) + '/' + dts_fpath, "rb" ) )
                    rdict = {}
                    rdict['C'+str(ncycle)] = c_dts
                    for key, value in self.parameters.items():
                        rdict[key]= value['value']
                    pickle.dump( rdict, open( dts_fpath, "wb" ) )
                    os.system('rm -f {}'.format( './c' + str(ncycle) + '/' + dts_fpath))
        reload_inv_path = 'inv_dts.p'
        discharge_inv_path = 'dinv_dts.p'
        pickle.dump( self.reload_inventory, open( reload_inv_path, "wb" ) )
        pickle.dump( self.discharge_inventory, open( discharge_inv_path, "wb" ) )
        os.system('rm -f mcyc_exp.dep')
        os.chdir(pwd)
        print('{} calculation is done at {}!'.format(self.name,os.getcwd()))
        print('End of Single-Cycle Run ... Duration: {} s'.format(time.time()-start))
        gc.collect()  
        print('finished collecting garbage...')
        print('exiting evaluate...')
        return()

    def evaluate_mcycle(self):
        """
        Creates the input deck, runs the calculation and retrieves the results and the cost.

        Multi-Cycle problem with 3 cycles.

        Parameters: 
        loc: String - Directory of execution
        fname: String - File name

        Written by Gregory Delipei 7/29/2023
        """
        
        print('Start of Multi-Cycle Run')

        # Pre-processing steps for all cycles

        pwd = Path(os.getcwd())

        if not os.path.exists(self.name):
            os.makedirs(self.name)
        else:
            shutil.rmtree(self.name, ignore_errors=True)
            os.makedirs(self.name)
        
        cdir = self.library
        shutil.copyfile(cdir + '/' + 'mcyc_exp_quarter.dep', self.name +"/" + 'mcyc_exp.dep')
        os.chdir(self.name)

        ## Get initial fuel assembly types and burnup

        fuel_locations = list(self.core_dict['fuel']['C1'].keys())
        cyc0_assemblies = []
        txt = Path('mcyc_exp.dep').read_text()
        txt_dep = txt.split(' Assembly Type Layout')[1].split('========')[0]
        txt_dep=txt_dep.split('\n')
        txt_dep=list(filter(lambda a: a != '', txt_dep))
        txt_dep=list(filter(lambda a: a != ' ', txt_dep))
        for i in range(len(txt_dep)):
            line_txt = txt_dep[i].split()
            for j in range(len(line_txt)):
                val = int(line_txt[j])
                if val != 10:
                    cyc0_assemblies.append(val)

        nfa=len(cyc0_assemblies)
        bu_file = pwd / self.name
        bu_file = bu_file / 'mcyc_exp.dep'
        bu_2d, bu_3d=self.get_burnup(bu_file)

        ## Compute the reactivity of each fuel assembly and create initial reload inventory
        reac = self.get_reac(cyc0_assemblies,bu_2d)
        self.reload_inventory ={}
        self.discharge_inventory ={}
        self.reload_inventory_counter = 0
        for i in range(len(bu_2d)):
            idict = {}
            idict['BU2D']=bu_2d[i]
            idict['BU3D']=bu_3d[i,:]
            idict['REAC']=reac[i]
            idict['TYPE']=cyc0_assemblies[i]
            idict['LOC1']='OUT'
            idict['LOC2']='OUT'
            idict['LOC3']='OUT'
            self.reload_inventory_counter +=1 
            self.reload_inventory['ASB'+str(self.reload_inventory_counter)] = idict
            
        ## Group assemblies by reactivity

        c1_ass_groups = self.rank_assb()


        # Cycle 1 Calculation

        self.lp_dict = {}
        self.lp_dict['C1']={}
        c1_genome = self.genome[0:56]
        c1_lp = self.getlp_from_genome(c1_genome,c1_ass_groups)

        ## update inventory
        for key, value in self.reload_inventory.items():
            if key in c1_lp:
                id = c1_lp.index(key)
                iloc = fuel_locations[id]
                value['LOC1'] = iloc
            else:
                pass

        floc = 0
        for i in range(len(c1_lp)):
            cyc_tag = 'C1' 
            self.lp_dict[cyc_tag][fuel_locations[i]]=c1_lp[i]
            self.core_dict['fuel'][cyc_tag][fuel_locations[i]]['Value']=c1_lp[i]
            self.core_dict['core'][cyc_tag][fuel_locations[i]]['Value']=c1_lp[i]

        self.full_core = self.get_full_core()
        
        if self.map_size == 'quarter':
            self.core_lattice = self.get_quarter_lattice()
        else:
            self.core_lattice = self.get_full_lattice()
        
        cycle_dir = 'c1'
        cyc1_assemblies = []
        for i in range(len(c1_lp)):
            iasb = c1_lp[i]
            if iasb in self.reload_inventory:
                itag = self.reload_inventory[iasb]['TYPE']
            else:
                itag = int(iasb[2:])
            cyc1_assemblies.append(itag)

        self.run_cycle(fuel_locations,cyc1_assemblies,c1_lp,cycle_dir,ncyc=1)
        if 'initial' in self.name and self.cycle_parameters['C1']["max_boron"] > 5000:
            print('Re-run initial case due to non-convergence')
            self.generate_initial(self.settings['genome']['chromosomes'])
            os.chdir(pwd)
            self.evaluate()
            return()
        
        if 'design_limits' in self.settings['genome']:
            if 'reload_burnup' in self.settings['genome']['design_limits']:
                max_bu = self.settings['genome']['design_limits']['reload_burnup']
                self.update_discharge_inventory(max_bu)
        
        # Cycle 2 Calculation


        c2_ass_groups = self.rank_assb()

        self.lp_dict = {}
        self.lp_dict['C2']={}
        c2_genome = self.genome[56:112]
        c2_lp = self.getlp_from_genome(c2_genome,c2_ass_groups)
        ## update inventory
        for key, value in self.reload_inventory.items():
            if 'LOC2' not in value:
                value['LOC2']='OUT'
            if key in c2_lp:
                id = c2_lp.index(key)
                iloc = fuel_locations[id]
                value['LOC2'] = iloc
            else:
                pass

        floc = 0
        for i in range(len(c1_lp)):
            cyc_tag = 'C2' 
            self.lp_dict[cyc_tag][fuel_locations[i]]=c2_lp[i]
            self.core_dict['fuel'][cyc_tag][fuel_locations[i]]['Value']=c2_lp[i]
            self.core_dict['core'][cyc_tag][fuel_locations[i]]['Value']=c2_lp[i]
        
        self.full_core = self.get_full_core()
        
        if self.map_size == 'quarter':
            self.core_lattice = self.get_quarter_lattice()
        else:
            self.core_lattice = self.get_full_lattice()

        cycle_dir = 'c2'
        cyc2_assemblies = []
        for i in range(len(c2_lp)):
            iasb = c2_lp[i]
            if iasb in self.reload_inventory:
                itag = self.reload_inventory[iasb]['TYPE']
            else:
                itag = int(iasb[2:])
            cyc2_assemblies.append(itag)
        
        self.run_cycle(fuel_locations,cyc2_assemblies,c2_lp,cycle_dir,ncyc=2)
        if 'initial' in self.name and self.cycle_parameters['C2']["max_boron"] > 5000:
            print('Re-run initial case due to non-convergence')
            self.generate_initial(self.settings['genome']['chromosomes'])
            os.chdir(pwd)
            self.evaluate()
            return()

        if 'design_limits' in self.settings['genome']:
            if 'reload_burnup' in self.settings['genome']['design_limits']:
                max_bu = self.settings['genome']['design_limits']['reload_burnup']
                self.update_discharge_inventory(max_bu)

        # Cycle 3 Calculation

        
        c3_ass_groups = self.rank_assb()

        self.lp_dict = {}
        self.lp_dict['C3']={}
        c3_genome = self.genome[112:]
        c3_lp = self.getlp_from_genome(c3_genome,c3_ass_groups)
        for key, value in self.reload_inventory.items():
            if 'LOC1' not in value:
                value['LOC1']=None
            if 'LOC3' not in value:
                value['LOC3']='OUT'
            if key in c3_lp:
                id = c3_lp.index(key)
                iloc = fuel_locations[id]
                value['LOC3'] = iloc
            else:
                pass

        floc = 0
        for i in range(len(c1_lp)):
            cyc_tag = 'C3' 
            self.lp_dict[cyc_tag][fuel_locations[i]]=c3_lp[i]
            self.core_dict['fuel'][cyc_tag][fuel_locations[i]]['Value']=c3_lp[i]
            self.core_dict['core'][cyc_tag][fuel_locations[i]]['Value']=c3_lp[i]
        
        self.full_core = self.get_full_core()
        
        if self.map_size == 'quarter':
            self.core_lattice = self.get_quarter_lattice()
        else:
            self.core_lattice = self.get_full_lattice()

        cycle_dir = 'c3'
        cyc3_assemblies = []
        for i in range(len(c2_lp)):
            iasb = c3_lp[i]
            if iasb in self.reload_inventory:
                itag = self.reload_inventory[iasb]['TYPE']
            else:
                itag = int(iasb[2:])
            cyc3_assemblies.append(itag)

        self.run_cycle(fuel_locations,cyc3_assemblies,c3_lp,cycle_dir,ncyc=3)
        if 'initial' in self.name and self.cycle_parameters['C3']["max_boron"] > 5000:
            print('Re-run initial case due to non-convergence')
            self.generate_initial(self.settings['genome']['chromosomes'])
            os.chdir(pwd)
            self.evaluate()

        if 'design_limits' in self.settings['genome']:
            if 'reload_burnup' in self.settings['genome']['design_limits']:
                max_bu = self.settings['genome']['design_limits']['reload_burnup']
                self.update_discharge_inventory(max_bu)

        # Get MultiCycle Results
        
        self.parameters['max_boron']['value']=np.max(np.array([self.cycle_parameters['C1']['max_boron'],
                                                    self.cycle_parameters['C2']['max_boron'],
                                                    self.cycle_parameters['C3']['max_boron']]))
        self.parameters['PinPowerPeaking']['value']=np.max(np.array([self.cycle_parameters['C1']['PinPowerPeaking'],
                                                            self.cycle_parameters['C2']['PinPowerPeaking'],
                                                            self.cycle_parameters['C3']['PinPowerPeaking']]))
        self.parameters['FDeltaH']['value']=np.max(np.array([self.cycle_parameters['C1']['FDeltaH'],
                                                    self.cycle_parameters['C2']['FDeltaH'],
                                                    self.cycle_parameters['C3']['FDeltaH']]))
        self.parameters['cycle1_length']['value'] =  self.cycle_parameters['C1']['cycle_length']    
        self.parameters['cycle2_length']['value'] =  self.cycle_parameters['C2']['cycle_length']  
        self.parameters['cycle3_length']['value'] =  self.cycle_parameters['C3']['cycle_length']   
        self.get_lcoe()

        # Store Results 
        if 'options' in self.settings:
            if 'store' in self.settings['options']:
                opt = self.settings['options']['store']
                dts_fpath = "ldts.p" if opt=='light' else "hdts.p" 
                if dts_fpath in os.listdir('./c1/') and dts_fpath in os.listdir('./c2/') and dts_fpath in os.listdir('./c3/'):
                    c1_dts = pickle.load(open( './c1/' + dts_fpath, "rb" ) )
                    c2_dts = pickle.load(open( './c2/' + dts_fpath, "rb" ) )
                    c3_dts = pickle.load(open( './c3/' + dts_fpath, "rb" ) )
                    rdict = {}
                    rdict['C1'] = c1_dts
                    rdict['C2'] = c2_dts
                    rdict['C3'] = c3_dts
                    for key, value in self.parameters.items():
                        rdict[key]= value['value']
                    pickle.dump( rdict, open( dts_fpath, "wb" ) )
                    os.system('rm -f {}'.format('./c1/' + dts_fpath))
                    os.system('rm -f {}'.format('./c2/' + dts_fpath))
                    os.system('rm -f {}'.format('./c3/' + dts_fpath))
        os.system('rm -f mcyc_exp.dep')
        os.chdir(pwd)
        print('{} calculation is done at {}!'.format(self.name,os.getcwd()))
        print('End of Multi-Cycle Run')
        gc.collect()  
        print('finished collecting garbage...')
        print('exiting evaluate...')
        return()