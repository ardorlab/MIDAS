import os
import sys
import time
from turtle import color, pd
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
from sqlalchemy import false
import copy

class Simulate3_Core_17_17():
    """
    Class for the Simulate3 core loading pattern solutions.

    Parameters: Dictionary with the following structure
        'Symmetry': String specifying Quarter or Octant Symmetry
        'Symmetry_Axes': Tuple of 3 tuples indicating the 3 symmetry coordinates
        'Inventory': Dictionary with the available fuel assembly inventory
        'Inventory_Groups': Dictionary with the groups of each assembly.
    Written by Gregory Delipei 7/12/2022
    """
    def __init__(self,core_param):
        self.nrow=17
        self.ncol=17
        self.symmetry=core_param["Symmetry"]
        self.symmetry_axes=core_param["Symmetry_Axes"]
        self.core_dict={}
        self.core_dict['core_map'], self.core_dict['core_id'] = self.generate_core()
        self.core_dict['Inventory'] = core_param['Inventory']
        self.core_dict['Inventory_Groups'] = core_param['Inventory_Groups']
        self.random_design()

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
        core_map = ["R0000","R0001","R0002","R0003","R0004","R0005","R0006","R0007","R0008","R0009","R0010","R0011","R0012","R0013","R0014","R0015","R0016",
                    "R0100","R0101","R0102","R0103","R0104", "A05" , "A06" , "A07" , "A08" , "A09" , "A10" , "A11" ,"R0112","R0113","R0114","R0115","R0116",
                    "R0200","R0201","R0202", "B03" , "B04" , "B05" , "B06" , "B07" , "B08" , "B09" , "B10" , "B11" , "B12" , "B13" ,"R0214","R0215","R0216",
                    "R0300","R0301", "C02" , "C03" , "C04" , "C05" , "C06" , "C07" , "C08" , "C09" , "C10" , "C11" , "C12" , "C13" , "C14" ,"R0315","R0316",
                    "R0400","R0401", "D02" , "D03" , "D04" , "D05" , "D06" , "D07" , "D08" , "D09" , "D10" , "D11" , "D12" , "D13" , "D14" ,"R0415","R0416",
                    "R0500", "E01" , "E02" , "E03" , "E04" , "E05" , "E06" , "E07" , "E08" , "E09" , "E10" , "E11" , "E12" , "E13" , "E14" , "E15" ,"R0516",
                    "R0600", "F01" , "F02" , "F03" , "F04" , "F05" , "F06" , "F07" , "F08" , "F09" , "F10" , "F11" , "F12" , "F13" , "F14" , "F15" ,"R0616",
                    "R0700", "G01" , "G02" , "G03" , "G04" , "G05" , "G06" , "G07" , "G08" , "G09" , "G10" , "G11" , "G12" , "G13" , "G14" , "G15" ,"R0716",
                    "R0800", "H01" , "H02" , "H03" , "H04" , "H05" , "H06" , "H07" , "H08" , "H09" , "H10" , "H11" , "H12" , "H13" , "H14" , "H15" ,"R0816",
                    "R0900", "I01" , "I02" , "I03" , "I04" , "I05" , "I06" , "I07" , "I08" , "I09" , "I10" , "I11" , "I12" , "I13" , "I14" , "I15" ,"R0916",
                    "R1000", "J01" , "J02" , "J03" , "J04" , "J05" , "J06" , "J07" , "J08" , "J09" , "J10" , "J11" , "J12" , "J13" , "J14" , "J15" ,"R1016",
                    "R1100", "K01" , "K02" , "K03" , "K04" , "K05" , "K06" , "K07" , "K08" , "K09" , "K10" , "K11" , "K12" , "K13" , "K14" , "K15" ,"R1116",
                    "R1200","R1201", "L02" , "L03" , "L04" , "L05" , "L06" , "L07" , "L08" , "L09" , "L10" , "L11" , "L12" , "L13" , "L14" ,"R1215","R1216",
                    "R1300","R1301", "M02" , "M03" , "M04" , "M05" , "M06" , "M07" , "M08" , "M09" , "M10" , "M11" , "M12" , "M13" , "M14" ,"R1315","R1316",
                    "R1400","R1401","R1402", "N03" , "N04" , "N05" , "N06" , "N07" , "N08" , "N09" , "N10" , "N11" , "N12" , "N13" ,"R1414","R1415","R1416",
                    "R1500","R1501","R1502","R1503","R1504", "O05" , "O06" , "O07" , "O08" , "O09" , "O10" , "O11" ,"R1512","R1513","R1514","R1515","R1516",
                    "R1600","R1601","R1602","R1603","R1604","R1605","R1606","R1607","R1608","R1609","R1610","R1611","R1612","R1613","R1614","R1615","R1616"]
        core_map = np.array(core_map).reshape((self.nrow,self.ncol))
        core_id = []
        for i in range(self.nrow-1,-1,-1):
            for j in range(self.ncol):
                core_id.append((i-8,j-8))
        core_id=np.array(core_id).reshape((self.nrow,self.ncol,2))
        if self.symmetry == 'Quarter':
            self.core_dict['core'], self.core_dict['fuel'] = self.quarter_core(core_map,core_id)
        elif self.symmetry == 'Octant':
            self.core_dict['core'], self.core_dict['fuel'] = self.octant_core(core_map,core_id)
        else:
            raise ValueError(
                f"The selected symmetry ({self.symmetry}) is not valid."
            )
        return(core_map, core_id)
        
    def quarter_core(self,core_map,core_id):
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
           - core_dict: a core dictionary.

        Written by Gregory Delipei 7/12/2022
        """
        fuel_dict={}
        for key, value in core_dict.items():
            if key[0]=="R":
                continue
            else:
                fuel_dict[key]=value
        return(fuel_dict)

    def random_design(self):
        """
        Generates a random design following the constraints in the inventory.

        Parameters: None

        Written by Gregory Delipei 7/13/2022
        """

        # Initialization

        avail_locations=list(self.core_dict['fuel'].keys())
        assembly_types = list(self.core_dict['Inventory'].keys())

        for iass in assembly_types:
            self.core_dict['Inventory'][iass]['In_Desing'] = 0

        assembly_types_group = copy.deepcopy(assembly_types)
        avail_locations_group = copy.deepcopy(avail_locations)
        maxiter=10000 # maximum iterations to avoid infinite loop.
        nfuel, nfuel_sym, nrefl, nrefl_sym = self.compute_fa_number()
        total_fuel = nfuel
        total=0

        # Select randomly the fuel assemblies with exact limits in the inventory groups
        for key, value in self.core_dict['Inventory_Groups'].items():
            total_group=0
            if value['Limit']=='Exact':
                niter=0
                while total_group != value['Limit_Value'] and niter<maxiter:
                    avail_choices = []
                    for iloc in avail_locations_group:
                        for iass in value['Values']:
                            proposed_choice = (iloc, iass)
                            symmetry_multiplier = len(self.core_dict['fuel'][iloc]['Symmetric_Assemblies'])+1
                            if (total_group + symmetry_multiplier <= value['Limit_Value']) and (self.core_dict['Inventory'][iass]['In_Design'] <= self.core_dict['Inventory'][iass]['Max_Limit']-symmetry_multiplier):
                                avail_choices.append(proposed_choice)
                    # Re-iterate if the selected random filling strategy cannot meet the inventory limits.
                    if len(avail_choices)==0:
                        total_group=0
                        niter +=1
                        avail_locations_group = copy.deepcopy(avail_locations)
                        for iass in assembly_types_group:
                            self.core_dict['Inventory'][iass]['In_Design'] = 0
                        print(f"New Exact Filling Strategy - {niter}")
                        continue

                    sampled_choice = random.choice(avail_choices)
                    sloc = sampled_choice[0]
                    sass = sampled_choice[1]
                    self.core_dict['fuel'][sloc]['Value']=sass
                    self.core_dict['core'][sloc]['Value']=sass
                    symmetry_multiplier = len(self.core_dict['fuel'][sloc]['Symmetric_Assemblies'])+1
                    self.core_dict['Inventory'][sass]['In_Design']+= symmetry_multiplier
                    avail_locations_group.remove(sloc)
                    total_group+=symmetry_multiplier

                # Update the remaining available quantities for the next inventory groups
                for iass in value['Values']:
                    assembly_types_group.remove(iass)
                avail_locations = copy.deepcopy(avail_locations_group)
                assembly_types=copy.deepcopy(assembly_types_group)
                total+=total_group
        
        # Select randomly the fuel assemblies without exact limits in the inventory groups
        niter=0
        while total != total_fuel and niter<maxiter:
            avail_choices = []
            for iloc in avail_locations:
                for key, value in self.core_dict['Inventory_Groups'].items(): 
                    if value['Limit']!='Exact':
                        subtotal = 0
                        for iass in value["Values"]:
                            subtotal+=self.core_dict['Inventory'][iass]['In_Design']
                        for iass in value["Values"]:
                            proposed_choice = (iloc, iass)
                            symmetry_multiplier = len(self.core_dict['fuel'][iloc]['Symmetric_Assemblies'])+1
                            if (self.core_dict['Inventory'][iass]['In_Design'] <= self.core_dict['Inventory'][iass]['Max_Limit']-symmetry_multiplier) and (subtotal <=value['Limit_Value']-symmetry_multiplier):
                                avail_choices.append(proposed_choice)
      
            # Re-iterate if the selected random filling strategy cannot meet the inventory limits.
            if len(avail_choices)==0:
                total=total_group
                niter +=1
                avail_locations = copy.deepcopy(avail_locations_group)
                for iass in assembly_types:
                    self.core_dict['Inventory'][iass]['In_Desing'] = 0
                print(f"New Filling Strategy - {niter}")
                continue
            sampled_choice = random.choice(avail_choices)
            sloc = sampled_choice[0]
            sass = sampled_choice[1]
            self.core_dict['fuel'][sloc]['Value']=sass
            self.core_dict['core'][sloc]['Value']=sass
            symmetry_multiplier = len(self.core_dict['fuel'][sloc]['Symmetric_Assemblies'])+1
            self.core_dict['Inventory'][sass]['In_Design']+= symmetry_multiplier
            avail_locations.remove(sloc)
            total+=symmetry_multiplier

        return
    
    def action(self,act):
        """
        Performs an action on the current design and updates it.

        Parameters: 
            - act: Dictionary with the action options.

        Written by Gregory Delipei 7/14/2022
        """
        pass

        return

    def get_actions(self):
        """
        Extracts all possible actions in a dictionary.

        Parameters: None

        Written by Gregory Delipei 7/14/2022
        """

        # Compute exchange type of actions
        exchange_act = {}
        for key, value in self.core_dict["fuel"].items():
            avail_choices=[]
            key_symmetry = len(self.core_dict['fuel'][key]['Symmetric_Assemblies'])+1
            for key_ex,value_ex in  self.core_dict["fuel"].items():
                key_ex_symmetry = len(self.core_dict['fuel'][key_ex]['Symmetric_Assemblies'])+1
                if key_symmetry == key_ex_symmetry:
                    selected_choice=key_ex
                    avail_choices.append(selected_choice)
            exchange_act[key]=avail_choices

        # Compute change type of actions
        change_act = {}
        for key, value in self.core_dict["fuel"].items():
            avail_choices=[]
            loc_symmetry = len(self.core_dict['fuel'][key]['Symmetric_Assemblies'])+1
            loc_value = self.core_dict['fuel'][key]['Value']
            loc_group = self.get_inventory_group(loc_value)
            loc_group_limit = self.core_dict['Inventory_Groups'][loc_group]['Limit']
            for key_group,value_group in  self.core_dict['Inventory_Groups'].items(): 
                key_group_limit = self.core_dict['Inventory_Groups'][key_group]['Limit']
                for iass in value_group["Values"]:
                    if loc_group_limit=='Exact' and key_group==loc_group:
                        iass_limit = self.core_dict['Inventory'][iass]['Max_Limit']
                        iass_indesign = self.core_dict['Inventory'][iass]['In_Design']
                        if loc_symmetry + iass_indesign <= iass_limit:
                            selected_choice = iass
                            avail_choices.append(selected_choice)
                    elif loc_group_limit=='Max' and key_group_limit=='Max':
                        iass_limit = self.core_dict['Inventory'][iass]['Max_Limit']
                        iass_indesign = self.core_dict['Inventory'][iass]['In_Design']
                        group_limit = value_group['Limit_Value']
                        group_indesign = self.get_group_indesign(key_group)
                        if (loc_symmetry + iass_indesign <= iass_limit) and (loc_symmetry + group_indesign <= group_limit):
                            selected_choice = iass
                            avail_choices.append(selected_choice)
            change_act[key]=avail_choices

        act_dict={'Exchange': exchange_act,
                  'Change': change_act}
        return(act_dict)

    def get_inventory_group(self,iass):
        """
        Get in which group an assembly belongs to.

        Parameters: None

        Written by Gregory Delipei 7/14/2022
        """
        igroup = None
        for key,value in self.core_dict['Inventory_Groups'].items():
            if iass in value['Values']:
                igroup = key
        return(igroup)
    
    def get_group_indesign(self,group):
        """
        Get in which group an assembly belongs to.

        Parameters: None

        Written by Gregory Delipei 7/14/2022
        """
        sum_in = 0
        for iass in self.core_dict['Inventory_Groups'][group]['Values']:
            sum_in += self.core_dict['Inventory'][iass]['In_Design']
        return(sum_in)

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

    def compute_fa_number(self):
        """
        Computes the total number of fuel/reflector assemblies in the core with and without symmetry.

        Parameters: None

        Written by Gregory Delipei 7/14/2022
        """
        nfuel=0
        nfuel_sym=0
        nrefl = 0
        nrefl_sym = 0
        for key, value in self.core_dict['core'].items():
            symmetry_multiplier = len(self.core_dict['core'][key]['Symmetric_Assemblies'])+1
            if key[0]=='R':
                nrefl_sym+=1
                nrefl+=symmetry_multiplier
            else:
                nfuel_sym+=1
                nfuel+=symmetry_multiplier
        return(nfuel,nfuel_sym,nrefl,nrefl_sym)

        return


if __name__ == "__main__":
    """ Example of this class usage

     Define inventory information through nested dictionaries:
        Inventory: Dictionary with
                -keys the name of the fuel assembly type
                -values dictionaries with the following options:
                    -Max_Limit: Float for the maximum allowed number of this assembly in the design.
                    -In Design: Intiger for the number of this fuel assembly type in the current design.
                    -Cost: Float for the manufacturing cost of this assembly type.
                    -Tag: String for the associated tag for plotting or for creating code input decks.
        Inventory_Groups: Dictionary with
                -keys the name of the various inventory groups.
                -values dictionaries with the following options:
                    -Values: List of strings for the fuel assembly types belonging to this group.
                    -Limit: String for selecting this group to have an Exact limit to be reached in each design.
                    The Exact strings indicates this while the Max indicates a maximum limit constraint.
                    -Limit_Value: Float defining the limit value.   
        Core Parameters: Dictionary with
                -keys the various parameters requested by the class.
                    - Symmetry: String indicating the symmetry option. Quarter and Octant are available. 
                    - Symmetry_Axes: Three tuple coordinates indicating the symmetry axes. For Quarter 
                    symmetry the first coordinate is the center, the second coordinate is the limit in the 
                    horizontal direction, while the third coordinate is the limit on the vertical direction.
                    For Octant symmetry only the bottom-right octant is available so far so the coordinates are
                    forced to the center, diagonal right limit and vertical bottom limit.
                    - Inventory: Dictionary with the inventory information.
                    - Inventory_Groups: Dictionary with the inventory groups information.
                   
    """


    core_inventory={'FE500': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'E500'},
                    'FE450Gd8': {'Max_Limit':20, 'In_Design':0, 'Cost':0, 'Tag':'FE450Gd8'},
                    'FE450Bp8': {'Max_Limit':20, 'In_Design':0, 'Cost':0, 'Tag':'FE450Bp8'},
                    'FE340Bp4': {'Max_Limit':48, 'In_Design':0, 'Cost':0, 'Tag':'FE340Bp4'},
                    'FE340Gd4': {'Max_Limit':48, 'In_Design':0, 'Cost':0, 'Tag':'FE340Gd4'},
                    'RB1': {'Max_Limit':45, 'In_Design':0, 'Cost':0, 'Tag':'RB1'},
                    'RB2': {'Max_Limit':45, 'In_Design':0, 'Cost':0, 'Tag':'RB2'},
                    'RB3': {'Max_Limit':40, 'In_Design':0, 'Cost':0, 'Tag':'RB3'}
                    }
    core_inventory_groups={'Fresh': {'Values':['FE500','FE450Gd8','FE450Bp8','FE340Bp4','FE340Gd4'], 'Limit': 'Exact', 'Limit_Value':84},
                        'Reload': {'Values':['RB1','RB2','RB3'], 'Limit': 'Max', 'Limit_Value':200 }
                        }

    core_param={'Symmetry': 'Octant',
                'Symmetry_Axes': ((8,8),(16,16),(16,8)),
                'Inventory': core_inventory,
                'Inventory_Groups': core_inventory_groups}

    at = Simulate3_Core_17_17(core_param)
    at.plot_design('core_plot.png')

    sum_fresh=0
    for iass in core_inventory_groups['Fresh']['Values']:
        sum_fresh+=at.core_dict["Inventory"][iass]['In_Design']
    sum_reload=0
    for iass in core_inventory_groups['Reload']['Values']:
        sum_reload+=at.core_dict["Inventory"][iass]['In_Design']
    print(f"Sum of Fresh: {sum_fresh}")
    print(f"Sum of Reload: {sum_reload}")

    loc = 'O11'
    act=at.get_actions()
    sum_act=len(act['Exchange'][loc]) + len(act['Change'][loc])
    print(f"Sum of Actions for {loc}: {sum_act}")