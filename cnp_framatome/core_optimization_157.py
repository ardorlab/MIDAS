import os
import sys
import time
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path
from simulate import Extractor

class Simulate3_Core_157():
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
        self.action_lower=-1
        self.action_upper=1
        self.symmetry=core_param["Symmetry"]
        self.symmetry_axes=core_param["Symmetry_Axes"]
        self.core_dict={}
        self.core_dict['core_map'], self.core_dict['core_id'] = self.generate_core()
        self.core_dict['Inventory'] = core_param['Inventory']
        self.core_dict['Inventory_Groups'] = core_param['Inventory_Groups']
        self.core_dict['Parameters'] = core_param['Parameters']
        self.core_dict['Objectives'] = core_param['Objectives']
        self.random_design()
        invent=list(self.core_dict['Inventory'].keys())
        nass = len(invent)
        self.cmap={}
        for i in range(nass):
            val= i-int(np.floor(nass/2))
            self.cmap[invent[i]]=val

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
        core_map = [  None ,  None ,  None ,  None ,  None ,  None ,"R0006","R0007","R0008","R0009","R0010",  None ,  None ,  None ,  None ,  None ,  None ,
                      None ,  None ,  None ,  None ,"R0104","R0105","R0106", "A07" , "A08" , "A09" ,"R0110","R0111","R0112",  None ,  None ,  None ,  None ,
                      None ,  None ,  None ,"R0203","R0204", "B05" , "B06" , "B07" , "B08" , "B09" , "B10" , "B11" ,"R0212","R0213",  None ,  None ,  None ,
                      None ,  None ,"R0302","R0303", "C04" , "C05" , "C06" , "C07" , "C08" , "C09" , "C10" , "C11" , "C12" ,"R0313","R0314",  None ,  None ,
                      None ,"R0401","R0402", "D03" , "D04" , "D05" , "D06" , "D07" , "D08" , "D09" , "D10" , "D11" , "D12" , "D13" ,"R0414","R0415",  None ,
                      None ,"R0501", "E02" , "E03" , "E04" , "E05" , "E06" , "E07" , "E08" , "E09" , "E10" , "E11" , "E12" , "E13" , "E14" ,"R0515",  None ,
                    "R0600","R0601", "F02" , "F03" , "F04" , "F05" , "F06" , "F07" , "F08" , "F09" , "F10" , "F11" , "F12" , "F13" , "F14" ,"R0615","R0616",
                    "R0700", "G01" , "G02" , "G03" , "G04" , "G05" , "G06" , "G07" , "G08" , "G09" , "G10" , "G11" , "G12" , "G13" , "G14" , "G15" ,"R0716",
                    "R0800", "H01" , "H02" , "H03" , "H04" , "H05" , "H06" , "H07" , "H08" , "H09" , "H10" , "H11" , "H12" , "H13" , "H14" , "H15" ,"R0816",
                    "R0900", "I01" , "I02" , "I03" , "I04" , "I05" , "I06" , "I07" , "I08" , "I09" , "I10" , "I11" , "I12" , "I13" , "I14" , "I15" ,"R0916",
                    "R1000","R1001", "J02" , "J03" , "J04" , "J05" , "J06" , "J07" , "J08" , "J09" , "J10" , "J11" , "J12" , "J13" , "J14" ,"R1015","R1016",
                      None ,"R1101", "K02" , "K03" , "K04" , "K05" , "K06" , "K07" , "K08" , "K09" , "K10" , "K11" , "K12" , "K13" , "K14" ,"R1115",  None ,
                      None ,"R1201","R1202", "L03" , "L04" , "L05" , "L06" , "L07" , "L08" , "L09" , "L10" , "L11" , "L12" , "L13" ,"R1214","R1215",  None ,
                      None ,  None ,"R1302","R1303", "M04" , "M05" , "M06" , "M07" , "M08" , "M09" , "M10" , "M11" , "M12" ,"R1313","R1314",  None ,  None ,
                      None ,  None ,  None ,"R1403","R1404", "N05" , "N06" , "N07" , "N08" , "N09" , "N10" , "N11" ,"R1412","R1413",  None ,  None ,  None ,
                      None ,  None ,  None ,  None ,"R1504","R1505","R1506", "O07" , "O08" , "O09" ,"R1510","R1511","R1512",  None ,  None ,  None ,  None ,
                      None ,  None ,  None ,  None ,  None ,  None ,"R1606","R1607","R1608","R1609","R1610",  None ,  None ,  None ,  None ,  None ,  None ]
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

    def get_full_core(self):
        """
        Generates the 17x17 full fuel core from symmetry.

        Parameters: None
    
        Written by Gregory Delipei 7/24/2022
        """
        full_core  = {}
        for key, value in self.core_dict['fuel'].items():
            full_core[key]=value['Value']
            for skey in value['Symmetric_Assemblies']:
               full_core[skey]=value['Value'] 
        return(full_core)

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
            if key is None:
                continue
            elif key[0]=="R":
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
        avail_actions = self.get_actions()
        action_type = act['Type']
        action_location = act['Location']
        action = act['Value']

        if action in avail_actions[action_type][action_location]:
            if action_type =='Exchange':
                loc_value = self.core_dict['core'][action_location]['Value']
                action_value = self.core_dict['core'][action]['Value']
                self.core_dict['core'][action_location]['Value'] = action_value
                self.core_dict['fuel'][action_location]['Value'] = action_value
                self.core_dict['core'][action]['Value'] = loc_value
                self.core_dict['fuel'][action]['Value'] = loc_value
            elif action_type=='Change':
                loc_value = self.core_dict['core'][action_location]['Value']
                loc_symmetry=len(self.core_dict['fuel'][action_location]['Symmetric_Assemblies'])+1
                self.core_dict['core'][action_location]['Value'] = action
                self.core_dict['fuel'][action_location]['Value'] = action
                self.core_dict['Inventory'][loc_value]['In_Design']-=loc_symmetry
                self.core_dict['Inventory'][action]['In_Design']+=loc_symmetry
        else:
            raise ValueError(
                f"The selected action is not valid."
            )
        return

    def mapaction(self,mact):
        """
        Performs an action on the current design and updates it.

        Parameters: 
            - act: Dictionary with the action options.

        Written by Gregory Delipei 7/14/2022
        """
        avail_actions = self.get_actions()
        action_mvalue = mact['Value']
        action_location = mact['Location']
        loc_actions=avail_actions['Map'][action_location]
        for key,value in loc_actions.items():
            bounds = value['Bounds']
            if bounds[0]<=action_mvalue<bounds[1]:
                action_type = value['Type']
                action = value['Value']
            if action_mvalue==bounds[1]==1:
                action_type = value['Type']
                action = value['Value']
        
        if action in avail_actions[action_type][action_location]:
            if action_type =='Exchange':
                loc_value = self.core_dict['core'][action_location]['Value']
                action_value = self.core_dict['core'][action]['Value']
                self.core_dict['core'][action_location]['Value'] = action_value
                self.core_dict['fuel'][action_location]['Value'] = action_value
                self.core_dict['core'][action]['Value'] = loc_value
                self.core_dict['fuel'][action]['Value'] = loc_value
            elif action_type=='Change':
                loc_value = self.core_dict['core'][action_location]['Value']
                loc_symmetry=len(self.core_dict['fuel'][action_location]['Symmetric_Assemblies'])+1
                self.core_dict['core'][action_location]['Value'] = action
                self.core_dict['fuel'][action_location]['Value'] = action
                self.core_dict['Inventory'][loc_value]['In_Design']-=loc_symmetry
                self.core_dict['Inventory'][action]['In_Design']+=loc_symmetry
        else:
            raise ValueError(
                f"The selected action is not valid."
            )
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
                    elif loc_group_limit=='Max' and key_group_limit=='Max' and key_group!=loc_group:
                        iass_limit = self.core_dict['Inventory'][iass]['Max_Limit']
                        iass_indesign = self.core_dict['Inventory'][iass]['In_Design']
                        group_limit = value_group['Limit_Value']
                        group_indesign = self.get_group_indesign(key_group)
                        if (loc_symmetry + iass_indesign <= iass_limit) and (loc_symmetry + group_indesign <= group_limit):
                            selected_choice = iass
                            avail_choices.append(selected_choice)
                    elif loc_group_limit=='Max' and key_group_limit=='Max' and key_group==loc_group:
                        iass_limit = self.core_dict['Inventory'][iass]['Max_Limit']
                        iass_indesign = self.core_dict['Inventory'][iass]['In_Design']
                        if (loc_symmetry + iass_indesign <= iass_limit):
                            selected_choice = iass
                            avail_choices.append(selected_choice)
            change_act[key]=avail_choices

        # Create mapping from [0,1] to action
        map_act = {}
        for key, value in self.core_dict["fuel"].items():
            nex_act = len(exchange_act[key])
            nch_act = len(change_act[key])
            nact = nex_act + nch_act
            act_bounds = np.linspace(self.action_lower,self.action_upper,nact+1)
            mdict={}
            it=0
            for i in range(nch_act):
                it+=1
                adict={}
                adict['Bounds'] = np.array([act_bounds[it-1],act_bounds[it]])
                adict['Type'] = 'Change'
                adict['Value'] = change_act[key][i]
                mdict['Act'+str(it)] = adict
            for i in range(nex_act):
                it+=1
                adict={}
                adict['Bounds'] = np.array([act_bounds[it-1],act_bounds[it]])
                adict['Type'] = 'Exchange'
                adict['Value'] = exchange_act[key][i]
                mdict['Act'+str(it)] = adict
            map_act[key]=mdict
                
        act_dict={'Exchange': exchange_act,
                  'Change': change_act,
                  'Map': map_act}
        return(act_dict)

    def get_mapstate(self):
        """
        Gets the current state in a normalized format.

        Parameters: None

        Written by Gregory Delipei 7/14/2022
        """
        mstate=np.zeros(len(self.core_dict['fuel'].keys()),dtype=np.int8)
        it=0
        for key, value in self.core_dict['fuel'].items():
            mstate[it]=self.cmap[value['Value']]
            it+=1
        return(mstate)

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
            if key is None:
                continue
            elif key[0]=='R':
                nrefl_sym+=1
                nrefl+=symmetry_multiplier
            else:
                nfuel_sym+=1
                nfuel+=symmetry_multiplier
        return(nfuel,nfuel_sym,nrefl,nrefl_sym)

        return

    def evaluate(self, ldir='./run', fname='solution'):
        """
        Creates the SIMULATE input deck, runs SIMULATE and retrieves the results and the cost.

        Parameters: 
        loc: String - Directory of execution
        fname: String - File name

        Written by Gregory Delipei 7/29/2022
        """

        full_core = self.get_full_core()

        # Create SIMULATE INPUT DECK

        pwd = Path(os.getcwd())
        if not os.path.exists(ldir):
            os.makedirs(ldir)
        else:
            shutil.rmtree(ldir, ignore_errors=True)
            os.makedirs(ldir)
        
        shutil.copyfile(self.core_dict['Parameters']['Restart'], ldir +"/" + self.core_dict['Parameters']['Restart'])
        shutil.copyfile(self.core_dict['Parameters']['CASMO_XS'], ldir +"/" + self.core_dict['Parameters']['CASMO_XS'])
        os.chdir(ldir)
        filename = fname+'.inp'
        with open(filename, 'a') as ofile:
            str_res = "\'RES\' " + "\'" + self.core_dict['Parameters']['Restart'] + "\' " +  "0.0/ \n"
            ofile.write(str_res)
            ofile.write("\n")
            
            irow = int(self.nrow/2)
            icol = int(self.ncol/2)
            for i in range(int(np.ceil(self.nrow/2))):
                str_lp = "\'FUE.TYP\'  " + str(i+1) + ", " 
                for j in range(int(np.ceil(self.nrow/2))):
                    loc=self.core_dict['core_map'][irow+i,icol+j]
                    if loc is None:
                        tag = "0"
                    elif loc[0]=="R":
                        tag = "1"
                    else:
                        val = full_core[loc]
                        tag = self.core_dict['Inventory'][val]['Tag']
                    str_lp = str_lp + tag
                    if(j==int(np.ceil(self.nrow/2)-1)):
                         str_lp = str_lp + '/\n'
                    else:
                        str_lp = str_lp + " "
                ofile.write(str_lp)

            ofile.write("\n")
            ofile.write("\n")
            str_cope = f"\'COR.OPE\' " + str(self.core_dict['Parameters']['Thermal_Power']) + ", " + str(self.core_dict['Parameters']['Core_Flow'])+ ", "  + str(self.core_dict['Parameters']['Pressure']) + "/ \n" 
            ofile.write(str_cope)
            str_ctin = f"\'COR.TIN\' " + str(self.core_dict['Parameters']['Inlet_Temperature']) + "/ \n" 
            ofile.write(str_ctin)
            ofile.write("\n")

            str_depsta = "\'DEP.STA\' \'AVE\' 0.0 0.5 1 2 -1 20 / \n"
            ofile.write(str_depsta)
            str_com = "\'COM\' The following performs an automated search on cycle length at a boron concentration of 10 ppm \n"
            ofile.write(str_com)
            str_itsrc= "\'ITE.SRC\' \'SET\' \'EOLEXP\' , , 0.02, , , , , , \'MINBOR\' 10., , , , , 4, 4, , , / \n"
            ofile.write(str_itsrc)
            str_sta = "\'STA\'/ \n"
            ofile.write(str_sta)
            str_end = "\'END\'/ \n"
            ofile.write(str_end)

        # Run SIMULATE INPUT DECK
        sim3_maxruns = 1000
        s3cmd = "simulate3 " + filename
        for i in range(sim3_maxruns):
            try:
                os.system(s3cmd)
                # Get Results
                ofile = fname + '.out'
                res=self.get_sim3_results(ofile)
                self.core_dict["Results"] = res
                fit=self.get_fitness(self.core_dict["Results"])
                self.core_dict["Results"]["Fitness"] = fit
                break
            except:
                print(f"Failed SIMULATE3 run {i+1}")
                os.remove(ofile)

        os.chdir(pwd)

    def get_sim3_results(self,ofile):
        res = {}
        file_ = open(ofile)
        file_lines = file_.readlines()
        file_.close()

        EFPD_list = Extractor.efpd_list(file_lines)
        res["Cycle_Length"] = EFPD_list[-1]
        
        FDH_list = Extractor.FDH_list(file_lines)
        res["Fdh"] = max(FDH_list)
           
        peak_list = Extractor.pin_peaking_list(file_lines)
        res["Fq"] = max(peak_list)
            
        boron_list = Extractor.boron_list(file_lines)
        res["Max_Boron"] = max(boron_list)

        return(res)
    
    def get_fitness(self,res):
        
        fit = 0 
        fit += res["Cycle_Length"]*self.core_dict['Objectives']["Cycle_Length"]
        fit += max(0,res["Max_Boron"]-1300)*self.core_dict['Objectives']["Max_Boron"]
        fit += max(0,res["Fdh"]-1.48)*self.core_dict['Objectives']["Fdh"]
        fit += max(0,res["Fq"]-2.10)*self.core_dict['Objectives']["Fq"]
    
        return(fit)
    
    def set_state(self,state):
        
        for key,value in self.core_dict['Inventory'].items():
            nvalue = value['In_Design']=0
            self.core_dict['Inventory'][key] = value

        for key, value in state.items():
            symmetry_multiplier = len(self.core_dict['fuel'][key]['Symmetric_Assemblies'])+1
            self.core_dict['fuel'][key]['Value']=value
            self.core_dict['core'][key]['Value']=value
            self.core_dict['Inventory'][value]['In_Design']+=symmetry_multiplier
        return

    def get_state(self):
        state={}
        for key, value in self.core_dict['fuel'].items():
            state[key]=value['Value']
        return(state)
    

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

