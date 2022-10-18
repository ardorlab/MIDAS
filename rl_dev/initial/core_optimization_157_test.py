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
import core_optimization_157  


core_inventory={'FE200': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'2'},
                'FE250': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'3'},
                'FE250Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'4'},
                'FE320': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'5'},
                'FE320Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'6'},
                }
core_inventory_groups={'E200': {'Values':['FE200'], 'Limit': 'Max', 'Limit_Value':88},
                        'E250': {'Values':['FE250','FE250Bp'], 'Limit': 'Max', 'Limit_Value':56},
                        'E320': {'Values':['FE320','FE320Bp'], 'Limit': 'Max', 'Limit_Value':64} 
                    }
pwr_param = {'Restart':'cycle1.res',
                'CASMO_XS':'cms.pwr-all.lib',
                'Thermal_Power': 100.0,
                'Core_Flow': 100.0,
                'Pressure': 2250.0,
                'Inlet_Temperature': 550.0 }

objectives = {'Cycle_Length':1.0,
                'Fdh': -400,
                'Fq': -400,
                'Max_Boron': -1}

core_param={'Symmetry': 'Octant',
            'Symmetry_Axes': ((8,8),(16,16),(16,8)),
            'Inventory': core_inventory,
            'Inventory_Groups': core_inventory_groups,
            'Parameters':pwr_param,
            'Objectives': objectives}

at = Simulate3_Core_157(core_param)
at.plot_design('core_plot.png')

sum_e200=0
for iass in core_inventory_groups['E200']['Values']:
    sum_e200+=at.core_dict["Inventory"][iass]['In_Design']
sum_e250=0
for iass in core_inventory_groups['E250']['Values']:
    sum_e250+=at.core_dict["Inventory"][iass]['In_Design']
sum_e320=0
for iass in core_inventory_groups['E320']['Values']:
    sum_e320+=at.core_dict["Inventory"][iass]['In_Design']
print(f"Sum of E200: {sum_e200}")
print(f"Sum of E250: {sum_e250}")
print(f"Sum of E320: {sum_e320}")

loc = 'N11'
act=at.get_actions()
sum_act=len(act['Exchange'][loc]) + len(act['Change'][loc])
print(f"Sum of Actions for {loc}: {sum_act}")

loc = 'K10'
sact={'Type': 'Change',
        'Location': loc,
        'Value': act['Change'][loc][0]}

loc_value = at.core_dict['core'][loc]['Value']
print(f"Location In Design before action: {at.core_dict['Inventory'][loc_value]['In_Design']}")
print(f"Action In Design before action: {at.core_dict['Inventory'][act['Change'][loc][0]]['In_Design']}")
at.action(sact)
at.plot_design('core_plot_act.png')
print(f"Location In Design after action: {at.core_dict['Inventory'][loc_value]['In_Design']}")
print(f"Action In Design after action: {at.core_dict['Inventory'][act['Change'][loc][0]]['In_Design']}")
act=at.get_actions()
sum_e200=0
for iass in core_inventory_groups['E200']['Values']:
    sum_e200+=at.core_dict["Inventory"][iass]['In_Design']
sum_e250=0
for iass in core_inventory_groups['E250']['Values']:
    sum_e250+=at.core_dict["Inventory"][iass]['In_Design']
sum_e320=0
for iass in core_inventory_groups['E320']['Values']:
    sum_e320+=at.core_dict["Inventory"][iass]['In_Design']
print(f"Sum of E200: {sum_e200}")
print(f"Sum of E250: {sum_e250}")
print(f"Sum of E320: {sum_e320}")
