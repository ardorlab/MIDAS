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
from core_optimization_157 import Simulate3_Core_157



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

S3Game = Simulate3_Core_157(core_param)
    

# at.evaluate()
# print(at.core_dict['Results'])

state = {"H08": "FE200",
         "I08": "FE200",
         "I09": "FE200",
         "J08": "FE200",
         "J09": "FE200",
         "J10": "FE250Bp",
         "K08": "FE250Bp",
         "K09": "FE250",
         "K10": "FE250",
         "K11": "FE200",
         "L08": "FE200",
         "L09": "FE200",
         "L10": "FE200",
         "L11": "FE250",
         "L12": "FE200",
         "M08": "FE250",
         "M09": "FE250",
         "M10": "FE320Bp",
         "M11": "FE320Bp",
         "M12": "FE320",
         "N08": "FE320Bp",
         "N09": "FE320",
         "N10": "FE320",
         "N11": "FE320Bp",
         "O08": "FE200",
         "O09": "FE320Bp"}

S3Game.set_state(state)
S3Game.plot_design('core_plot.png')