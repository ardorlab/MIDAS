import os
import sys
import time
from turtle import update
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path
from simulate import Extractor
from core_optimization_157 import Simulate3_Core_157
import yaml


def OutInM1(S3Game,res_dict,gameid=0):
    S0 = copy.deepcopy(S3Game)
    state_order=["O08","O09","N11","M12","L12","M11","N10","N09",
                 "N08","M08","M09","M10","L11","K11","L10","L09",
                 "L08","K08","K09","K10","J10","J09","J08","I08",
                 "I09","H08"]
    for i in range(len(state_order)):
        loc = state_order[i]
        act=S0.get_actions()
        act_rnd = random.uniform(-1,1)
        action={'Location': loc,
                'Value': act_rnd}
        S0.mapaction(action)
        S0.evaluate()
        res_dict['G'+str(gameid) + 'M' + str(i+1)]= copy.deepcopy(S0.core_dict)
    return(res_dict)

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
state0_file = 'state0.yml'
with open(state0_file) as f:
    state0 = yaml.safe_load(f)
S3Game.set_state(state0)
S3Game.evaluate()
gameid=0
moveid=0
res_dict={'G'+str(gameid) + 'M' + str(moveid): copy.deepcopy(S3Game.core_dict)}

# t0=time.time()
# updated_res=OutInM1(S3Game,res_dict,gameid=1)
# t1=time.time()
# print(f'Game time = {t1-t0}s')
# best_key = 'G0M0'
# best_fit = updated_res["G0M0"]["Results"]['Fitness'] 
# for key,value in updated_res.items():
#     if value["Results"]['Fitness'] > best_fit:
#         best_fit = value["Results"]['Fitness']
#         best_key = key
# print("The best loading pattern is: \n")
# print(updated_res[best_key]["Results"])