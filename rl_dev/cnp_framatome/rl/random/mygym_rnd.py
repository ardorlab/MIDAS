import sys
import gym
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.monitor import Monitor
import numpy as np
from stable_baselines3 import SAC 
import time
import os
import random
sys.path.append('/home/gkdelipe/codes/mof/MOF/')

from rl_dev.cnp_framatome.games.gym_pwr_157 import Simulate3_Core_157 as PWR157

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
pwr_param = {   'lib':'../../lib',
                'Restart':'cycle1.res',
                'CASMO_XS':'cms.pwr-all.lib',
                'Thermal_Power': 100.0,
                'Core_Flow': 100.0,
                'Pressure': 2250.0,
                'Inlet_Temperature': 550.0 }

objectives = {'Cycle_Length':1.0,
                'Fdh': -400,
                'Fq': -400,
                'Max_Boron': -1}

start = 'state0.yml'
state_order=["O08","O09","N11","M12","L12","M11","N10","N09",
             "N08","M08","M09","M10","L11","K11","L10","L09",
             "L08","K08","K09","K10","J10","J09","J08","I08",
             "I09","H08"]

core_param={'Symmetry': 'Octant',
            'Symmetry_Axes': ((8,8),(16,16),(16,8)),
            'Inventory': core_inventory,
            'Inventory_Groups': core_inventory_groups,
            'Parameters':pwr_param,
            'Objectives': objectives,
            'Start':start,
            'Order': state_order}

log_dir = "./monitor_rnd/"
os.makedirs(log_dir, exist_ok=True)

env = PWR157(core_param)
with open(log_dir+'random.txt', "w") as ofile:
    ofile.write('Rewards\n')

for i in range(26*5000):
    action = random.uniform(-1,1)
    obs, reward, done, info = env.step(action)
    if done:
        with open(log_dir+'random.txt', "a") as ofile:
            ofile.write(str(reward) + "\n")
        env.reset()
