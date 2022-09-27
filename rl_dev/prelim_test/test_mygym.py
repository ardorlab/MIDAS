from stable_baselines3.common.env_checker import check_env
from gym_core_optimization_157_2 import Simulate3_Core_157 as PWR157
import numpy as np
import random

core_inventory={'FE200': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'2'},
                'FE250': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'3'},
                'FE250Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'4'},
                'FE320': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'5'},
                'FE320Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'6'},
                }
core_inventory_groups={'E200': {'Values':['FE200'], 'Limit': 'Exact', 'Limit_Value':88},
                        'E250': {'Values':['FE250','FE250Bp'], 'Exact': 'Max', 'Limit_Value':56},
                        'E320': {'Values':['FE320','FE320Bp'], 'Exact': 'Max', 'Limit_Value':64} 
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

env = PWR157(core_param)
env.evaluate()
#check_env(env, warn=True)
obs = env.reset()
env.render()

print(env.observation_space)
print(env.action_space)
print(env.action_space.sample())

n_steps = 30
with open('mygym.txt', 'w') as f:
    for step in range(n_steps):
        f.write("Step {}".format(step + 1))
        obs, reward, done, info = env.step(random.uniform(-1,1))
        f.write('\n    obs=' + str(obs) +  '\n    reward=' + str(reward) + '\n    done=' + str(done) + "\n")
        env.render()
        if done:
            f.write("\nGoal reached!" + "\n    reward=" + str(reward))
            break