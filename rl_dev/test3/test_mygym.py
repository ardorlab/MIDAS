from stable_baselines3.common.env_checker import check_env
from gym_core_optimization_157_4 import Simulate3_Core_157 as PWR157
import numpy as np
import random

core_inventory={'FE200': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'2'},
                'FE250': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'3'},
                'FE250Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'4'},
                'FE320': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'5'},
                'FE320Bp': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'6'},
                }
core_inventory_groups={'E200': {'Values':['FE200'], 'Limit': 'Exact', 'Limit_Value':11},
                        'E250': {'Values':['FE250','FE250Bp'], 'Limit': 'Exact', 'Limit_Value':7},
                        'E320': {'Values':['FE320','FE320Bp'], 'Limit': 'Exact', 'Limit_Value':8} 
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

for i in ["a","b","c"]:
    for j in [1,2,3]:
        p=j
        if j==2:
            continue
    
# env = PWR157(core_param)
# env.evaluate()
# #check_env(env, warn=True)
# obs = env.reset()
# env.render()

# print(env.observation_space)
# print(env.action_space)
# print(env.action_space.sample())

# n_steps = 30
# # with open('mygym.txt', 'w') as f:
# #     for step in range(n_steps):
# #         f.write("Step {}".format(step + 1))
# #         obs, reward, done, info = env.step(random.uniform(-1,1))
# #         f.write('\n    obs=' + str(obs) +  '\n    reward=' + str(reward) + '\n    done=' + str(done) + "\n")
# #         env.render()
# #         if done:
# #             f.write("\nGoal reached!" + "\n    reward=" + str(reward))
# #             break

at = PWR157(core_param)
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

loc = 'L09'
sact={'Type': 'Change',
        'Location': loc,
        'Value': act['Change'][loc][0]}

loc_value = at.core_dict['core'][loc]['Value']
at.action(sact)
at.plot_design('core_plot_act2.png')



mact = {'Location':'J09',
        'Value':0.6}
at.mapaction(mact)
at.plot_design('core_plot_act3.png')