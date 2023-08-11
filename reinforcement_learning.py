import os
import sys
import time
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import shutil
import copy
import fitness
from pathlib import Path
from metrics import Optimization_Metric_Toolbox
import gym
from gym import spaces
import yaml
from fileinput import filename
import gym
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.monitor import Monitor
from stable_baselines3 import SAC, PPO, A2C, DQN
import torch as th

class Cycle1_Gym_Env(gym.Env):
    """
    Class for wrapper adapted for Gym environment

    Written by Gregory Delipei 7/12/2022
    """
    metadata = {'render.modes': ['console']}
    # Define constants for clearer code

    def __init__(self,solution,file_settings,fitness):
        super(Gym_Env, self).__init__()
        self.solution = solution
        self.best_solution = solution
        self.nrow=solution.nrow
        self.fitness= fitness
        self.ncol=solution.ncol
        self.action_lower=-1
        self.action_upper=1
        self.symmetry=solution.symmetry
        self.core_dict={}
        self.core_dict['core_map'] = solution.core_dict['core_map']
        self.core_dict['core_id'] = solution.core_dict['core_id']
        self.core_dict['Inventory'] = solution.core_dict['Inventory']
        self.order = file_settings['optimization']['order']
        self.restart = file_settings['optimization']['restart']
        self.fitness_const = file_settings['optimization']['selection']['fitness_constraint']
        state0_file = file_settings['optimization']['start']
        
      #  self.random_design()
        invent=list(self.core_dict['Inventory'].keys())
        nass = len(invent)
        self.counter=0
        self.total_run=0
        self.maxcounter = len(self.order)
        with open(state0_file) as f:
            self.start = yaml.safe_load(f)
        self.solution.set_state(self.start)
        self.action_type= file_settings['optimization']['stable_baselines3_options']['action_space']
        self.observation_type = file_settings['optimization']['stable_baselines3_options']['observation_space']
        if self.action_type=='discrete':
            self.cmap={}
            cmap_range=np.arange(0,nass)
            for i in range(nass):
                self.cmap[invent[i]]=cmap_range[i]
            self.action_space = spaces.Discrete(nass)
        else:
            self.cmap={}
            cmap_range=np.linspace(-1,1,nass)
            for i in range(nass):
                self.cmap[invent[i]]=cmap_range[i]
            self.action_space = spaces.Box(low=-1, high=1,
                                            shape=(1,), dtype=np.float32)
        
        if self.observation_type=='multi_discrete':
            self.observation_space = spaces.MultiDiscrete([nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,31])
        else:
            self.observation_space = spaces.Box(low=-1.0, high=1.0, shape=(2,len(self.start)), dtype=np.float32)
        

        self.trackfile='trackfile.txt'
        with open(self.trackfile, "w") as ofile:
            ofile.write('PWR core optimization with MOF')

    def reset(self):
        """
        Important: the observation must be a numpy array
        :return: (np.array) 
        """
        self.counter=0
        if self.restart:
            new_state = {}
            for key, value in self.best_solution.core_dict['fuel'].items():
                new_state[key]=value['Value']
            with open('state1.yml', 'w') as outfile:
                yaml.dump(new_state, outfile, default_flow_style=False)
            self.solution.set_state(new_state)
        else:
            self.solution.set_state(self.start)
        
        if self.observation_type=="continuous":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            loc_state = np.zeros(cstate.shape[0])
            loc_state[self.counter]=1
            rstate=np.r_[[cstate],[loc_state]]
        elif self.observation_type=="multi_discrete":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            cstate[-1] = self.counter
            rstate=cstate
        return(rstate)

    def step(self, action):
        loc = self.order[self.counter]
        act=self.solution.get_actions()
        sact={'Location': loc,
                'Value': action,
                'Space':self.action_type,
                'Action_Map':self.cmap}
        self.solution.mapaction(sact)

        if self.counter==len(self.start)-1:
            mid = 0
            self.total_run +=1
            self.solution.name = "solution_{}".format(self.total_run)
            self.solution.evaluate()
            solList = self.fitness.calculate([self.solution])
            self.solution=solList[0]
            if self.fitness_const:
                self.solution.fitness=self.solution.get_fitness()
            reward = self.solution.fitness
            info={'cycle_length':self.solution.parameters["cycle_length"]['value'],
                  'max_boron':self.solution.parameters["max_boron"]['value'],
                  'FDeltaH':self.solution.parameters["FDeltaH"]['value'],
                  'PinPowerPeaking':self.solution.parameters["PinPowerPeaking"]['value'],
                  'State': self.solution.get_mapstate(self.cmap,self.observation_type)}
            all_values = open('all_value_tracker.txt','a')
            all_values.write(f"{self.total_run},   {self.solution.name},   ")
            for param in self.solution.parameters:
                all_values.write(f"{self.solution.parameters[param]['value']},    ")
            all_values.write(f"{self.solution.fitness},    ")
            all_values.write('\n')
            all_values.close()
            if self.total_run ==1:
                self.best_solution=copy.deepcopy(self.solution)
            else:
                if self.best_solution.fitness<self.solution.fitness:
                    self.best_solution = copy.deepcopy(self.solution)
        else:
            mid = self.counter + 1
            reward=0
            info={}
        self.counter+=1
        done = bool(self.counter==len(self.order))
        
        if self.observation_type=="continuous":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            loc_state = np.zeros(cstate.shape[0])
            loc_state[mid]=1
            rstate=np.r_[[cstate],[loc_state]]
        elif self.observation_type=="multi_discrete":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            cstate[-1] = mid
            rstate=cstate
        return((rstate, reward, done, info))
    
    def render(self, mode='console'):
        if mode != 'console':
            mode='console'
        with open(self.trackfile, 'a') as ofile:
            ofile.write("\nStep {}".format(self.counter))
            ofile.write('\n    obs=' + str(self.solution.get_mapstate()) +  '\n    reward=' + str(self.solution.fitness) + "\n")

    def close(self):
        pass

class Mcycle_Gym_Env(gym.Env):
    """
    Class for wrapper adapted for Gym environment

    Written by Gregory Delipei 7/12/2022
    """
    metadata = {'render.modes': ['console']}
    # Define constants for clearer code

    def __init__(self,solution,file_settings,fitness):
        super(Gym_Env, self).__init__()
        self.solution = solution
        self.best_solution = solution
        self.nrow=solution.nrow
        self.fitness= fitness
        self.ncol=solution.ncol
        self.action_lower=-1
        self.action_upper=1
        self.symmetry=solution.symmetry
        self.core_dict={}
        self.core_dict['core_map'] = solution.core_dict['core_map']
        self.core_dict['core_id'] = solution.core_dict['core_id']
        self.core_dict['Inventory'] = solution.core_dict['Inventory']
        self.order = file_settings['optimization']['order']
        self.restart = file_settings['optimization']['restart']
        self.fitness_const = file_settings['optimization']['selection']['fitness_constraint']
        state0_file = file_settings['optimization']['start']
        
      #  self.random_design()
        invent=list(self.core_dict['Inventory'].keys())
        nass = len(invent)
        self.counter=0
        self.total_run=0
        self.maxcounter = len(self.order)
        with open(state0_file) as f:
            self.start = yaml.safe_load(f)
        self.solution.set_state(self.start)
        self.action_type= file_settings['optimization']['stable_baselines3_options']['action_space']
        self.observation_type = file_settings['optimization']['stable_baselines3_options']['observation_space']
        if self.action_type=='discrete':
            self.cmap={}
            cmap_range=np.arange(0,nass)
            for i in range(nass):
                self.cmap[invent[i]]=cmap_range[i]
            self.action_space = spaces.Discrete(nass)
        else:
            self.cmap={}
            cmap_range=np.linspace(-1,1,nass)
            for i in range(nass):
                self.cmap[invent[i]]=cmap_range[i]
            self.action_space = spaces.Box(low=-1, high=1,
                                            shape=(1,), dtype=np.float32)
        
        if self.observation_type=='multi_discrete':
            self.observation_space = spaces.MultiDiscrete([nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,nass,nass,nass,nass,nass,
                                                      nass,31])
        else:
            self.observation_space = spaces.Box(low=-1.0, high=1.0, shape=(2,len(self.start)), dtype=np.float32)
        

        self.trackfile='trackfile.txt'
        with open(self.trackfile, "w") as ofile:
            ofile.write('PWR core optimization with MOF')

    def reset(self):
        """
        Important: the observation must be a numpy array
        :return: (np.array) 
        """
        self.counter=0
        if self.restart:
            new_state = {}
            for key, value in self.best_solution.core_dict['fuel'].items():
                new_state[key]=value['Value']
            with open('state1.yml', 'w') as outfile:
                yaml.dump(new_state, outfile, default_flow_style=False)
            self.solution.set_state(new_state)
        else:
            self.solution.set_state(self.start)
        
        if self.observation_type=="continuous":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            loc_state = np.zeros(cstate.shape[0])
            loc_state[self.counter]=1
            rstate=np.r_[[cstate],[loc_state]]
        elif self.observation_type=="multi_discrete":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            cstate[-1] = self.counter
            rstate=cstate
        return(rstate)

    def step(self, action):
        loc = self.order[self.counter]
        act=self.solution.get_actions()
        sact={'Location': loc,
                'Value': action,
                'Space':self.action_type,
                'Action_Map':self.cmap}
        self.solution.mapaction(sact)

        if self.counter==len(self.start)-1:
            mid = 0
            self.total_run +=1
            self.solution.name = "solution_{}".format(self.total_run)
            self.solution.evaluate()
            solList = self.fitness.calculate([self.solution])
            self.solution=solList[0]
            if self.fitness_const:
                self.solution.fitness=self.solution.get_fitness()
            reward = self.solution.fitness
            info={'cycle_length':self.solution.parameters["cycle_length"]['value'],
                  'max_boron':self.solution.parameters["max_boron"]['value'],
                  'FDeltaH':self.solution.parameters["FDeltaH"]['value'],
                  'PinPowerPeaking':self.solution.parameters["PinPowerPeaking"]['value'],
                  'State': self.solution.get_mapstate(self.cmap,self.observation_type)}
            all_values = open('all_value_tracker.txt','a')
            all_values.write(f"{self.total_run},   {self.solution.name},   ")
            for param in self.solution.parameters:
                all_values.write(f"{self.solution.parameters[param]['value']},    ")
            all_values.write(f"{self.solution.fitness},    ")
            all_values.write('\n')
            all_values.close()
            if self.total_run ==1:
                self.best_solution=copy.deepcopy(self.solution)
            else:
                if self.best_solution.fitness<self.solution.fitness:
                    self.best_solution = copy.deepcopy(self.solution)
        else:
            mid = self.counter + 1
            reward=0
            info={}
        self.counter+=1
        done = bool(self.counter==len(self.order))
        
        if self.observation_type=="continuous":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            loc_state = np.zeros(cstate.shape[0])
            loc_state[mid]=1
            rstate=np.r_[[cstate],[loc_state]]
        elif self.observation_type=="multi_discrete":
            cstate = self.solution.get_mapstate(self.cmap,self.observation_type)
            cstate[-1] = mid
            rstate=cstate
        return((rstate, reward, done, info))
    
    def render(self, mode='console'):
        if mode != 'console':
            mode='console'
        with open(self.trackfile, 'a') as ofile:
            ofile.write("\nStep {}".format(self.counter))
            ofile.write('\n    obs=' + str(self.solution.get_mapstate()) +  '\n    reward=' + str(self.solution.fitness) + "\n")

    def close(self):
        pass


class Reinforcement_Learning(object):
    """
    Class for performing optimization through reinforcement learning 
    using stable_baselines3

    Parameters:
        solution: Class
            Contains the genome, fitness, and optimization
            objective scores for solutions to the optimization problem.
        population: Class
            Class that contains the population size and stores the current
            solutions in the parent and child populations.
        generation: Class
            Keeps track of the current and total number of generations that
        file_settings: Dictionary
            The settings file read into the optimization. Carried through because
            some information needed to be carried along, but pickling all of the
            information didn't seem like a good way to carrty it thorugh the optimization.
    """

    def __init__(self, solution,
                 population,
                 generation,
                 fitness,
                 num_procs,
                 file_settings):

        self.solution = solution
        self.population = population
        self.generation = generation
        self.fitness= fitness
        self.num_procs = num_procs
        self.file_settings = file_settings


    def main_in_serial(self):
        """
        Performs optimization using a reinforcement learning in serial.
        Done because for some reason the neural network stuff seems to be
        breaking with parallel.

        Parameters: None

        Written by Brian Andersen 1/9/2020
        """
        opt = Optimization_Metric_Toolbox()

        track_file = open('optimization_track_file.txt', 'w')
        track_file.write("Beginning Optimization \n")
        track_file.close()

        all_values = open('all_value_tracker.txt','w')
        all_values.close()


        loading_pattern_tracker = open("loading patterns.txt", 'w')
        loading_pattern_tracker.close()

        log_dir = self.file_settings['optimization']['stable_baselines3_options']['logdir']
        os.makedirs(log_dir, exist_ok=True)

        objectves = self.file_settings['optimization']['objectives']
        
        info_kwd=list(objectves.keys())
        info_kwd.append('State')
        info_kwd = tuple(info_kwd)
        foo = self.solution()
        foo.name = "solution"
        foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        foo.add_additional_information(self.file_settings)
        Custom_Env = globals()[self.file_settings['optimization']['environment']]
        env = Monitor(Custom_Env(foo,self.file_settings,self.fitness),log_dir,info_keywords=info_kwd)
        env = DummyVecEnv([lambda: env])
        net1 = self.file_settings['optimization']['stable_baselines3_options']['policy_net']
        net2 = self.file_settings['optimization']['stable_baselines3_options']['qvalue_net']
        tens_log = self.file_settings['optimization']['stable_baselines3_options']['tensorboard_log']
        if tens_log is False:
            tens_log = None
        model_save = self.file_settings['optimization']['stable_baselines3_options']['model_save']
        gam = self.file_settings['optimization']['stable_baselines3_options']['gamma']
        sb_algo = self.file_settings['optimization']['stable_baselines3_options']['algorithm']
        if sb_algo == 'SAC':
            policy_kwargs = dict(net_arch=dict(pi=net1, qf=net2)) # pi is the actor network, while qf is the critic network
            model=SAC('MlpPolicy', env, verbose=1, 
            learning_rate=self.file_settings['optimization']['stable_baselines3_options']['learning_rate'], 
            buffer_size=self.file_settings['optimization']['stable_baselines3_options']['buffer_size'], 
            learning_starts=self.file_settings['optimization']['stable_baselines3_options']['learning_starts'], 
            batch_size=self.file_settings['optimization']['stable_baselines3_options']['batch_size'], 
            tau=self.file_settings['optimization']['stable_baselines3_options']['tau'], 
            gamma=self.file_settings['optimization']['stable_baselines3_options']['gamma'], 
            train_freq=self.file_settings['optimization']['stable_baselines3_options']['train_freq'], 
            gradient_steps=self.file_settings['optimization']['stable_baselines3_options']['gradient_steps'],
            tensorboard_log=tens_log,policy_kwargs=policy_kwargs)
        elif sb_algo == 'PPO':
            policy_kwargs = dict(net_arch=[dict(pi=net1, vf=net2)])
            model=PPO('MlpPolicy', env, verbose=1, 
            learning_rate=self.file_settings['optimization']['stable_baselines3_options']['learning_rate'], 
            n_steps=self.file_settings['optimization']['stable_baselines3_options']['n_steps'], 
            batch_size=self.file_settings['optimization']['stable_baselines3_options']['batch_size'], 
            n_epochs=self.file_settings['optimization']['stable_baselines3_options']['n_epochs'], 
            gamma=self.file_settings['optimization']['stable_baselines3_options']['gamma'], 
            gae_lambda=self.file_settings['optimization']['stable_baselines3_options']['gae_lambda'], 
            clip_range=self.file_settings['optimization']['stable_baselines3_options']['clip_range'],
            tensorboard_log=tens_log,policy_kwargs=policy_kwargs)
        elif sb_algo == 'A2C':
            policy_kwargs = dict(net_arch=[dict(pi=net1, vf=net2)])
            model=A2C('MlpPolicy', env, verbose=1, 
            learning_rate=self.file_settings['optimization']['stable_baselines3_options']['learning_rate'], 
            n_steps=self.file_settings['optimization']['stable_baselines3_options']['n_steps'], 
            gamma=self.file_settings['optimization']['stable_baselines3_options']['gamma'], 
            gae_lambda=self.file_settings['optimization']['stable_baselines3_options']['gae_lambda'], 
            ent_coef=self.file_settings['optimization']['stable_baselines3_options']['ent_coef'], 
            vf_coef=self.file_settings['optimization']['stable_baselines3_options']['vf_coef'],
            tensorboard_log=tens_log,policy_kwargs=policy_kwargs)
        elif sb_algo == 'DQN':
            policy_kwargs = dict(net_arch=net2)
            model = DQN("MlpPolicy", env, verbose=1,
                    train_freq=self.file_settings['optimization']['stable_baselines3_options']['train_freq'],
                    gradient_steps=self.file_settings['optimization']['stable_baselines3_options']['gradient_steps'],
                    gamma=self.file_settings['optimization']['stable_baselines3_options']['gamma'],
                    exploration_fraction=self.file_settings['optimization']['stable_baselines3_options']['exploration_fraction'],
                    exploration_final_eps=self.file_settings['optimization']['stable_baselines3_options']['exploration_final_eps'],
                    target_update_interval=self.file_settings['optimization']['stable_baselines3_options']['target_update_interval'],
                    learning_starts=self.file_settings['optimization']['stable_baselines3_options']['learning_starts'],
                    buffer_size=self.file_settings['optimization']['stable_baselines3_options']['buffer_size'],
                    batch_size=self.file_settings['optimization']['stable_baselines3_options']['batch_size'],
                    learning_rate=self.file_settings['optimization']['stable_baselines3_options']['learning_rate'],
                    tensorboard_log=tens_log,policy_kwargs=policy_kwargs)

        games_numbers = self.generation.total
        steps_per_game = len(self.file_settings['optimization']['order'])
        model.learn(total_timesteps=steps_per_game*games_numbers)
        model.save(model_save)
        obs = env.reset()


        track_file = open('optimization_track_file.txt','a')
        track_file.write("End of Optimization \n")
        track_file.close()

        opt.plotter()