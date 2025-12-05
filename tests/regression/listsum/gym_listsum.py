"""
Optimization problem: list summation 

Given a list size N and possible values for its elements in the range [0,9]
the objective is to find the list that maximizes the sum of its elements, which is [9,9,9,...,9]

There are 10^N possible solutions. 


From simple studies it seems we report the minimum, median, and maximum fitness values 
for different M number of random trials of a problem size of N=10:
   M = 1: 22/45/65
   M = 10: 47/60/81
   M = 100: 59/67/79
   M = 1000: 68/73/79
   M = 10000: 72/77/81
   M = 100000: 78/81/86
"""

import gymnasium as gym
from gymnasium import spaces
from stable_baselines3 import SAC, PPO, A2C, DQN
import os
import matplotlib.pyplot as plt
import copy
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.monitor import Monitor
from typing import Any, Generic, Optional, TypeVar, Union

import numpy as np
import logging
import os

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, input):
    """
    This function will compute the sum of the list elements write them into a dictionary in the solution.parameters object

    Written by Gregory Delipei. 09/10/2025
    """
    
    #Create results directory if it doesn't exist
    results_path = os.path.join(input.results_dir_name, solution.name)
    os.makedirs(results_path, exist_ok=True)
    
    #Create separate dictionary with parameters
    new_dict = {}
    new_dict["list_sum"] = list_summation(np.array(solution.chromosome, dtype=int))
   
    #Only give the optimizer the parameters which were included in the input file
    for key in new_dict:
        if key in input.objectives:
            solution.parameters[key]["value"] = new_dict[key]
        else:
            logger.info(f'Parameter {key} is not used')
    return solution

def list_summation(list_solution):
    return np.sum(list_solution)

def random_best_solution(M):
    current_best_solution = None
    current_best_fitness = -np.inf
    for _ in range(M):
        rnd_solution = np.random.randint(0, 10, N)
        rnd_fitness = list_summation(rnd_solution)
        if rnd_fitness > current_best_fitness:
            current_best_solution = rnd_solution
            current_best_fitness = rnd_fitness
    return current_best_solution, current_best_fitness

class Gym_Listsum_Env(gym.Env):
    """
    Class for wrapper adapted for Gym environment
    """
    metadata = {'render.modes': ['console']}
    # Define constants for clearer code

    def __init__(self,N):
        self.list_size = N
        self.action_lower=-1
        self.action_upper=1
        self.start = [1]*self.list_size
        self.current = copy.deepcopy(self.start)
        
        self.counter=0
        self.total_run=0
        self.cmap={}
        cmap_range=np.arange(0,self.list_size)
        for i in range(self.list_size):
            self.cmap[str(i)]=cmap_range[i]
        
        self.action_space = spaces.Discrete(self.list_size)
        
        self.observation_space = spaces.Box(low=-1.0, high=1.0, shape=(2,len(self.start)), dtype=np.float32)
        self.reset()
        
    def get_mapstate(self,solution):
        """
        Gets the current state in a normalized format.
        """

        mstate=np.zeros(self.list_size,dtype=np.int8)

        for i in range(len(solution)):
            mstate[i]=self.cmap[str(solution[i])]
        return(mstate)

    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        """
        Important: the observation must be a numpy array
        :return: (np.array) 
        """
        if seed is not None:
            super().reset(seed=seed)
        self.counter=0
        self.current = copy.deepcopy(self.start)

        cstate = self.get_mapstate(self.start)
        loc_state = np.zeros(cstate.shape[0])
        loc_state[self.counter]=1
        rstate=np.r_[[cstate],[loc_state]]
        return(rstate)

    def step(self, action):
        loc = self.counter
        self.current[loc] = int(action)

        if self.counter==len(self.start)-1:
            mid = 0
            self.total_run +=1
            reward = list_summation(np.array(self.current, dtype=int))
            info = {'list_sum':reward,
                    'State':self.current}
            if self.total_run ==1:
                self.best_solution=copy.deepcopy(self.current)
                self.best_reward = reward
            else:
                if self.best_reward<reward:
                    self.best_solution = copy.deepcopy(self.current)
                    self.best_reward = reward
        else:
            mid = self.counter + 1
            reward=0
            info={}
        self.counter+=1
        terminated = bool(self.counter==len(self.start))
        truncated = False
        
        cstate = self.get_mapstate(self.current)
        loc_state = np.zeros(cstate.shape[0])
        loc_state[mid]=1
        rstate=np.r_[[cstate],[loc_state]]
        return((rstate, reward, terminated, truncated, info))
    
    def render(self, mode='console'):
        if mode != 'console':
            mode='console'

    def close(self):
        pass


if __name__ == "__main__":
    N=10
    population_size=1 
    number_of_generations=1000
    policy_net= [16, 32, 16]
    qvalue_net= [16, 32, 16]
    train_freq= 10
    gradient_steps= 1
    exploration_fraction= 0.5
    exploration_final_eps= 0.02
    target_update_interval= 1
    learning_starts= 1
    buffer_size= 10
    batch_size= 32
    learning_rate= 0.0001
    gamma= 0.99
    tensorboard_log= False
    model_save= 'listsum_rl'
    action_space= 'discrete'
    observation_space= 'continuous'

    log_dir = 'rl_results'
    os.makedirs(log_dir, exist_ok=True)

    
    info_kwd=['list_sum','State']
    info_kwd = tuple(info_kwd)
    
    env = Gym_Listsum_Env(N)
    #env = gym.make("CartPole-v1", render_mode="rgb_array")
    env = Monitor(env,log_dir,info_keywords=info_kwd)
    env = DummyVecEnv([lambda: env])
    net1 = policy_net
    net2 = qvalue_net
    tens_log = tensorboard_log
    if tens_log is False:
        tens_log = None
    
    
    policy_kwargs = dict(net_arch=net2)
    model = DQN("MlpPolicy", env, verbose=1,
            train_freq=train_freq,
            gradient_steps=gradient_steps,
            gamma=gamma,
            exploration_fraction=exploration_fraction,
            exploration_final_eps=exploration_final_eps,
            target_update_interval=target_update_interval,
            learning_starts=learning_starts,
            buffer_size=buffer_size,
            batch_size=batch_size,
            learning_rate=learning_rate,
            tensorboard_log=tens_log,policy_kwargs=policy_kwargs)

    games_numbers = number_of_generations
    steps_per_game = N
    model.learn(total_timesteps=steps_per_game*games_numbers)
    model.save(model_save)
    obs = env.reset()

    vec_env = model.get_env()
    obs = vec_env.reset()
    for i in range(N):
        action, _state = model.predict(obs, deterministic=True)
        if i == N-1:
            last_solution = obs[0,0,:] 
            last_solution[-1]=int(action)
        obs, reward, done, info = vec_env.step(action)
        vec_env.render("human")
        # VecEnv resets automatically
        # if done:
        #   obs = vec_env.reset()

    print('The obtained reward is {:.2f}'.format(reward[0]))
    print('The obtained solution is: {}'.format(last_solution))
    



