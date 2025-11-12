import gymnasium as gym
from stable_baselines3 import SAC, PPO, A2C, DQN
import os
import matplotlib.pyplot as plt
import copy
from stable_baselines3.common.vec_env import DummyVecEnv
from stable_baselines3.common.monitor import Monitor
from midas.utils.metrics import Optimization_Metric_Toolbox

# env = gym.make("CartPole-v1", render_mode="rgb_array")

# model = A2C("MlpPolicy", env, verbose=1)
# model.learn(total_timesteps=10_000)

# vec_env = model.get_env()
# obs = vec_env.reset()
# for i in range(1000):
#     action, _state = model.predict(obs, deterministic=True)
#     obs, reward, done, info = vec_env.step(action)
#     vec_env.render("human")
#     # VecEnv resets automatically
#     # if done:
#     #   obs = vec_env.reset()

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