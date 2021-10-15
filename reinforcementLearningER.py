import os
import os.path
from os import path
import pathlib
import sys
import copy
import math
import numpy as np
import random 
from solution_types import evaluate_function
from metrics import Reinforcement_Learning_Metric_Toolbox
try:
   import torch
   import torch.nn as nn
   import torch.nn.functional as F
   from torch.autograd import Variable
   import torch.optim as optim
except:
   pass
else:
   pass

"""
This file contains a few of the classes and methods involved with 
performing an optimization via reinforcement learning.

Written by Gregory Delipei based on Alex Chrysler code. 9/30/2020
"""

class RL_dqn(object):
    """
    Object for reinforcement learning
    """

    def __init__(self, chromosomes,id_map,epsilon,return_to_state,hidden_units,learning_rate,discount):
        self.state = None
        self.reward = None
        self.state_genome_dict={}
        self.id_map = id_map
        self.map_genome2state(chromosomes)
        self.repeat_state = return_to_state
        self.neural_netPN = FeedForwardNeuralNet(hidden_units)
        self.neural_netTN = copy.deepcopy(self.neural_netPN)
        self.loss_fn = nn.MSELoss()
        self.learning_rate = learning_rate
        self.discount = discount
        self.epsilon = epsilon
        self.Q = None
        self.ind = None
        self.current_q = None
        self.optimizer = optim.SGD(self.neural_netPN.parameters(), lr=learning_rate) # learning rate in optimizer is different than RL learning rate

    def get_genome_from_state(self,state):
        """
        Returns the genome from the state
        """
        genome=[]
        state_id=0
        for i in self.id_map:
           if i ==1:
              for gene in self.state_genome_dict.keys():
                 if self.state_genome_dict[gene] == state[state_id]:
                   genome.append(gene)
              state_id +=1
           elif i==0:
              genome.append("Reflector")
        return(genome)


    def map_genome2state(self, chromosomes):
        """
        Converts the state values to the qtable indices
            - state is an array of fuel types configuration
            - the state is returned as a tuple
        """

        cid=0
        for chromosome in chromosomes:
           if 'map_state' in chromosomes[chromosome]:
              if chromosomes[chromosome]['map_state']:
                self.state_genome_dict[chromosome]=cid
                cid +=1

    def state_converter(self, csolution):
        """
        Converts the state values to the qtable indices
            - state is an array of fuel types configuration
            - the state is returned as a tuple
        """

        genome=csolution.genome
        genome_state=[]
        for gene in genome:
           if gene in self.state_genome_dict.keys():
              genome_state.append(self.state_genome_dict[gene])
           else:
              pass
        return(genome_state)


    def impose_state(self, csolution):
        """
        Imposes a desired state and updates the current reward
            - state is an array of fuel types configuration
        """
        self.state = self.state_converter(csolution)
        self.current_reward = csolution.fitness




    def map_net_output(self, net_output):
        """
        Neural Network Function:
        Maps the neural network output to the action with the highest Q value
        """
        _, ind = torch.max(net_output, -1)
        # possible range of ind is now 1 to max instead of 0 to max-1
        ind = torch.Tensor.item(ind) + 1
        print("Index from Q values: {}".format(ind))
        max_ind = len(net_output)
        num_assemblies = len(self.state)
        num_fuel_types = len(self.state_genome_dict.keys())
        new_act = [0, 0]
        fuel_types  = [x for x in self.state_genome_dict.values()]
        current_state = self.state
        print("Current State: {}".format(current_state))
        for i in range(num_assemblies):
            if ind > (i*max_ind)/num_assemblies and ind <= ((i+1)*max_ind)/num_assemblies:
                ass_to_change = i  # assembly 0 1 2 or 3 for toy ex
                print("Assembly to change: {}".format(ass_to_change))
                # reduces ind to one of [0 1 2 3 4 5] for toy ex
                ind = int(ind - (((i * max_ind) / num_assemblies) + 1))
                break
        if self.repeat_state == False:
            # removes the number of the fuel type currently present
            print("Fuel types: {}".format(fuel_types))
            print("Fuel type to remove: {}".format(
                current_state[ass_to_change]))
            fuel_types.remove(current_state[ass_to_change])
            print("Remaining fuel types: {}".format(fuel_types))
        new_act[0] = fuel_types[ind] # place 0 1 2 3 or 4 for fuel type
        new_act[1] = ass_to_change  # into assembly 0 1 2 or 3
        print("Action: {}".format(new_act))
        return new_act

    def random_act(self, net_output):
        """
        It randomly generates a state that will be passed through predict state 
        and update state such that if a random number is less than the Epsilon 
        Greedy value, that shall be taken as the next state and evaluated. 
        """
        print("Epsilon-Greedy Action")
        new_act = [0, 0]
        fuel_types  = [x for x in self.state_genome_dict.values()]
        num_assemblies = len(self.state)
        print(num_assemblies)
        current_state = copy.deepcopy(self.state)
        current_state = np.array(current_state)
        new_act[0] = fuel_types[random.randint(0,4)]
        new_act[1] = random.randint(0,num_assemblies-1)
        print("Action: {}" .format(new_act))
        #return the randomized state array to be fed into predict_state if Epsilon Greedy condition is chosen
        return new_act


    def predict_state(self,current_solution):
        """
        It performs the RL policy by updating the state based on an action
            - lr is the LEARNING_RATE
            - dc is the DISCOUNT

        First, the action is selected by taking the maximum of the qvalue
        for the current state.

        Second, the state is updated and the reward is calculated

        Third, the current qvalue before the action and the maximum future
        qvalue after the action are computed and combined with the reward to
        calculate the new qvalue for the initial state and action.
        """
        self.impose_state(current_solution)
        current_state = copy.deepcopy(self.state)
        #run a random number for the Epsilon Greedy condtition
        Erand = random.random()
        current_state = np.array(current_state)
        Q = self.neural_netPN(Variable(torch.from_numpy(current_state).type(torch.FloatTensor)))
        if Erand > self.epsilon:
            current_act = self.map_net_output(Q)
        else:
            current_act = self.random_act(Q)
        self.epsilon *= self.discount
        # Get the q values for each action in the current state
        current_q, ind = torch.max(Q, -1)
        current_q = torch.Tensor.item(current_q)
        new_state = copy.deepcopy(self.state)
        new_state[current_act[1]] = current_act[0]
        new_genome = self.get_genome_from_state(new_state)
        self.current_q = current_q
        self.Q = Q
        self.ind = ind
        return(new_genome)

    def update_state(self,new_solution):
        """
        It performs the RL policy by updating the state based on an action
            - lr is the LEARNING_RATE
            - dc is the DISCOUNT

        First, the action is selected by taking the maximum of the qvalue
        for the current state.

        Second, the state is updated and the reward is calculated

        Third, the current qvalue before the action and the maximum future
        qvalue after the action are computed and combined with the reward to
        calculate the new qvalue for the initial state and action.
        """     
        self.impose_state(new_solution)
        new_state = np.array(self.state)
        Q1 = self.neural_netPN(Variable(torch.from_numpy(new_state).type(torch.FloatTensor)))
        max_future_q, _ = torch.max(Q1, -1)
        max_future_q = torch.Tensor.item(max_future_q)
        qvalue_update = (1 - self.learning_rate) * self.current_q + self.learning_rate * (self.current_reward + self.discount * max_future_q)
        o_stdout = sys.stdout
        with open("rewards.txt","a") as f:
                sys.stdout = f
                print(str(self.current_reward))
        sys.stdout = o_stdout        
        # Creating Q_target for training neural network
        Q_target = self.Q.clone()
        Q_target = Variable(Q_target.data)
        Q_target[self.ind] = qvalue_update
        # Calculate loss
        loss = self.loss_fn(self.Q, Q_target)
        # Update neural network
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()


class FeedForwardNeuralNet(nn.Module):
    def __init__(self, hidden_units):
        super(FeedForwardNeuralNet, self).__init__()
        #Define layers of the neural network as l0, l1, ..., lx where x is the total number of layers minus 1

        num_layers = len(hidden_units)-1
        self.hidden = []
        for i in range(num_layers):
            self.hidden.append(nn.Linear(hidden_units[i], hidden_units[i+1], bias=False))
            self.add_module("l"+str(i), self.hidden[-1])

    def forward(self, x):
        for j in range(len(self.hidden)-1):
            x = F.relu(self.hidden[j](x))

        x = self.hidden[-1](x)
        print(x)
        return x

class ReinforcementLearning(object):
    def __init__(self,solution,
                      population,
                      generation,
                      parameters,
                      fitness,
                      file_settings
                ):

        self.solution = solution
        self.stored_solutions = []
        self.stored_fitness = []
        self.best_solution = None
        self.pouplation_size= population.size
        #self.episodes = parameters["episodes"]
        self.Eend = parameters["Eend"]
        self.epochs = generation
        #self.Nanneal = parameters["Nanneal"]
        #self.C = parameters["C"]
        self.type = parameters["name"]
        self.architecture = parameters["architecture"]
        self.epsilon= parameters["egreedy"]
        self.learning_rate = parameters["learning_rate"]
        self.discount = parameters["discount"]
        self.return_to_state = parameters["return_to_state"]
        self.fitness = fitness
        self.file_settings = file_settings

           
    def main_in_serial(self ):#I believe instead of self we may need things like TN and PN but I'm not sure
         """
         Performs optimization using reinforcement learning for a population size of one
         """ 
         opt = Reinforcement_Learning_Metric_Toolbox()

         track_file = open('optimization_track_file.txt','w')
         track_file.write("Beginning Optimization \n")
         track_file.close()

         # defining the active solution 
         active = self.solution()
         active.name = "initial_solution"
         active.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
         active.add_additional_information(self.file_settings)
         if active.fixed_genome:
            active.generate_initial_fixed(self.file_settings['genome']['chromosomes'],
                                          self.file_settings['optimization']['fixed_groups'])
         else:
            active.generate_initial(self.file_settings['genome']['chromosomes'])
         #print(active.genome)
         active.evaluate()
         one = []
         one.append(active)
         one = self.fitness.calculate(one)
         self.stored_solutions.append(one[0].genome)
         self.best_solution=one[0]
         self.stored_fitness.append(one[0].fitness)
         rl = RL_dqn(self.file_settings["genome"]["chromosomes"],
                     self.file_settings["genome"]["chromosomes"]["Assembly_One"]["map"], self.epsilon,
                     self.return_to_state,self.architecture,self.learning_rate,self.discount)
         #B = np.empty()
         for epoch in range(self.epochs):
                previous_solution = one[0]
                for iteration in range(self.pouplation_size):
                    o_stdout = sys.stdout
                    if path.exists("rewards.txt"):
                        with open("rewards.txt","a") as f:
                                sys.stdout = f
                                f.write(str(epoch))
                    else:
                        with open ("rewards.txt", "w") as f:
                                sys.stdout = f
                                f.write(str(epoch))    
                    sys.stdout = o_stdout
               #    import pdb;pdb.set_trace()
                    new_genome = rl.predict_state(previous_solution)
                    nsolution = self.solution()
                    nsolution.genome = new_genome
                    nsolution.name = "solution_{}_{}".format(epoch,iteration)
                    nsolution.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
                    nsolution.add_additional_information(self.file_settings)
                    nsolution.evaluate()
                    asolution = [nsolution]
                    asolution = self.fitness.calculate(asolution)
                    new_solution = asolution[0]
                    rl.update_state(new_solution)
                    #determining which solution makes the next generation 
                    self.stored_solutions.append(new_solution.genome)
                    self.stored_fitness.append(new_solution.fitness)
                    #temp = np.array(self.)
                    #B = np.concatenate((B, a2, ...))
                    if(new_solution.fitness>self.best_solution.fitness):
                        self.best_solution=new_solution
                    #store stuff into B class in here I believe    
                    #if mod(iteration,F_train)==0:
                        #sample memories here
                        #if epochs+1==self.epochs:
                            #y = B
                        #else:
                           # y = B + TN#.calcQ(B.State_next,B.action_next)
                        #I think there might be a GD in here already, but not sure
                    #if iteration % self.C==0:
                     #   self.neural_netTN = copy.deepcopy(self.neural_netPN)
                    #if iteration>self.Nanneal:
                     #   self.epsilon=self.Eend
                        

                opt.record_best_and_new_solution(self.best_solution,new_solution,rl)
                previous_solution = new_solution
               

         track_file = open('optimization_track_file.txt','a')
         track_file.write("End of Optimization \n")
         track_file.close()




