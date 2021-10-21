import os
import sys
import copy
import math
import numpy
import random
from multiprocessing import Pool
from solution_types import evaluate_function
from metrics import Simulated_Annealing_Metric_Toolbox

"""
This file contains a few of the classes and methods involved with 
performing an optimization via simulated aneealing.
NOTE: This file does not include the means for reproduction, but is 
written as imf the class is included. 

Written by Johnny Klemes. 3/19/2020
"""

class Cooling_Schedule(object):
    """
    Class for Simulated Annealing cooling schedules. 
    
    THe cooling schedule sets the tolerance for accepting new solutions. 
    A high initial temperature indicates accepting a more new designs even if 
    they have a less favorable objective function. Thus the logarithmic cooling 
    schedule is favorable for this problem.
    
    There are two cooling schedules defined below. In both the temperature is 
    determined by the current era of a lifetime, represented by a piecewise 
    function. THe second cooling schedule is identical to the first but with
    two cycles.
   
    """
    def __init__(self,generation):
        self.generation = generation

    @staticmethod
    def piecewise(gen,lifet):
       if gen <= 3/12*lifet:
          temp = gen
          y = 2*numpy.log(temp/(4/12*lifet))
       elif gen <= 4/12*lifet and gen > 3/12*lifet:
          y = 1000
       elif gen > 4/12*lifet and gen <= 6/12*lifet:
          temp = gen - 4/12*lifet
          y = numpy.log(0.8*temp/(3/12*lifet))
       elif gen > 6/12*lifet and gen <= 7/12*lifet:
          y = 1000
       elif gen > 7/12*lifet and gen <= 10/12*lifet:
          temp = gen - 7/12*lifet
          y = 3*numpy.log(0.8*temp/(4/12*lifet))
       else:
          y = 1000

       return y

    @staticmethod
    def twicewise(gen,lifet):
       if gen <= 3/24*lifet:
          temp = gen
          y = 2*numpy.log(temp/(4/24*lifet))
       elif gen <= 4/24*lifet and gen > 3/24*lifet:
          y = 1000
       elif gen > 4/24*lifet and gen <= 6/24*lifet:
          temp = gen - 4/24*total
          y = numpy.log(0.8*temp/(3/24*lifet))
       elif gen > 6/24*lifet and gen <= 7/24*lifet:
          y = 1000
       elif gen > 7/24*lifet and gen <= 10/24*lifet:
          temp = gen - 7/24*lifet
          y = 3*numpy.log(0.8*temp/(4/24*lifet))
       elif gen > 10/24*lifet and gen <= 1/2*lifet:
          y = 1000
       elif gen > 1/2*lifet and gen <= 15/24*lifet:
          temp = gen - 1/2*lifet
          y = 2*numpy.log(temp/(4/24*lifet))
       elif gen <= 16/24*lifet and gen > 15/24*lifet:
          y = 1000
       elif gen <= 19/24*lifet and gen > 16/24*lifet:
          temp = gen - 16/24*lifet
          y = 3*numpy.log(0.8*temp/(4/24*lifet))
       else:
          y = 1000

       return y

class Exponential_Decreasing_Cooling_Schedule(object):
   """
   Logarithmic cooling schedule for simulated annealing. Implementing this because Johnny Klemes cooling schedules seem broken.
   His implementation that I pulled from Github is broken at the least, and there isn't sufficient documentation available to 
   understand how to fix it. 

   This cooling schedule is about as simple as it can get. 

   T = T0*alpha

   Where 0.9 < alpha < 1.0 and 1 < T0 < 10
   """
   def __init__(self,alpha,temperature):
      self.alpha = 1. - alpha
      self.temperature = temperature

   def update(self):
      """
      Updates the temperature of the cooling schedule.
      """
      if self.temperature <= 0.0001:
         self.temperature = 0.0001
      else:
         self.temperature *= self.alpha
      
class SimulatedAnnealing(object):
   def __init__(self,solution,
                      population,
                      generation,
                      cooling_schedule,
                      fitness,
                      mutation,
                      file_settings
                ):

      self.solution = solution
      self.population = population
      self.generation = generation
      self.cooling_schedule = cooling_schedule
      self.fitness = fitness
      self.mutation = mutation
      self.file_settings = file_settings
      self.number_generations_post_cleanup = 300 #Arbitrarily chosen default.
      if 'cleanup' in file_settings['optimization']:
          if file_settings['optimization']['cleanup']['perform']:
              self.perform_cleanup = True
          else: 
              self.perform_cleanup = False
      else:
          self.perform_cleanup = False

   def main_in_serial(self):
         """
         Performs optimization using simulated annealing for a population size of one
         """
         opt = Simulated_Annealing_Metric_Toolbox()

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
         # will store the fit for active solution
         one = []
         one.append(active)

         one = self.fitness.calculate(one)
         for self.generation.current in range(self.generation.total):
            for number in range(self.population.size):
               challenge = self.solution()
               challenge.genome = self.mutation.reproduce(active.genome)
               challenge.name = "solution_{}_{}".format(self.generation.current,number)
               challenge.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
               challenge.add_additional_information(self.file_settings)
               challenge.evaluate()

               test = [challenge]
               test = self.fitness.calculate(test)
   
               #determining which solution makes the next generation 
               opt.record_best_and_new_solution(active,challenge,self.cooling_schedule)
               acceptance = numpy.exp(-1*(challenge.fitness-active.fitness)/self.cooling_schedule.temperature)
               if challenge.fitness < active.fitness:
                   active = copy.deepcopy(challenge)
               elif random.uniform(0,1) < acceptance:
                   active = challenge
            self.cooling_schedule.update()
   
         track_file = open('optimization_track_file.txt','a')
         track_file.write("End of Optimization \n")
         track_file.close()

         # plot the parameters over time
         opt.plotter()
         self.cleanup()

   def cleanup(self):
      """
      Deletes solution results that are no longer relevant to the optimization.

      Written by Brian Andersen 10/28/2020.
      """
      if self.perform_cleanup:
            if self.generation.current <= self.number_generations_post_cleanup:
                pass
            else:
               file_name = f'initial_solution'
               if os.path.isdir(file_name):
                   delete_file = True
                   for parent,child in zip(self.population.parents,self.population.children):
                       if parent.name == file_name:
                           delete_file = False
                           break
                       elif child.name == file_name:
                           delete_file = False
                           break
                   if delete_file:
                       os.rmdir(file_name)
               for clgen in range(self.generation.current):
                   for clpop in range(self.population.size):
                       file_name = f"solution_{clgen}_{clpop}"
                       if os.path.isdir(file_name):
                           delete_file = True
                           for parent,child in zip(self.population.parents,self.population.children):
                               if parent.name == file_name:
                                   delete_file = False
                                   break
                               elif child.name == file_name:
                                   delete_file = False
                                   break
                           if delete_file:
                               os.rmdir(file_name)   
