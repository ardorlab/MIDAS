import os
import sys
import copy
import math
import numpy
import random

"""
This file contains a few of the classes and methods involved with 
performing an optimization via simulated aneealing.
NOTE: This file does not include the means for reproduction, but is 
written as imf the class is included. 

Written by Johnny Klemes. 3/19/2020
"""

class Exponential_Decreasing_Cooling_Schedule(object):
   """
   Exponentially decreasing cooling schedule for use with the simulated annealing
   optimization algorithm. 

   Written by Brain Andersen. 6/1/2020
   """
   def __init__(self,alpha,temperature,no_improvement_counter):
      self.alpha = 1. - alpha
      self.temperature = temperature
      self.no_improvement_limit = no_improvement_counter

   def update(self):
      """
      Updates the temperature of the cooling schedule.

      Written by Brian Andersen. 6/1/2020
      """
      if self.temperature <= 0.0001:
         self.temperature = 0.0001
      else:
         self.temperature *= self.alpha
      