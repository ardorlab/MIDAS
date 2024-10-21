import numpy as np 
from skopt import Optimizer
from skopt import gp_minimize
from skopt.space import Real, Integer, Categorical
from skopt.plots import plot_convergence
from multiprocessing import Pool
import matplotlib.pyplot as plt 

