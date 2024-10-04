## Import Block ##
import numpy as np


## Classes ##
class Randomizer():
    """
    Class that generates random solutions for the selected applications.
       - Fully random solutions from scratch.
       - Random perturbations for n decision variables based on a given.
       - Problem constraints are preserved. 
    
    Placeholder for the moment.
    Written by Gregory Delipei. 10/03/2024
    """
    def __init__(self, inp_lines):
        self.input = inp_lines
    
    def generate(self):
        """
        Function that generates a random solution from scratch.
        Preservation of problem constraints.
        """
        new_solution = np.random.randn()

        return(new_solution)
    
    def perturb(self,current_solution,n):
        """
        Function that perturbs an existing solution in 
        n different random locations.
        Preservation of problem constraints.
        """
        new_solution = current_solution + n*np.random.randn()

        return(new_solution)
    
    
 