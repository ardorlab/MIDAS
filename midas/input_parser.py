## Import Block ##
import yaml
import logging
"""
Classes for parsing and cleansing input data from the user-specified '.yaml' MIDAS input file.

Created by Nicholas Rollins. 09/11/2024
"""


## Classes ##
class Input_Parser():
    """
    Centralized class for parsing user-supplied input arguments from the 
    MIDAS '.yaml' input file
    
    Written by Nicholas Rollins. 09/11/2024
    """
    def __init__(self, num_procs, inp_file):
        self.num_procs = int(num_procs)
        with open(inp_file) as f:
            self.file_settings = yaml.safe_load(f)
        
        self.parse_input_data()
        self.validate_input()
        
    def parse_input_data(self):
        """
        Interpret parsed input data.
        
        Written by Nicholas Rollins. 09/11/2024
        """
        ## Optimization Block ##
        info = self.file_settings['optimization']
        
        self.methodology = info['methodology'].lower()
        self.calculation_type = info['data_type'].lower() #!or external_model?
        self.results_dir_name = info['results_directory_name']
        self.population_size = info['population_size']
        self.num_generations = info['number_of_generations']
        self.parameters = info['objectives']
        
        self.genome = self.file_settings['decision_variables']['parameters'] #!TODO: add genome validation.
        
        ## Algorithm Block ##
        info = self.file_settings['algorithm']
        
        self.fitness = info['fitness'].lower()
        self.selection = info['selection'].lower()
        self.reproducer = info['reproducer'].lower()
        self.mutation_type = info['mutation_type'].lower()
        self.mutation_rate = float(info['mutation_rate'])
        
        ## Calculation Block ##
        info = self.file_settings['parcs_data']
        
        self.xs_lib = info['xs_library_path']
        self.power = float(info['power'])
        self.flow = float(info['flow'])
        self.inlet_temp = float(info['inlet_temperature'])
        self.num_assemblies = float(info['number_assemblies'])
        
        #!TODO: add core_dict and other values for PARCS/Simulate3
        
        return
    
    def validate_input(self):
        """
        Verify parsed input data. Input is sanitized to ensure proper entries and 
        formatting before being provided to the Optimizer.
        
        Written by Nicholas Rollins. 10/03/2024
        """
        if self.methodology not in ["genetic_algorithm"]:
            raise ValueError("Requested methodology '" + self.methodology + "' invalid.")
        
        try:
            if self.population_size.lower() in ["calculate_from_genes"]:
                self.population_size = self.population_size.lower()
            else:
                raise ValueError("Population size may be a positive number, or 'calculate_from_genes'.")
        except AttributeError:
                self.population_size = int(self.population_size)
        
        try:
            if self.num_generations.lower() in ["calculate_from_genes"]:
                self.num_generations = self.num_generations.lower()
            else:
                raise ValueError("Number of generations may be a positive number, or 'calculate_from_genes'.")
        except AttributeError:
                self.num_generations = int(self.num_generations)
        
        return