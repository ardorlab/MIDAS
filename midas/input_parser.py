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
    
    def yaml_linereader(data,keyword,default):
        """
        #!TODO: write docstring.
        
        Written by Nicholas Rollins. 10/03/2024
        """
        try:
            parsed_val = data[keyword]
            #!TODO: validatation of value based on keyword
        except:
            parsed_val = default
        return parsed_val
        
    def parse_input_data(self):
        """
        Interpret parsed input data.
        
        Written by Nicholas Rollins. 09/11/2024
        """
        ## Optimization Block ##
        info = self.file_settings['optimization']
        
        self.methodology = yaml_line_reader(info, 'methodology', None)
        self.calculation_type = yaml_line_reader(info, 'data_type', None) #!or external_model?
        self.results_dir_name = yaml_line_reader(info, 'results_directory_name', None)
        self.population_size = yaml_line_reader(info, 'population_size', 1)
        self.num_generations = yaml_line_reader(info, 'number_of_generations', 1)
        self.parameters = yaml_line_reader(info, 'objectives', None)
        
        self.genome = yaml_line_reader(self.file_settings['decision_variables'], 'parameters', None) #!TODO: add genome validation.
        
        ## Algorithm Block ##
        info = self.file_settings['algorithm']
        
        self.fitness = yaml_line_reader(info, 'fitness', None)
        self.selection = yaml_line_reader(info, 'selection', None)
        self.reproducer = yaml_line_reader(info, 'reproducer', None)
        self.mutation_type = yaml_line_reader(info, 'mutation_type', None)
        self.mutation_rate = yaml_line_reader(info, 'mutation_rate', None)
        
        ## Calculation Block ##
        info = self.file_settings['parcs_data']
        
        self.xs_lib = yaml_line_reader(info, 'xs_library_path', None)
        self.power = yaml_line_reader(info, 'power', 100)
        self.flow = yaml_line_reader(info, 'flow', None)
        self.inlet_temp = yaml_line_reader(info, 'inlet_temperature', None)
        self.num_assemblies = yaml_line_reader(info, 'number_assemblies', None)
        
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