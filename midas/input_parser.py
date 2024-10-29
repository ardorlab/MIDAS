## Import Block ##
import yaml
import logging
from pathlib import Path
"""
Classes for parsing and cleansing input data from the user-specified '.yaml' MIDAS input file.

Created by Nicholas Rollins. 09/11/2024
"""


## Classes ##
def yaml_line_reader(data,keyword,default):
    """
    Parse the data of a given keyword from the '.yaml' input file data.
    If the keyword is not found, revert to a provided default value. Then,
    validate the parsed value to ensure it is formatted correctly and/or 
    is a supported option for that keyword.
    
    Written by Nicholas Rollins. 10/03/2024
    """
    try:
        parsed_val = data[keyword]
    except:
        parsed_val = default
    parsed_val = validate_input(keyword, parsed_val)
    return parsed_val

def validate_input(keyword, value):
    """
    Verify parsed input data. Input is sanitized to ensure proper entries and 
    formatting before being provided to the Optimizer.
    
    Written by Nicholas Rollins. 10/03/2024
    """
    ## Initialize logging for the present file
    logger = logging.getLogger("MIDAS_logger")
    
## Optimization Block ##
    if keyword == 'methodology':
        value = str(value).lower().replace(' ','_')
        if value not in ["genetic_algorithm","bayesian_optimization"]:
            raise ValueError("Requested methodology '" + self.methodology + "' invalid.")
    
    elif keyword == 'code_type':
        value = str(value).lower()
        if value not in ["parcs342"]:
            raise ValueError("Code types currently supported: PARCS342.")
    
    elif keyword == 'data_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["single_cycle","eq_cycle"]:
            raise ValueError("Data type not supported.")
    
    elif keyword == 'results_directory_name':
        value = Path(str(value).replace(' ','_'))
    
    elif keyword == 'population_size':
        try:
            value = value.lower().replace(' ','_')
            if value not in ["calculate_from_genes"]:
                raise ValueError("Population size may be a positive number, or 'calculate_from_genes'.")
        except AttributeError:
                value = int(value)
    
    elif keyword == 'number_of_generations':
        try:
            value = value.lower().replace(' ','_')
            if value not in ["calculate_from_genes"]:
                raise ValueError("Number of generations may be a positive number, or 'calculate_from_genes'.")
        except AttributeError:
                value = int(value)

    elif keyword == 'batch_size':
        value = int(value)
    
    elif keyword == 'solution_symmetry':
        value = str(value).lower().replace(' ','_')
        if value not in ['octant','quarter','full']:
            raise ValueError("Symmetry of the solution list must be octant, quarter, or full.")
    
    elif keyword == 'objectives':
        if isinstance(value, dict):
            new_dict = {}
            #check objectives/constraints
            for key, item in value.items():
                new_key = str(key).lower().replace(' ','_')
                if new_key not in ['max_boron',
                                   'pinpowerpeaking',
                                   'fdeltah',
                                   'cycle_length']:
                    raise ValueError(f"Requested objective/constraint '{key}' not supported.")
                new_item = {}
                if isinstance(item, dict):
                    #check goals, weights, and targets
                    for subkey, subitem in item.items():
                        new_subkey = str(subkey).lower()
                        if new_subkey == 'goal':
                            new_subitem = str(subitem).lower().replace(' ','_')
                            if new_subitem not in ['maximize',
                                                   'minimize',
                                                   'meet_target',
                                                   'greater_than_target',
                                                   'less_than_target']:
                                raise ValueError(f"Requested objective/constraint goal '{subitem}' not supported.")
                        elif new_subkey == 'weight':
                            new_subitem = float(subitem)
                            if not 0.0 < new_subitem:
                                raise ValueError(f"Requested weight for {key} must be a positive non-zero real number.")
                        elif new_subkey == 'target':
                            new_subitem = float(subitem)
                        new_item[new_subkey] = new_subitem #save modified parameter
                    #check parameters logic
                    if 'goal' not in new_item:
                        raise ValueError(f"'Goal' parameter missing for {key}.")
                    if 'weight' not in new_item:
                        raise ValueError(f"'Weight' parameter missing for {key}.")
                    if new_item['goal'] in ['maximize','minimize']:
                        if 'target' in new_item:
                            logger.warning(f"Target provided for {key} with requested goal '{subitem}'. Target will be ignored.")
                    else:
                        if 'target' not in new_item:
                            raise ValueError(f"'Target' parameter missing for {key}.")
                else:
                    raise ValueError("Requested objective/constraint must be nested with its applicable parameters.")
                new_dict[new_key] = new_item #save modified objective/constraint
            return new_dict
        else:
            raise ValueError("'Objectives' must be nested with objectives/constraints and their parameters.")
    
## Algorithm Block ##
    elif keyword == 'selection':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key == 'fitness':
                    new_item = str(item).lower().replace(' ','_')
                    if new_item not in ['weighted']:
                        raise ValueError(f"Requested fitness type '{item}' not supported.")
                elif new_key =='method':
                    new_item = str(item).lower()
                    if new_item not in ['tournament','roulette']:
                        raise ValueError(f"Requested selection method '{item}' not supported.")
                new_dict[new_key] = new_item
            return new_dict
        else:
            raise ValueError("'Selection' must be nested with its parameters.")
    
    elif keyword == 'reproducer':
        value = str(value).lower().replace(' ','_')
        if value not in ["standard"]:
            raise ValueError("Reproducer type not supported.")
    
    elif keyword == 'mutation_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["mutate_by_gene"]:
            raise ValueError("Mutation type not supported.")
    
    elif keyword == 'mutation_rate':
        new_value = str(value).replace(' ','').split(',')
        new_dict = {}
        if len(new_value) >= 2:
            new_dict['initial_rate'] = float(new_value[0])
            new_dict['final_rate'] = float(new_value[1])
            if len(new_value) > 2:
                logger.warning("Only two entries expected for 'mutation_rate' (initial_rate, final_rate). Other entries are ignored.")
        else:
            new_dict['initial_rate'] = float(new_value[0])
            new_dict['final_rate'] = float(new_value[0])
        return new_dict
    
    elif keyword == 'acquisition_function':
        value = str(value).lower()
        if value == 'expected improvement':
            value = str("EI")
        elif value == 'probability of improvement':
            value = str("PI")
        elif value == 'lower confidence bound':
            value = str("LCB")
        elif value == 'gp hedge':
            value = str("gp_hedge")
        if value not in ["EI","PI","LCB","gp_hedge"]:
            raise ValueError("Acquisition function not supported.")
## Fuel Assembly Block ##
    elif keyword == 'assembly_options':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                #check reflector options
                if new_key == 'reflectors':
                    new_item = {}
                    if isinstance(item, dict):
                        #check assembly type
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types and cross sections
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'refl_type':
                                        new_subsubitem = str(subsubitem).lower()
                                        if new_subsubitem not in ['all','radial','top','bot']:
                                            raise ValueError(f"Reflector type for '{subsubkey}' must be radial, top, bottom, or all.")
                                    elif new_subsubkey == 'serial':
                                        new_subsubitem = str(subsubitem)
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested reflector missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Reflectors option missing reflectors and their parameters.")
                elif new_key == 'blankets':
                    new_item = {}
                    if isinstance(item, dict):
                        #check assembly type
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check cross sections
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'serial':
                                        new_subsubitem = str(subsubitem)
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested blanket missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Blankets option missing blankets and their parameters.")
                elif new_key == 'fuel':
                    new_item = {}
                    if isinstance(item, dict):
                        #check assembly type
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types and cross sections
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'type':
                                        new_subsubitem = int(subsubitem)
                                    elif new_subsubkey == 'serial':
                                        new_subsubitem = str(subsubitem)
                                    elif new_subsubkey == 'blanket':
                                        new_subsubitem = str(subsubitem)
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested fuel type missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Fuels option missing fuels and their parameters.")
            #check parameters logic
            if 'fuel' not in new_dict or len(new_dict['fuel']) == 0:
                raise ValueError("Assembly options must include fuel types.")
            if 'reflectors' not in new_dict or len(new_dict['reflectors']) == 0:
                raise ValueError("Assembly options must include at least one reflector type under key 'reflectors'.")
            if 'blankets' in new_dict and len(new_dict['blankets']) == 0:
                raise ValueError("Blanket options must include at least one blanket type.")
            radial_exists = False
            for key, value in new_dict['reflectors'].items():
                if "refl_type" not in value and "serial" not in value:
                    raise ValueError(f"'refl_type' or 'serial' parameters missing from '{key}'.")
                elif value['refl_type'] == 'radial' or value['refl_type'] == 'all':
                    radial_exists = True
            if not radial_exists:
                raise ValueError("Assembly options must include at least one reflector labeled 'radial' or 'all'.")
            if 'blankets' in new_dict:
                for key, value in new_dict['blankets'].items():
                    if "serial" not in value:
                        raise ValueError(f"'serial' parameter missing from '{key}'.")
            list_unique_fuel_types = []
            for key, value in new_dict['fuel'].items():
                if "type" not in value and "serial" not in value:
                    raise ValueError(f"'type' or 'serial' parameters missing from '{key}'.")
                if 'blanket' in new_dict['fuel'][key]:
                    try:
                        if new_dict['fuel'][key]['blanket'] not in new_dict['blankets']:
                            raise ValueError(f"Requested blanket '{new_dict['fuel'][key]['blanket']}'" + \
                                              f" for fuel type '{key}' does not exist in blankets.")
                    except KeyError:
                        raise ValueError(f"Requested blanket '{new_dict['fuel'][key]['blanket']}'" + \
                                              f" for fuel type '{key}' does not exist in blankets.")
                if value['type'] not in list_unique_fuel_types:
                    list_unique_fuel_types.append(value['type'])
                else:
                    raise ValueError("All fuel types must have a unique 'type' value.")
            return new_dict
        else:
            raise ValueError("Assembly optinos must be nested with reflectors, fuels, and/or blankets with their parameters.")
    
## Genome Block ##
    elif keyword in ['parameters', 'batches']:
        new_dict = {}
        if isinstance(value, dict):
            for key, item in value.items():
                new_key = str(key)
                new_dict[new_key] = {}
                #check decision variable options
                if isinstance(value[key], dict):
                    for subkey, subitem in item.items():
                        new_subkey = str(subkey).lower()
                        if new_subkey == 'map':
                            new_dict[new_key][new_subkey] = subitem
                        elif new_subkey == 'constraint':
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types and values
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower()
                                    if new_subsubkey == 'type':
                                        new_subsubitem = str(subsubitem).lower().replace(' ','_')
                                        if new_subsubitem not in ['max_quantity','less_than_variable']:
                                            raise ValueError(f"Requested decision variable constraint type '{subsubitem}' not supported.")
                                    elif new_subsubkey == 'value':
                                        try:
                                            new_subsubitem = int(subsubitem)
                                        except ValueError:
                                            new_subsubitem = str(subsubitem)
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                if not subitem: #allow constraint option to be "None".
                                    new_subitem = None
                            new_dict[new_key][new_subkey] = new_subitem
                else:
                    raise ValueError(f"Decision variables '{key}' must be nested with its parameters.")
            #check decision variable logic
            for key, value in new_dict.items():
                if 'map' not in value:
                    raise ValueError(f"Decision variable '{key}' must include a 'map' parameter.")
                if 'constraint' in value:
                    if 'type' not in value['constraint']:
                        raise ValueError(f"Decision variable '{key}' includes a constraint but no 'type' parameter.")
                    elif 'value' not in value['constraint']:
                        raise ValueError(f"Decision variable '{key}' includes a constraint but no 'value' parameter.")
                    elif value['constraint']['type'] == 'max_quantity':
                        if not isinstance(value['constraint']['value'], int):
                            raise ValueError(f"Maximum quantity constraint for decision variable '{key}' must be an integer.")
                    elif value['constraint']['type'] == 'less_than_variable':
                        if not isinstance(value['constraint']['value'], str):
                            raise ValueError(f"'Less Than' constraint for decision variable '{key}' must be a valid decision variable name.")
                        elif value['constraint']['value'] not in new_dict.keys():
                            raise ValueError(f"'Less Than' constraint for decision variable '{key}' must be a valid decision variable name.")
                        elif value['constraint']['value'] == key:
                            raise ValueError(f"'Less Than' constraint for decision variable '{key}' may not be '{key}'.")
                else:
                    new_dict[key]['constraint'] = None
                if keyword == 'batches': #!TODO: add a check to make sure batch 0 exists, and batches are numbered in order.
                    new_key = key.lower().replace(' ','_').split('_')
                    try:
                        batch_num = int(new_key[-1])
                    except ValueError:
                        raise ValueError("Please restrict the 'batches' names to the form e.g. 'batch 0' (zero-indexed).")
            return new_dict
        else:
            if not value:
                return None
            else:
                raise ValueError(f"Decision variable '{keyword}' must be nested with parameter options and their parameters.")
    
## Calculation Block ##
    elif keyword == 'num_rows':
        value = int(value)
    
    elif keyword == 'num_cols':
        value = int(value)
    
    elif keyword == 'number_assemblies':
        value = int(value)
    
    elif keyword == 'core_symmetry':
        value = str(value).lower()
        if value not in ['full','quarter']:
            raise ValueError("Requested core symmetry (used for printing) not valid.")
    
    elif keyword == 'xs_library_path':
        value = Path(str(value))
    
    elif keyword == 'xs_extension':
        value = str(value).split('.')[-1] #this supports both e.g. ".exe" and "exe".
        if value: #skip if no extension
            value = "." + value
    
    elif keyword == 'power':
        value = float(value)
    
    elif keyword == 'flow':
        value = float(value)
    
    elif keyword == 'inlet_temperature':
        value = float(value)
    
    elif keyword == 'th_fdbk':
        value = bool(value)
    
    elif keyword == 'pin_power_recon':
        value = bool(value)
    
    elif keyword == 'num_axial_nodes':
        value = int(value)
    
    elif keyword == 'axial_nodes':
        value = [a.strip(' ') for a in value.split(',')]
        new_value = []
        for node in value:
            if "*" in node:
                new_value.append(str(node)) #shorthand for repeated floats (e.g. 15*25.739) are left as strings to be evaluated elsewhere.
            else:
                new_value.append(float(node))
        return new_value
    
    elif keyword == 'boc_core_exposure':
        try:
            value = float(value)
        except ValueError:
            raise ValueError("'boc_core_exposure' must be a real number.")
    
    return value


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
    
    def parse_input_data(self):
        """
        Interpret parsed input data.
        
        Written by Nicholas Rollins. 09/11/2024
        """
        #!TODO: add input validation, especially for inputs as dicts
        #!TODO: add defaults
    ## Optimization Block ##
        info = self.file_settings['optimization']
        
        self.methodology = yaml_line_reader(info, 'methodology', 'genetic_algorithm')
        self.code_interface = yaml_line_reader(info, 'code_type', 'PARCS342')
        self.calculation_type = yaml_line_reader(info, 'data_type', 'single_cycle')
        self.results_dir_name = yaml_line_reader(info, 'results_directory_name', 'output_files')
        self.population_size = yaml_line_reader(info, 'population_size', 1)
        self.num_generations = yaml_line_reader(info, 'number_of_generations', 1)
        self.batch_size = yaml_line_reader(info, 'batch_size', 1)
        self.symmetry = yaml_line_reader(info, 'solution_symmetry', 'octant')
        self.objectives = yaml_line_reader(info, 'objectives', None)
        
    ## Algorithm Block ##
        info = self.file_settings['algorithm']
        
        selection_default = {'fitness':'weighted','method':'roulette'}
        self.selection = yaml_line_reader(info, 'selection', selection_default)
        self.reproducer = yaml_line_reader(info, 'reproducer', 'standard')
        self.mutation_type = yaml_line_reader(info, 'mutation_type', 'mutate_by_gene')
        self.mutation_rate = yaml_line_reader(info, 'mutation_rate', 0.5)
        self.acquisition_function = yaml_line_reader(info, 'acquisition_function', 'gp_hedge')
        
    ## Fuel Assembly Block ##
        self.fa_options = yaml_line_reader(self.file_settings, 'assembly_options', None)
        
    ## Genome Block ##
        self.genome = yaml_line_reader(self.file_settings['decision_variables'], 'parameters', None)
        self.batches = yaml_line_reader(self.file_settings['decision_variables'], 'batches', None)
        #check that decision variable options are valid.
        for key, value in self.genome.items():
            if key not in self.fa_options['fuel']:
                raise ValueError(f"Decision variable option '{key}' not found in the list of fuel types under 'assembly_options'.")
        if self.calculation_type == 'eq_cycle' and not self.batches:
            raise ValueError(f"'Batches' must be specified in Decision Variables for the 'EQ Cycle' type.")
        
    ## Calculation Block ##
        info = self.file_settings['parcs_data'] #!TODO: this needs to be set by calculation_type
        
        self.nrow = yaml_line_reader(info['map'], 'num_rows', 17)
        self.ncol = yaml_line_reader(info['map'], 'num_cols', 17)
        self.num_assemblies = yaml_line_reader(info['map'], 'number_assemblies', 193)
        self.map_size = yaml_line_reader(info['map'], 'core_symmetry', 'full')
        self.xs_lib = yaml_line_reader(info, 'xs_library_path', '../../') #as a relative path, this assumes the needed cross sections are in the base directory for the job.
        self.xs_extension = yaml_line_reader(info, 'xs_extension', '')
        self.power = yaml_line_reader(info, 'power', 3800.0)
        self.flow = yaml_line_reader(info, 'flow', 18231.89)
        self.inlet_temp = yaml_line_reader(info, 'inlet_temperature', 565.0)
        self.th_fdbk = yaml_line_reader(info, 'th_fdbk', True)
        self.pin_power_recon = yaml_line_reader(info, 'pin_power_recon', True)
        self.number_axial = yaml_line_reader(info, 'num_axial_nodes', 19)
        self.axial_nodes = yaml_line_reader(info, 'axial_nodes', "16.12, 20.32, 15*25.739, 20.32, 16.12")
        self.boc_exposure = yaml_line_reader(info, 'boc_core_exposure', 0.0)
        
        return