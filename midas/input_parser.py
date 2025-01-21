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

## General Settings Block ##
    if keyword == 'debug_mode':
        value = bool(value)
    
    elif keyword == 'results_directory_name':
        value = Path(str(value).replace(' ','_'))
    
    elif keyword == 'set_seed':
        if value:
            value = int(value)
    
    elif keyword == 'clear_results':
        value = str(value).lower().replace(' ','_')
        #variations on "none" or "keep" are equivalent; "none" is the preferred variation.
        if value not in ["all", "all_but_best", "none", "keep", "keep_all"]:
            raise ValueError("Clear results request type is invalid or not supported.")
    
    elif keyword == 'optimizer':
        value = str(value).lower().replace(' ','_')
        if value not in ["genetic_algorithm","bayesian_optimization"]:
            raise ValueError("Requested methodology '" + value + "' invalid.")
    
    elif keyword == 'code_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["parcs342", "parcs343", "nuscale_database"]:
            raise ValueError("Code types currently supported: PARCS342, PARCS343, NuScale_Database.")
    
    elif keyword == 'calc_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["single_cycle","eq_cycle"]:
            raise ValueError("Data type not supported.")
    
## Optimization Block ##
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
    
    
    elif keyword == 'solution_symmetry':
        value = str(value).lower().replace(' ','_')
        if value not in ['octant','quarter','full']:
            raise ValueError("Symmetry of the solution list must be octant, quarter, or full.")
        
    elif keyword == 'termination_criteria':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key =='method':
                    new_item = str(item).lower()
                    if new_item not in ['consecutive','spearman','none']:
                        raise ValueError(f"Requested termination criteria method '{item}' not supported.")
                    if new_item == 'spearman':
                        logger.warning("Spearman rank requires burnup data for each assembly, currently only the NuScale look up table has this functionality")
                    if 'termination_generations' not in value.keys():
                        raise ValueError(f"'termination_generations' must be provided if termination criteria is requested in yaml file")
                    
                elif new_key =='termination_generations':
                    new_item = int(item)

                new_dict[new_key] = new_item
             
            return new_dict
    
    elif keyword == 'objectives':
        if isinstance(value, dict):
            new_dict = {}
            #check objectives/constraints
            for key, item in value.items():
                new_key = str(key).lower().replace(' ','_')
                if new_key not in ['max_boron',
                                   'pinpowerpeaking',
                                   'fdeltah',
                                   'cycle_length',
                                   'assembly_burnup',
                                   'cycle_cost',
                                   'av_fuelenrichment']:
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
                        elif new_subkey == 'target':
                            new_subitem = float(subitem)
                        elif new_subkey == 'settings':
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey = str(subsubkey).lower().replace(' ','_')
                                    new_subsubitem = str(subsubitem).lower().replace(' ','_')
                                    if new_subsubkey == 'scope':
                                        if new_subsubitem not in ['full_core','feed_batch_only'] and new_key == 'av_fuelenrichment':
                                            raise ValueError(f"Requested setting '{subitem}' not supported for objective '{key}'.")
                                new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError(f"Requested settings for objective '{key}' must be nested with its applicable parameters.")
                        new_item[new_subkey] = new_subitem #save modified parameter
                    #check parameters logic
                    if 'goal' not in new_item:
                        raise ValueError(f"'Goal' parameter missing for {key}.")
                    if 'weight' not in new_item:
                        new_item['weight'] = 1.0 #assume weight value
                    if new_item['goal'] in ['maximize','minimize']:
                        if 'target' in new_item:
                            logger.warning(f"Target provided for {key} with requested goal '{subitem}'. Target will be ignored.")
                    else:
                        if 'target' not in new_item:
                            raise ValueError(f"'Target' parameter missing for {key}.")  
                    if new_key == 'av_fuelenrichment':
                        if 'settings' in new_item:
                            if 'scope' not in new_item['settings']:
                                new_item['settings']['scope'] = 'full_core' #default value
                        else:
                            new_item['settings'] = {}
                            new_item['settings']['scope'] = 'full_core' #default value 
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
                    if new_item not in ['tournament','roulette','random','ktournament','truncation','sus']:
                        raise ValueError(f"Requested selection method '{item}' not supported.")         
                elif new_key == 'k':
                    new_item = int(item)
                    if new_item < 0:
                        raise ValueError("k parameter must be greater than 1 and less than population_size")             
                new_dict[new_key] = new_item

            if 'k' not in new_dict.keys() and new_dict['method'] == 'ktournament':
                new_key = 'k'
                new_item = 4
                new_dict[new_key] = new_item
                logger.warning("'k' parameter is missing from input while ktournament selection method is used, k is set to default value of 4.")  
                
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
    
    elif keyword == 'crossover':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key == 'method':
                    new_item = str(item).lower().replace(' ','_')
                    if new_item not in ['uniform', 'random_element', 'one_point', 'two_point']:
                        raise ValueError(f"Requested crossover method '{item}' not supported.")
                elif new_key =='crossover_rate':
                    new_item = float(item)
                    if item > 1:
                        raise ValueError("Crossover rate must be 1.0 or lower.")
                    
                elif new_key =='num_swaps':
                    new_item = int(item)
                    if item < 1:
                        raise ValueError("num_swaps must be 1 or higher and an integer.")

                new_dict[new_key] = new_item
            
            return new_dict
        else:
            raise ValueError("'crossover' must be nested with its parameters.")

    elif keyword == 'elites':
        value = float(value)
        if value > 1.0 and value.is_integer() == False:
            raise ValueError("'elites' value can either be between 0 and 1 to represent a percentage or an integer equal to or greater than 1 to represent a number of elites")

    elif keyword == 'acquisition_function':
        value = str(value).lower().replace(' ','_')
        if value in ['expected_improvement','ei']:
            value = 'EI'
        elif value in ['probability_of_improvement','pi']:
            value = 'PI'
        elif value in ['lower_confidence_bound','lcb']:
            value = 'LCB'
        elif value in ['upper_confidence_bound','ucb']:
            value = 'UCB'
        #Verification that valid acquisition function is specified
        if value not in ['EI','PI','LCB','UCB']:
            raise ValueError("Acquisition function not supported. Supported acquisition types are EI, PI, UCB, LCB.")
    
    elif keyword == 'exploration_exploitation_factor':
        value = float(value)
    
    elif keyword == 'kernel_smoothness_factor':
        value = float(value)
        
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
                                    elif new_subsubkey == 'enrichment':
                                        new_subsubitem = float(subsubitem)
                                        if new_subsubitem >= 1.0:
                                            new_subsubitem = new_subsubitem/100 #change weight percent to weight fraction
                                    elif new_subsubkey == 'hm_loading':
                                        new_subsubitem = float(subsubitem) # kg
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
                                    elif new_subsubkey == 'enrichment':
                                        new_subsubitem = float(subsubitem)
                                        if new_subsubitem >= 1.0:
                                            new_subsubitem = new_subsubitem/100 #change weight percent to weight fraction
                                    elif new_subsubkey == 'hm_loading':
                                        new_subsubitem = float(subsubitem) # kg
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
            if value:
                raise ValueError("Assembly options must be nested with reflectors, fuels, and/or blankets with their parameters.")
    
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
    elif keyword == 'exec_walltime':
        value = int(value)
        if value <= 0:
            raise ValueError("'exec_walltime' must be a positive number, measured in seconds.")
    
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
        value = [a.strip() for a in value.strip().replace(', ',',').replace(' ',',').split(',')]
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
    
    elif keyword=='depletion_steps':
        value = [a.strip() for a in value.strip().replace(', ',',').replace(' ',',').split(',')]
        new_value = []
        for step in value:
            if "*" in step:
                s_step = step.split('*')
                new_value.extend(int(s_step[0])*[float(s_step[1])])
            else:
                new_value.append(float(step))
        return new_value
    
    return value


class Input_Parser():
    """
    Centralized class for parsing user-supplied input arguments from the 
    MIDAS '.yaml' input file
    
    Written by Nicholas Rollins. 09/11/2024
    """
    def __init__(self, num_procs, inp_file):
        self.num_procs = int(num_procs)
        self.job_name = ".".join(inp_file.split('.')[:-1])
        with open(inp_file) as f:
            try:
                self.file_settings = yaml.safe_load(f)
            except yaml.parser.ParserError:
                raise yaml.parser.ParserError("Trouble reading the '.yaml' input file. Please check the integrity of the input, including the consistency of spaces!")
        
        self.parse_input_data()
    
    def parse_input_data(self):
        """
        Interpret parsed input data.
        
        Written by Nicholas Rollins. 09/11/2024
        """
    ## General Settings Block ##
        try:
            info = self.file_settings['general']
        except KeyError:
            info = None
        
        self.debug_mode = yaml_line_reader(info, 'debug_mode', False)
        self.results_dir_name = yaml_line_reader(info, 'results_directory_name', 'output_files')
        self.set_seed = yaml_line_reader(info, 'set_seed', None)
        self.clear_results = yaml_line_reader(info, 'clear_results', 'all_but_best')
        self.methodology = yaml_line_reader(info, 'optimizer', 'genetic_algorithm')
        self.code_interface = yaml_line_reader(info, 'code_type', 'PARCS342')
        self.calculation_type = yaml_line_reader(info, 'calc_type', 'single_cycle')
        
    ## Optimization Block ##
        try:
            info = self.file_settings['optimization']
        except KeyError:
            info = None
        
        self.population_size = yaml_line_reader(info, 'population_size', 1)
        self.num_generations = yaml_line_reader(info, 'number_of_generations', 1)
        self.symmetry = yaml_line_reader(info, 'solution_symmetry', 'octant')
        self.objectives = yaml_line_reader(info, 'objectives', None)
        termination_criteria_default = {'method':'None','termination_generations':0}
        self.termination_criteria = yaml_line_reader(info, 'termination_criteria', termination_criteria_default)
        
    ## Algorithm Block ##
        try:
            info = self.file_settings['algorithm']
        except KeyError:
            info = None
        
        selection_default = {'fitness':'weighted','method':'tournament'}
        self.selection = yaml_line_reader(info, 'selection', selection_default)
        self.reproducer = yaml_line_reader(info, 'reproducer', 'standard')
        self.mutation_type = yaml_line_reader(info, 'mutation_type', 'mutate_by_gene')
        self.mutation_rate = yaml_line_reader(info, 'mutation_rate', 0.5)
        crossover_default = {'method':'one_point','crossover_rate': 0.5, 'num_swaps': 1}
        self.crossover = yaml_line_reader(info, 'crossover', crossover_default)
        self.elites = yaml_line_reader(info, 'elites', 0)
        acquisition_function = yaml_line_reader(info, 'acquisition_function', 'LCB')
        #LCB and UCB acquisition functions will benefit more from an EE factor, and other acq functions may be better with a value of zero here
        if self.acquisition_function in ['LCB','UCB']:
            self.exploration_exploitation_factor = yaml_line_reader(info, 'exploration_exploitation_factor', 1.96)
        else:
            self.exploration_exploitation_factor = yaml_line_reader(info, 'exploration_exploitation_factor', 0)
        self.kernel_smoothness = yaml_line_reader(info, 'kernel_smoothness_factor', 0.5)
        
    ## Fuel Assembly Block ##
        self.fa_options = yaml_line_reader(self.file_settings, 'assembly_options', None)
        if not self.fa_options and self.code_interface not in ['nuscale_database']:
            raise ValueError("Assembly options must be nested with reflectors, fuels, and/or blankets with their parameters.")
        for param in ['cost_fuelcycle','av_fuelenrichment']:
            if param in self.objectives:
                for key in self.fa_options['fuel'].keys():
                    if not 'enrichment' in self.fa_options['fuel'][key] and \
                       not 'hm_loading' in self.fa_options['fuel'][key]:
                        raise ValueError(f"Entry for 'enrichment' or 'HM_loading' missing for fuel type '{key}'. This is required by the '{param}' objective.")
                for key in self.fa_options['blankets'].keys():
                    if not 'enrichment' in self.fa_options['blankets'][key] and \
                       not 'hm_loading' in self.fa_options['blankets'][key]:
                        raise ValueError(f"Entry for 'enrichment' or 'HM_loading' missing for blanket type '{key}'. This is required by the '{param}' objective.")
        
    ## Genome Block ##
        try:
            info = self.file_settings['decision_variables']
        except KeyError:
            info = None
        
        self.genome = yaml_line_reader(info, 'parameters', None)
        self.batches = yaml_line_reader(info, 'batches', None)
        #check that decision variable options are valid.
        if not self.genome:
            raise ValueError("'Parameters' must be specified in Decision Variables.")
        if self.calculation_type == 'eq_cycle' and not self.batches:
            raise ValueError("'Batches' must be specified in Decision Variables for the 'EQ Cycle' type.")
        for key, value in self.genome.items():
            if key not in self.fa_options['fuel']:
                raise ValueError(f"Decision variable option '{key}' not found in the list of fuel types under 'assembly_options'.")
        
    ## Calculation Block ##
        try:
            if self.code_interface in ["parcs342","parcs343"]:
                info = self.file_settings['parcs_data']
            elif self.code_interface == "nuscale_database":
                info = self.file_settings['nucale_data']
        except KeyError:
            info = None
        
        try:
            infomap = info['map']
        except:
            infomap = None
        
        self.code_walltime = yaml_line_reader(info, 'exec_walltime', 600)
        self.nrow = yaml_line_reader(infomap, 'num_rows', 17)
        self.ncol = yaml_line_reader(infomap, 'num_cols', 17)
        self.num_assemblies = yaml_line_reader(infomap, 'number_assemblies', 193)
        self.map_size = yaml_line_reader(infomap, 'core_symmetry', 'full')
        self.xs_lib = yaml_line_reader(info, 'xs_library_path', '../../') #as a relative path, this assumes the needed cross sections are in the base directory for the job.
        self.xs_extension = yaml_line_reader(info, 'xs_extension', '')
        self.power = yaml_line_reader(info, 'power', 3800.0)
        self.flow = yaml_line_reader(info, 'flow', 18231.89)
        self.inlet_temp = yaml_line_reader(info, 'inlet_temperature', 565.0)
        self.th_fdbk = yaml_line_reader(info, 'th_fdbk', True)
        self.pin_power_recon = yaml_line_reader(info, 'pin_power_recon', True)
        self.number_axial = yaml_line_reader(info, 'num_axial_nodes', 19)
        self.axial_nodes = yaml_line_reader(info, 'axial_nodes', [16.12, "15*25.739", 16.12])
        self.boc_exposure = yaml_line_reader(info, 'boc_core_exposure', 0.0)
        self.depl_steps = yaml_line_reader(info, 'depletion_steps', [1, 1, 30, 30, 30, 30, 30, 30])
        #NuScale database verification block
        if self.code_interface == 'nuscale_database':
            #Force octant symmetry for NuScale database
            if self.symmetry != 'octant':
                logging.warning(f'Core symmetry has been changed from {self.symmetry} to octant. NuScale database only supports octant symmetry.')
                self.symmetry == 'octant'
            
            #Verify assembly map length for each parameter in input file
            for parameter in self.genome:
                if len(self.genome[parameter]['map']) != 8:
                    raise ValueError(f'Parameter {parameter} has a map length of {len(parameter['map'])}, but needs length of 8')
            
            #Verify that the type parameter for each assembly is between 2-7, as these are the only assemblies available
            for assembly in self.fa_options['fuel']:
                if int(self.fa_options['fuel'][assembly]['type']) not in [2, 3, 4, 5, 6, 7]:
                    raise ValueError(f'Assembly {assembly} parameter "type" is incorrect. For NuScale database, types 2-7 exist')
                
        return