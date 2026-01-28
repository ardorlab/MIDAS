## Import Block ##
import yaml
import re
import logging
from pathlib import Path
from midas.utils import problem_preparation
"""
Classes for parsing and cleansing input data from the user-specified '.yaml' MIDAS input file.

Created by Nicholas Rollins. 09/11/2024
"""


## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

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
        if value not in ["genetic_algorithm","bayesian_optimization","simulated_annealing"]:
            raise ValueError("Requested methodology '" + value + "' invalid.")
    
    elif keyword == 'code_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["parcs342", "parcs343", "nuscale_database", "trace50p5", "polaris624","serpent","custom_function","styblinski_tang","listsum"]:
            raise ValueError("Code types currently supported: PARCS342, PARCS343, NuScale_Database, TRACE50p5.")
    
    elif keyword == 'calc_type':
        value = str(value).lower().replace(' ','_')
        if value not in ["single_cycle","eq_cycle", "lattice_physics", "continuous_variable"]:
            raise ValueError("Data type not supported.")
    
    elif keyword == 'input_template':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key =='apply':
                    new_item = item
                    if not isinstance(new_item, bool):
                        raise ValueError("'apply' flag for input template must be boolean")
                    if new_item:
                        logger.warning("Input template functionality currently only supports single cycle calculations") #TODO add other calculation functionality
                if new_key == 'loc':
                    new_item = Path(str(item))
                new_dict[new_key] = new_item

            if 'apply' in new_dict.keys() and new_dict['apply']:
                if 'loc' not in new_dict.keys():
                    raise ValueError("'apply' in input_template is set to true but path to template is not specified") 
            if 'loc' in new_dict.keys() and 'apply' not in new_dict.keys(): 
                logger.warning("path to input template is specified but 'apply' flag is not given. MIDAS will assume no input tempate in calculation")
                new_dict['apply'] = False
            return new_dict

    elif keyword == 'statistics_plots':
        value = str(value).lower().replace(' ','')
        if value == 'true':
            value = True
        elif value == 'false':
            value = False
        else:
            raise ValueError("statistics_plots only takes true or false as entries")
    
    elif keyword == 'convergence_plot':
        value = str(value).lower().replace(' ','')
        if value == 'true':
            value = True
        elif value == 'false':
            value = False
        else:
            raise ValueError("convergence_plot only takes true or false as entries")
    
    elif keyword == 'initial_population':
        if value:
            value = Path(str(value))
    
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
        if value not in ['octant','quarter','full', 'diagonal']:
            raise ValueError("Symmetry of the solution list must be octant, quarter, full, or diagonal.")
        
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
                                   'pxyz',
                                   'pxy',
                                   'cycle_length',
                                   'assembly_burnup',
                                   'cost_fuelcycle',
                                   'av_fuelenrichment',
                                   'maxcladtemp',
                                   'maxfueltemp',
                                   'maxgapq',
                                   'peak_reactivity',
                                   'max_critical_exposure',
                                   'keff_min',
                                   'keff_max',
                                   'keff_diff',
                                   'cpr',
                                   'chfr',
                                   'lhgr',
                                   'aplhgr',
                                   'max_doserate',
                                   'total_mass',
                                   'doppler_temperature_coefficient',
                                   'function_output',
                                   'list_sum']:
                    raise ValueError(f"Requested objective/constraint '{key}' not supported.")
                if new_key == 'aplhgr':
                    logger.warning("APLHGR requires 3d plotting of pin reconstruction.")
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
                        elif new_subkey == 'critical_power':
                            new_subitem = float(subitem)
                        elif new_subkey == 'linear_power':
                            new_subitem = float(subitem)
                        elif new_subkey == 'equation':
                            if new_key != 'function_output':
                                raise ValueError(f"'equation' option is only available under 'function_output' variable, but is under '{new_key}' instead")
                            new_subitem = str(subitem)
                        if new_key == 'cpr' and 'critical_power' not in item.keys():
                            raise ValueError(f"Critical power ratio is requested in objectives but the critical power is not provided.")
                        if new_key == 'lhgr' and 'linear_power' not in item.keys():
                            raise ValueError(f"Linear heat generation rate requested in objectives but the linear power density is not provided.")
                        if new_key == 'aplhgr' and 'linear_power' not in item.keys():
                            raise ValueError(f"Average palnar linear heat generation rate requested in objectives but the linear power density is not provided.")
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

            #check parameters logic
            if new_dict['method'] == 'ktournament' and 'k' not in new_dict.keys():
                new_dict['k'] = 4
                logger.warning("'k' parameter is missing from input while 'ktournament' selection method is used, k has been set to default value of 4.")
            elif new_dict['method'] == 'uniform' and 'crossover_rate' not in new_dict.keys():
                new_dict['crossover_rate'] = 0.50
                logger.warning("'crossover_rate' parameter is missing from input while 'uniform' selection method is used; crossover_rate has been set to default value of 0.50.")
            elif new_dict['method'] == 'random_element' and 'num_swaps' not in new_dict.keys():
                raise ValueError("'num_swaps' parameter is missing from input while 'random_element' selection method is used.")
                
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
    
    elif keyword == 'hyperparameter_convergence_criteria':
        value = float(value)
        if value <0:
            raise ValueError('Hyperparameter convergence criteria must be positive.')
    
    elif keyword == 'surrogate_off_generation':
        value = int(value)
        if value < 0:
            raise ValueError('Generation for turning the surrogate model fitting off must be greater than 0.')
        
    elif keyword == 'temperature':
        value = float(value)
        if value < 0.0:
            raise ValueError("Simulated Annealing intial temperature must be greater than 0")
        
    elif keyword == 'cooling_schedule':
        value = str(value).lower().replace(' ','_')
        if value not in ["exponential_decrease", "linear_update", "log_update", "lam"]:
            raise ValueError(f"Cooling schedule '{value}' not supported.")

    elif keyword == 'secondary_cooling_schedule':
        value = str(value).lower().replace(' ','_')
        if value not in ["exponential_decrease", "linear_update", "log_update", "none"]:
            raise ValueError(f"Secondary cooling schedule '{value}' not supported.")
        
    elif keyword == 'update_factor':
        value = float(value)
        if value < 0.0 or value > 1.0:
            raise ValueError("update factor for exponential cooling schedule must be 2.0 > alpha > 1.0")
        
    elif keyword == 'quality_factor':
        value = float(value)
        if value < 1.0 or value > 2.0:
            raise ValueError("quality factor for LAM cooling schedule must be 2.0 > qf > 1.0")
        
    elif keyword == 'scaling_factor':
        value = float(value)
        if value < 1.0 or value > 2.0:
            raise ValueError("scaling factor for LAM cooling schedule must be 2.0 > sf > 1.0")
        
    elif keyword == 'perturbation_type':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key == 'method':
                    new_item = str(item).lower().replace(' ','_')
                    if new_item not in ["perturb_by_gene"]:
                        raise ValueError(f"Requested perturbation method '{item}' not supported.")
                elif new_key =='num_perturbations':
                    new_item = int(item)
                    if item < 1:
                        raise ValueError("num_perturbations must be 1 or greater")
                new_dict[new_key] = new_item
            if 'num_perturbations' not in new_dict.keys():
                new_dict['num_perturbations'] = 1
            return new_dict
 
    elif keyword == 'buffer_size':
        value = int(value)
        if not isinstance(value, int):
            raise ValueError("bufer size must be an integer")
        elif isinstance(value, int) and value < 1:
            raise ValueError("buffer size must be greater than 1")


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
                                        if new_subsubitem not in ['all','radial','top','bottom']:
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
    
## Fuel Pin Parts Block ##
    elif keyword == 'rod_options':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                #check rod_geometries options
                if new_key == 'rod_geometries':
                    new_item = {}
                    if isinstance(item, dict):
                        #read rod types
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types, radii, and materials
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'type':
                                        new_subsubitem = str(subsubitem)
                                    elif new_subsubkey == 'radii':
                                        new_subsubitem = [float(x) for x in re.split(r'[, ]',str(subsubitem).strip('[]')) if x]
                                    elif new_subsubkey == 'materials':
                                        new_subsubitem = [str(x) for x in re.split(r'[, ]',str(subsubitem).strip('[]')) if x]
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested rod type missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Rod option missing rod_geometries and their parameters.")
                #check control_rods options
                elif new_key == 'control_rods':
                    new_item = {}
                    if isinstance(item, dict):
                        #read control rod types
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types, radii, and materials
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'type':
                                        new_subsubitem = str(subsubitem)[0]
                                    elif new_subsubkey == 'radii':
                                        new_subsubitem = [float(x) for x in re.split(r'[, ]',str(subsubitem).strip('[]')) if x]
                                    elif new_subsubkey == 'materials':
                                        new_subsubitem = [str(x) for x in re.split(r'[, ]',str(subsubitem).strip('[]')) if x]
                                    elif new_subsubkey == 'guide_tube':
                                        new_subsubitem = str(subsubitem).upper()
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested control rod type missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        if item: #optional argument can be 'None'.
                            raise ValueError("Rod option missing control_rods and their parameters.")
                #check compositions options
                elif new_key == 'compositions':
                    new_item = {}
                    if isinstance(item, dict):
                        #read compositions
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check types and values
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'type':
                                        new_subsubitem = str(subsubitem).upper().strip()
                                    elif new_subsubkey == 'values':
                                        new_subsubitem = [str(x) for x in re.split(r'[, ]',str(subsubitem).strip('[]')) if x]
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested composition type missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Rod option missing compositions and their parameters.")
                #check materials options
                elif new_key == 'materials':
                    new_item = {}
                    if isinstance(item, dict):
                        #read materials
                        for subkey, subitem in item.items():
                            new_subkey = str(subkey)
                            new_subitem = {}
                            if isinstance(subitem, dict):
                                #check comps, dens, and temps #!TODO: additional properties not supported?
                                for subsubkey, subsubitem in subitem.items():
                                    new_subsubkey =str(subsubkey).lower().replace(' ','_')
                                    if new_subsubkey == 'comp':
                                        new_subsubitem = str(subsubitem).lower().strip()
                                    elif new_subsubkey == 'dens':
                                        new_subsubitem = float(subsubitem)
                                    elif new_subsubkey == 'temp':
                                        new_subsubitem = float(subsubitem)
                                    elif new_subsubkey == 'fueltype':
                                        new_subsubitem = bool(subsubitem)
                                    new_subitem[new_subsubkey] = new_subsubitem
                            else:
                                raise ValueError("Requested composition type missing parameters.")
                            new_item[new_subkey] = new_subitem
                        new_dict[new_key] = new_item
                    else:
                        raise ValueError("Rod option missing compositions and their parameters.")
            #check parameters logic
            if 'materials' in new_dict:
                for key, value in new_dict['materials'].items():
                    if not 'fueltype' in value:
                        new_dict['materials'][key]['fueltype'] = False
            #!TODO: make sure each rod type has type, radii, materials
            #!TODO: make sure each composition has type, values
            #!TODO: make sure each material has comp and density (temp is optional)
            return new_dict

## Genome Block ##
    elif keyword in ['parameters']:
        logger.warning(f"Handeling continuous variables in MIDAS is still in development and continuous variable optimization may not work as indented. ")
        new_dict = {}
        if isinstance(value, dict):
            for key, item in value.items():
                new_key = str(key)
                new_dict[new_key] = {}
                #check decision variable options
                if isinstance(value[key], dict):
                    for subkey, subitem in item.items():
                        new_subkey = str(subkey).lower()
                        if new_subkey in ['continuous_range','discrete_range']:
                            if not isinstance(subitem, list):
                                raise ValueError(f"Entry '{new_subkey}' under decision variable '{new_key}' must be a list of two numbers.")
                            for rangebound in subitem: 
                                if (not isinstance(rangebound, float) and not isinstance(rangebound, int)) or isinstance(rangebound, bool): 
                                    raise ValueError(f"Entry '{new_subkey}' values under decision variable '{new_key}' must be two numeric values in ascending order.")
                            if len(subitem) != 2: 
                                raise ValueError(f"Entry '{new_subkey}' list under decision variable '{new_key}' must contain two numeric values in ascending order.")
                            if subitem[0] > subitem[1]: 
                                raise ValueError(f"Entry '{new_subkey}' list under decision variable '{new_key}' must contain two numeric values in ascending order.")
                            new_dict[new_key][new_subkey] = subitem
                        if new_subkey == "increment":
                            try:
                                if isinstance(subitem, list):
                                    if len(subitem) == 0:
                                        raise ValueError(f"The {new_subkey} entry is a list of length 0. It must either be a list with one or more entries for a non-uniform range, or a single number for a uniform range")
                                    for index in range(0, len(subitem)):
                                        subitem[index] = float(subitem[index])
                                        if subitem[index] < 0: 
                                            raise ValueError("Continuous variable 'increment' entries must be greater than 0")
                                    new_dict[new_key][new_subkey] = subitem
                                else:
                                    new_dict[new_key][new_subkey] = float(subitem)
                                    if new_dict[new_key][new_subkey] < 0: 
                                        raise ValueError("Continuous variable 'increment' must be greater than 0")
                            except TypeError:
                                raise ValueError(f"Subkey {new_subkey} has entry of type {type(subitem)} but only accepts a list of integers/floats or a single integer/float")
                        if new_subkey == "index":
                            new_dict[new_key][new_subkey] = int(subitem)
                    if "discrete_range" not in item.keys() and "continuous_range" not in item.keys():
                        raise ValueError(f"{item} variable must include either a 'discrete_range' or 'continuous_range' entry.")
                    if "discrete_range" in item.keys() and "increment" not in item.keys():
                        raise ValueError(f"Increment for {item} variable with 'discrete_range' is not provided.")
                    elif "discrete_range" in item.keys(): 
                        if (item["increment"] / (item["discrete_range"][1] - item["discrete_range"][0]))*100 >= 10:
                            logger.warning(f"Continuous variable increment for '{new_key}' is large relative to the range. Is this intentional?")
        return new_dict
    
    elif keyword in ['lattice_parameters', 'assembly_parameters', 'batches']:
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
    ## PARCS DATA ##
    elif keyword == 'core_type':
        value = str(value).upper()
        if value not in ["PWR", "BWR"]:
            raise ValueError(f"Requested core type '{value}' not supported. \n          Supported values include: PWR, BWR")
        if value == "bwr":
            logger.warning("functionality for BWR optimization is still under development")
            #TODO BWR parcs input generator
            #TODO core flow optimization
            #TODO control blade sequence optimization
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

    elif keyword == 'assembly_pitch':
        value = float(value)
    
    elif keyword == 'core_symmetry':
        value = str(value).lower()
        if value not in ['full','quarter']:
            raise ValueError("Requested core symmetry (used for printing) not valid.")
    
    elif keyword == 'xs_library_path':
        value = Path('../../') / Path(str(value))
    
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

    elif keyword == 'coolant_density':
        value = float(value)

    elif keyword == 'fuel_temperature':
        value = float(value)
    
    elif keyword == 'th_fdbk':
        if isinstance(value, dict):
            new_dict = {}
            for key, item in value.items():
                new_key = str(key).lower()
                if new_key =='apply':
                    new_item = item
                    if not isinstance(new_item, bool):
                        raise ValueError("'apply' flag for th_fdbk must be boolean")
                if new_key == 'loc':
                    new_item = str(item).lower()
                    if new_item == 'none':
                        new_item = None
                    else:
                        new_item = Path(str(item))
                new_dict[new_key] = new_item

            if 'apply' in new_dict.keys() and new_dict['apply']:
                if 'loc' not in new_dict.keys():
                    logger.warning("'apply' in th_fdbk is set to true but path to paths input is not specified. PARCS internal mass/energy balance solver is assumed.")
                    new_dict['loc'] = None
            if 'loc' in new_dict.keys() and 'apply' not in new_dict.keys(): 
                logger.warning("Path to PATHS input is specified but 'apply' flag is not given. MIDAS will assume mass/energy balance solver in calculation")
                new_dict['apply'] = False
            return new_dict
    
    elif keyword == 'pin_power_recon':
        value = bool(value)

    elif keyword == 'assembly_pins':
        value = int(value)
    
    elif keyword == 'assembly_guide_tubes':
        value = int(value)

    elif keyword == 'pin_dimensions':
        value = str(value)
    
    elif keyword == 'num_axial_nodes':
        value = int(value)
    
    elif keyword == 'axial_nodes':
        if isinstance(value, str):
            value = [a.strip() for a in value.strip().replace(', ',',').replace(' ',',').split(',')]
        new_value = []
        for node in value:
            if "*" in str(node):
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
        value = [x.strip() for x in re.split(r'[, ]',str(value).strip('[]')) if x]
        if isinstance(value, list):
            new_value = []
            for step in value:
                if "*" in str(step):
                    s_step = step.split('*')
                    new_value.extend(int(s_step[0])*[float(s_step[1])])
                else:
                    new_value.append(float(step))
        return new_value
    
    elif keyword == 'equilibrium_cycles':
        value = int(value)

    ## TRACE DATA ##
    elif keyword == 'initialize_code':
        value = str(value).lower().replace(' ','_')
        if value not in ["parcs343"]:
            raise ValueError("Requested code for initializing TRACE calculation not supported. Supported codes include: PARCS343.")
    
    elif keyword == 'ss_input_file':
        if value:
            value = Path(str(value))
            if not value.exists():
                raise ValueError(f"Could not locate TRACE steady-state input file: '{value}'.")
        else:
            value = None
    
    elif keyword == 'tr_input_file':
        if value:
            value = Path(str(value))
            if not value.exists():
                raise ValueError(f"Could not locate TRACE transient input file: '{value}'.")
        else:
            value = None
    
    elif keyword == 'maptab_file':
        if value:
            value = Path(str(value))
            if not value.exists():
                raise ValueError(f"Could not locate MAPTAB file for TRACE-PARCS coupling: '{value}'.")
        else:
            value = None
    
    elif keyword == 'ss_power_fraction':
        value = float(value)
        if value < 0.0:
            value *= -1
        if value <= 1.0:
            value *= 100 #change from fraction to percent
    
    ## POLARIS DATA ##
    elif keyword == 'exec_walltime':
        value = int(value)
        if value <= 0:
            raise ValueError("'exec_walltime' must be a positive number, measured in seconds.")
    
    elif keyword == 'system_type':
        value = str(value).upper().strip()
    
    elif keyword == 'xs_library':
        value = str(value).upper().strip()

    elif keyword == 'num_rows':
        value = int(value)
    
    elif keyword == 'num_cols':
        value = int(value)

    elif keyword == 'pin_pitch':
        value = float(value)

    elif keyword == 'lattice_symmetry':
        value = str(value).upper().strip()
        if value not in ['SE', 'FULL', 'DIAGONAL']:
            raise ValueError("'lattice_symmetry' must be either 'SE', 'FULL', or 'DIAGONAL'.")
        
    elif keyword == 'box':
        value = str(value)

    elif keyword == 'hgap':
        value = str(value)
    
    elif keyword == 'power':
        value = float(value)
    
    elif keyword == 'bulk_temperatures':
        value = float(value)
    
    elif keyword == 'fuel_temperatures':
        value = float(value)
    
    elif keyword == 'controlrods_inserted':
        value = bool(value)
    
    elif keyword == 'borated_material':
        if value:
            matname, ppm = map(str,[x for x in re.split(r'[ ,]',str(value).strip('[]')) if x])
            value = [str(matname),int(ppm)] #material name, boron concentration in ppm
    
    elif keyword == 'num_mesh_rings':
        value = int(value)
    
    elif keyword=='depletion_steps':
        value = [float(x) for x in re.split(r'[, ]',str(value).strip('[]')) if x]
    
    return value

def parcs343_template_check(self):
    """
    Checks to ensure that necessary flags are present if a template input file is provided for a code interface. 
    The check is purposefully minimal so that users are able to run the code exactly as they wish.
    Ensuring that the template is exactly as intended is up to the user.
    Written by Jake Mikouchi. 3/24/2025
    """
    if self.code_interface.lower() == 'parcs343':
        necessary_flags = {'caseid': False, 'cntl': False, 'param': False, 'geom': False, 'fdbk': False, 
                        'th': False, 'depl':False}
        with open(self.input_template['loc'], "r") as file:
            lines = file.readlines()  
            for line in lines:
                if "caseid" in line.lower() and "!" not in line.lower():
                    necessary_flags['caseid'] = True
                if "cntl" in line.lower() and "!" not in line.lower():
                    necessary_flags['cntl'] = True
                if "param" in line.lower() and "!" not in line.lower():
                    necessary_flags['param'] = True
                if "geom" in line.lower() and "!" not in line.lower():
                    necessary_flags['geom'] = True
                if "fdbk" in line.lower() and "!" not in line.lower():
                    necessary_flags['fdbk'] = True
                if "th" in line.lower() and "!" not in line.lower():
                    necessary_flags['th'] = True
                if "depl" in line.lower() and "!" not in line.lower():
                    necessary_flags['depl'] = True

    if self.code_interface.lower() == 'parcs342':
        necessary_flags = {'caseid': False, 'cntl': False, 'param': False, 'geom': False, 'fdbk': False, 
                        'th': False, 'depl':False}
        with open(self.input_template['loc'], "r") as file:
            lines = file.readlines()  
            for line in lines:
                if "caseid" in line.lower() and "!" not in line.lower():
                    necessary_flags['caseid'] = True
                if "cntl" in line.lower() and "!" not in line.lower():
                    necessary_flags['cntl'] = True
                if "param" in line.lower() and "!" not in line.lower():
                    necessary_flags['param'] = True
                if "geom" in line.lower() and "!" not in line.lower():
                    necessary_flags['geom'] = True
                if "fdbk" in line.lower() and "!" not in line.lower():
                    necessary_flags['fdbk'] = True
                if "th" in line.lower() and "!" not in line.lower():
                    necessary_flags['th'] = True
                if "depl" in line.lower() and "!" not in line.lower():
                    necessary_flags['depl'] = True

    if self.code_interface.lower() == "nuscale_database":
        raise ValueError(f'input templates are not supported for nuscale_database')

    for key, value in necessary_flags.items():
        if not necessary_flags[key]:
            raise ValueError(f'{self.code_interface} input template is missing {key} flag') 



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
        self.code_interface = yaml_line_reader(info, 'code_type', 'PARCS343')
        template_default = {'apply':False,'loc':''}
        self.input_template = yaml_line_reader(info, 'input_template', template_default)
        self.calculation_type = yaml_line_reader(info, 'calc_type', 'single_cycle')
        self.statistics_plots = yaml_line_reader(info, 'statistics_plots', True)
        self.convergence_plots = yaml_line_reader(info, 'convergence_plots', True)
        self.initial_population = yaml_line_reader(info, 'initial_population', None)
        if self.input_template['apply'] and self.code_interface == 'parcs343':
            parcs343_template_check(self)
    ## Optimization Block ##
        try:
            info = self.file_settings['optimization']
        except KeyError:
            info = None
        
        self.population_size = yaml_line_reader(info, 'population_size', 1)
        if self.methodology == 'simulated_annealing' and self.num_procs <= 1:
            self.population_size = 1
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
        self.acquisition_function = yaml_line_reader(info, 'acquisition_function', 'LCB')
        #LCB and UCB acquisition functions will benefit more from an EE factor, and other acq functions may be better with a value of zero here
        if self.acquisition_function in ['LCB','UCB']:
            self.exploration_exploitation_factor = yaml_line_reader(info, 'exploration_exploitation_factor', 1.96)
        else:
            self.exploration_exploitation_factor = yaml_line_reader(info, 'exploration_exploitation_factor', 0)
        self.kernel_smoothness = yaml_line_reader(info, 'kernel_smoothness_factor', 0.5)
        self.kernel_hyperparam_conv = yaml_line_reader(info, 'hyperparameter_convergence_criteria', 0.01)
        self.surrogate_fitting_off = yaml_line_reader(info, 'surrogate_off_generation', int(self.num_generations/2))
        self.initial_temperature = yaml_line_reader(info, 'temperature', 100)
        if self.num_procs > 1:
            self.cooling_schedule = yaml_line_reader(info, 'cooling_schedule', 'lam')
        if self.num_procs <= 1:
            self.cooling_schedule = yaml_line_reader(info, 'cooling_schedule', 'exponential_decrease')
            if self.cooling_schedule == "lam":
                raise ValueError(f"Cooling schedule '{self.cooling_schedule}' only available for parallel simulated annealing.")
        self.secondary_cooling_schedule = yaml_line_reader(info, 'secondary_cooling_schedule', 'exponential_decrease')
        self.update_factor = yaml_line_reader(info, 'update_factor', 0.95)
        self.quality_factor = yaml_line_reader(info, 'quality_factor', 1.1)
        self.scaling_factor = yaml_line_reader(info, 'scaling_factor', 1.5)
        perturbation_default = {'method':'perturb_by_gene','num_perturbations':1}
        self.perturbation_type = yaml_line_reader(info, 'perturbation_type', perturbation_default)
        self.buffer_size = yaml_line_reader(info, 'buffer_size', 10)

        
    ## Fuel Assembly Block ##   
        self.fa_options = yaml_line_reader(self.file_settings, 'assembly_options', None)
        if not self.fa_options and self.code_interface not in ['nuscale_database','polaris624','serpent','custom_function','styblinski_tang']:
            raise ValueError("Assembly options must be nested with reflectors, fuels, and/or blankets with their parameters.")
        if self.calculation_type in ['single_cycle','eq_cycle']:
            for param in ['cost_fuelcycle','av_fuelenrichment']:
                if param in self.objectives:
                    for key in self.fa_options['fuel'].keys():
                        if not 'enrichment' in self.fa_options['fuel'][key] and \
                           not 'hm_loading' in self.fa_options['fuel'][key]:
                            raise ValueError(f"Entry for 'enrichment' or 'HM_loading' missing for fuel type '{key}'. This is required by the '{param}' objective.")
                    if 'blankets' in self.fa_options:
                        for key in self.fa_options['blankets'].keys():
                            if not 'enrichment' in self.fa_options['blankets'][key] and \
                               not 'hm_loading' in self.fa_options['blankets'][key]:
                                raise ValueError(f"Entry for 'enrichment' or 'HM_loading' missing for blanket type '{key}'. This is required by the '{param}' objective.")
        
    ## Fuel Pin Parts Block ## (for lattice physics calcs)
        self.pin_options = yaml_line_reader(self.file_settings, 'rod_options', None)
        if not self.pin_options and self.calculation_type in ['lattice_physics']:
            raise ValueError("Fuel pin options must be nested with rod_geometries, compositions, and/or controls_rods.")
        
        
    ## Genome Block ##
        try:
            info = self.file_settings['decision_variables']
        except KeyError:
            info = None
        
        if self.calculation_type in ['single_cycle','eq_cycle','listsum']:
            self.genome = yaml_line_reader(info, 'assembly_parameters', None)
        elif self.calculation_type in ['lattice_physics']:
            self.genome = yaml_line_reader(info, 'lattice_parameters', None)
        elif self.calculation_type in ['continuous_variable']:
            logger.warning("'parameters' decision variable is reserved for continous variables")
            self.genome = yaml_line_reader(info, 'parameters', None)
            #Create a list of possible values a gene can take for discrete ranges
            self.genome = problem_preparation.Prepare_Problem_Values.prepare_discrete_range(self.genome)
            #Normalize all ranges for continuous variables
            self.genome = problem_preparation.Prepare_Problem_Values.normalize_continuous_variables(self.genome)

        self.batches = yaml_line_reader(info, 'batches', None)
        #check that decision variable options are valid.
        if not self.genome:
            raise ValueError("'assembly_parameters', 'lattice_parameters', or 'parameters' must be specified in Decision Variables.")
        if self.calculation_type == 'eq_cycle' and not self.batches:
            raise ValueError("'Batches' must be specified in Decision Variables for the 'EQ Cycle' type.")
        for key, value in self.genome.items():
            if self.fa_options:
                if key not in self.fa_options['fuel']:
                    raise ValueError(f"Decision variable option '{key}' not found in the list of fuel types under 'assembly_options'.")
            elif self.pin_options:
                if key not in self.pin_options['rod_geometries']:
                    raise ValueError(f"Decision variable option '{key}' not found in the list of rod types under 'rod_options'.")
        
    ## Calculation Block ## #!TODO: should each parameter set be nested under a code-specific object?
        info = None; THinfo = None; infomap = None # initialize info variables
        if self.code_interface not in ['custom_function','styblinski_tang','listsum']:
            try:
                if self.code_interface in ["parcs342","parcs343"]:
                    info = self.file_settings['parcs_data']
                elif self.code_interface == "nuscale_database":
                    info = self.file_settings['nuscale_data']
                elif self.code_interface == "serpent":
                    info = self.file_settings['serpent_data']
                elif self.code_interface == "trace50p5": #multiphysics calcs must first be initialized in neutronics code.
                    try:
                        info = self.file_settings['parcs_data']
                    except KeyError:
                        pass
                    THinfo = self.file_settings['trace_data']
                elif self.code_interface == "polaris624":
                    info = self.file_settings['polaris_data']
                try:
                    infomap = info['map']
                except KeyError:
                    pass
            except KeyError:
                pass
        
        # PARCS input block
        self.core_type = yaml_line_reader(info, 'core_type', "PWR")
        self.code_walltime = yaml_line_reader(info, 'exec_walltime', 600)
        self.nrow = yaml_line_reader(infomap, 'num_rows', 17)
        self.ncol = yaml_line_reader(infomap, 'num_cols', self.nrow)
        self.num_assemblies = yaml_line_reader(infomap, 'number_assemblies', 193)
        self.assembly_pitch = yaml_line_reader(infomap, 'assembly_pitch', 21.50)
        self.map_size = yaml_line_reader(infomap, 'core_symmetry', 'full')
        self.xs_lib = yaml_line_reader(info, 'xs_library_path', './') #!TODO: interpret this path relative to the MIDAS job base dir, not opt indv base dir.
        self.xs_extension = yaml_line_reader(info, 'xs_extension', '')
        self.power = yaml_line_reader(info, 'power', 3800.0)
        self.flow = yaml_line_reader(info, 'flow', 18231.89)
        self.inlet_temp = yaml_line_reader(info, 'inlet_temperature', 565.0)
        self.coolant_dens = yaml_line_reader(info, 'coolant_density', 0.740)
        self.fuel_temp = yaml_line_reader(info, 'fuel_temperature', 900.0)
        th_default = {'apply':True, 'loc': None}
        self.th_fdbk = yaml_line_reader(info, 'th_fdbk', th_default)
        self.pin_power_recon = yaml_line_reader(info, 'pin_power_recon', True)
        self.assembly_pins = yaml_line_reader(info, 'assembly_pins', 264)
        self.assembly_guide_tubes = yaml_line_reader(info, 'assembly_guide_tubes', 25)
        self.pin_dimensions = yaml_line_reader(info, 'pin_dimensions', '4.1 4.75 0.58 6.13')
        self.number_axial = yaml_line_reader(info, 'num_axial_nodes', 19)
        self.axial_nodes = yaml_line_reader(info, 'axial_nodes', [16.12, "15*25.739", 16.12])
        self.boc_exposure = yaml_line_reader(info, 'boc_core_exposure', 0.0)
        self.depl_steps = yaml_line_reader(info, 'depletion_steps', [1, 1, 30, 30, 30, 30, 30, 30])
        self.equilibrium_cycles = yaml_line_reader(info, "equilibrium_cycles", 10)
        if (not self.pin_power_recon and 'pinpowerpeaking' in self.objectives.keys()) or (not self.pin_power_recon and 'fdeltah' in self.objectives.keys()):
            logger.warning('Pin power reconstruction is turned off but pin peaking factors are requested in objectives.')
            
        self.active_cycles = yaml_line_reader(info, "active_cycles", 500)
        self.inactive_cycles = yaml_line_reader(info, "inactive_cycles", 50)
        self.particles_per_history = yaml_line_reader(info, "particles_per_history", 5000)
        
        # TRACE input block
        if self.code_interface == "trace50p5":
            TRACE_file_defaults = {'templatefile':'./trace_ss.inp', 'maptabfile':'./TRACE-PARCS.map'}
        else:
            TRACE_file_defaults = {'templatefile':None, 'maptabfile':None}
        
        self.init_code = yaml_line_reader(THinfo, 'initialize_code', 'PARCS343')
        self.inp_template_ss = yaml_line_reader(THinfo, 'ss_input_file', TRACE_file_defaults['templatefile'])
        self.inp_maptabfile = yaml_line_reader(THinfo, 'maptab_file', TRACE_file_defaults['maptabfile'])
        self.ss_powerfraction = yaml_line_reader(THinfo, 'ss_power_fraction', 100)
        self.inp_template_tr = yaml_line_reader(THinfo, 'tr_input_file', None)
        
        # POLARIS input block
        if self.code_interface == "polaris624":
            self.nrow = yaml_line_reader(info, 'num_rows', 17)
            self.ncol = yaml_line_reader(info, 'num_cols', self.nrow)
            self.map_size = yaml_line_reader(info, 'lattice_symmetry', 'SE')
        self.system_type =  yaml_line_reader(info, 'system_type', 'PWR')
        self.xs_library = yaml_line_reader(info, 'xs_library', 'fine_therm') # May suffer from / need similar changes as PARCS XS library handling
        self.pin_pitch = yaml_line_reader(info, 'pin_pitch', 1.26)
        self.box = yaml_line_reader(info, 'box', '0')
        self.hgap = yaml_line_reader(info, 'hgap', '0')
        self.powdens = yaml_line_reader(info, 'powdens', 36) #W/gIHM
        self.bulk_temps = yaml_line_reader(info, 'bulk_temperatures', 566.0) #K
        self.fuel_temps = yaml_line_reader(info, 'fuel_temperatures', 900.0) #K
        self.cr_inserted = yaml_line_reader(info, 'controlrods_inserted', False)
        self.boronmat = yaml_line_reader(info, 'borated_material', None) #str, ppm
        self.num_meshrings = yaml_line_reader(info, 'num_mesh_rings', 3)
        self.depl_steps = yaml_line_reader(info, 'depletion_steps', [1, 1, 30, 30, 30, 30, 30, 30])   
        #NuScale database verification block
        if self.code_interface == 'nuscale_database':
            #Force octant symmetry for NuScale database
            if self.symmetry != 'octant':
                logger.warning(f'Core symmetry has been changed from {self.symmetry} to octant. NuScale database only supports octant symmetry.')
                self.symmetry == 'octant'
            
            #Verify assembly map length for each parameter in input file
            for parameter in self.genome:
                if len(self.genome[parameter]['map']) != 8:
                    map_length = len(parameter['map'])
                    raise ValueError(f'Parameter {parameter} has a map length of {map_length}, but needs length of 8.')
            
            #Verify that the type parameter for each assembly is between 2-7, as these are the only assemblies available
            for assembly in self.fa_options['fuel']:
                if int(self.fa_options['fuel'][assembly]['type']) not in [2, 3, 4, 5, 6, 7]:
                    raise ValueError(f'Assembly {assembly} parameter "type" is incorrect. For NuScale database, types 2-7 exist.')
        return