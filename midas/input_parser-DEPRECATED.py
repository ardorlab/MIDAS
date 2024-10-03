import yaml
import numpy as np
from copy import deepcopy
from midas.applications.parcs_332 import PARCS_Geometry_Funcs #this needed for the repurposed add_additional_information
"""
Classes for parsing and cleansing input data from the user-specified '.yaml' MIDAS input file.

Created by Nicholas Rollins. 09/11/2024
"""


class input_parser(object):
    """
    Centralized class for parsing user-supplied input arguments from the 
    MIDAS '.yaml' input file
    
    Created by Nicholas Rollins. 09/11/2024
    """
    def __init__(self,yaml_file):
        with open(yaml_file) as f:
            yaml_data = yaml.safe_load(f)
        
        self.parse_input_data(yaml_data)
        self.add_additional_information(yaml_data) #!should this be combined into parse_input_data?
        self.file_settings = yaml_data #!temporary redundancy while some input parsing is still being handled elsewhere.
        return
    
    def parse_input_data(self,data):
        """
        Interpret and verify parsed input data.
        
        Created by Nicholas Rollins. 09/11/2024
        """
        ## Optimization Block ##
        self.opt_method = str(data['optimization']['methodology']).lower()
        self.ext_model = str(data['optimization']['external_model']).lower()
        self.calculation_type = str(data['optimization']['data_type']).lower()
        self.chromosomes = data['genome']['chromosomes'] #!TODO: add validation of chromosomes input
        self.population_size = int(data['optimization']['population_size'])
        self.mutation_settings = data['optimization']['mutation'] #!TODO: add validation of mutation settings
        self.reproducer = str(data['optimization']['reproducer']).lower()
        self.selection_settings = data['optimization']['selection'] #!TODO: add validation of selection settings
        
        if self.opt_method == "simulated_annealing":
            self.SA_temperature = float(data['optimization']['cooling_schedule']['temperature'])
            self.SA_alpha = float(data['optimization']['cooling_schedule']['alpha'])
        
        if "fitness" in self.selection_settings:
            self.fitness_method = str(self.selection_settings['fitness']).lower()
        else:
            self.fitness_method = str(data['optimization']['fitness']).lower()
        
        ## Validate Input Values ##
        if self.opt_method not in ["genetic_algorithm",
                                   "simulated_annealing",
                                   "reinforcement_learning",
                                   "random_solutions",
                                   "lava"]:
            raise ValueError("Invalid entry in methodology. Valid options are genetic_algorithm, simulated_annealing, \
                              reinforcement_learning, random_solutions, or lava.")
        if self.ext_model in ['parcs','parcs332']:
            self.ext_model = "PARCS332" #Normalize input string
        elif self.ext_model in ['ncsu_core']: #!TODO: add accepted Casmo/SIMULATE strings
            self.ext_model = "NCSU_Core"
        else:
            raise ValueError("Invalid entry in external_model. Valid options are PARCS or NCSU_Core.")
        if self.calculation_type not in ["pin_lattice",
                                         "single_assembly_simulate",
                                         "loading_pattern",
                                         "loading_pattern_cnn",
                                         "loading_pattern_parcs332",
                                         "mcycle_loading_pattern_parcs332",
                                         "mcycle_grouped_loading_pattern_parcs332",
                                         "mcycle_inventory_loading_pattern_parcs332",
                                         "loading_patternsimple_parcs332",
                                         "fixed_loading_pattern"]:
            raise ValueError("Invalid entry in data_type. Valid options are pin_lattice, single_assembly_simulate, \
                              loading_pattern, loading_pattern_cnn, loading_pattern_parcs332, \
                              mcycle_loading_pattern_parcs332, mcycle_grouped_loading_pattern_parcs332, \
                              mcycle_inventory_loading_pattern_parcs332, loading_patternsimple_parcs332, \
                              or fixed_loading_pattern.")
        
        try:
            self.num_generations = int(data['optimization']['number_of_generations'])
        except ValueError:
            if data['optimization']['number_of_generations'] == "calculate_from_genes":
                self.num_generations = "calculate_from_genes"
            else:
                raise ValueError("Number_of_generations must be integer or 'calculate_from_genes'.")

        if self.reproducer not in ["multiple_genes_per_chromosome",
                                   "fixed_problem",
                                   "unique_genes",
                                   "mcycle",
                                   "standard"]:
            raise ValueError("Invalid entry in methodology. Valid options are multiple_genes_per_chromosome, \
                              fixed_problem, unique_genes, mcycle, or standard.")

        if self.fitness_method not in ["weighted",
                                       "ranked",
                                       "binned",
                                       "acdpf",
                                       "weighted_positive",
                                       "adaptive",
                                       "quantum"]:
            raise ValueError("Invalid entry in fitness. Valid options are weighted, ranked, binned, \
                              ACDPF, weighted_positive, adaptive, or quantum.")
        return

    def add_additional_information(self,data): #!TODO: this needs to be moved to the input_parser class.
        """
        Adds information on reactor parameters.
        This function replaces all equivalent functions in the other PARCS332 classes.

        Parameters
            data: The input data dictionary for the parameters.
        
        Written by Gregory Delipei. 01/08/2022
        Modified by Nicholas Rollins. 09/12/2024
        """
    ## Initialize core values ##
        self.nrow=17 #!TODO: this shouldn't be hardcoded to allow for other core shapes.
        self.ncol=17
        self.core_dict = {}
        
    ## Parse loading pattern ##
        self.symmetry = str(data["genome"]["parcs_data"]["symmetry"]).lower()
        if self.symmetry not in ["quarter","octant","full"]:
            raise ValueError(f"The selected symmetry ({self.input.symmetry}) is not valid.")
        
        self.input = deepcopy(self) #!needed to use generate_core
        self.core_dict['core_map'], self.core_dict['core_id'] = PARCS_Geometry_Funcs.generate_core(self)
        
    ## Parse Fuel Assembly types ##
        if 'reflectors' in data["genome"]["chromosomes"]:
            self.core_dict['Reflectors'] = {}
            for key,value in data["genome"]["chromosomes"]["reflectors"].items():
                refl_it = {}
                if value['refl_type'].lower() == "all":
                    refl_it['Tag'] = str(value['type'])
                    refl_it['Cross_Section']=value['serial']
                    self.core_dict['Reflectors']['Radial'] = refl_it
                    self.core_dict['Reflectors']['Top'] = refl_it
                    self.core_dict['Reflectors']['Bottom'] = refl_it
                    break
                elif value['refl_type'].lower() == "radial":
                    refl_it['Tag'] = str(value['type'])
                    refl_it['Cross_Section']=value['serial']
                    self.core_dict['Reflectors']['Radial'] = refl_it
                    continue
                elif value['refl_type'].lower() == "top":
                    refl_it['Cross_Section']=value['serial']
                    self.core_dict['Reflectors']['Top'] = refl_it
                    continue
                elif value['refl_type'].lower() == "bottom":
                    refl_it['Cross_Section']=value['serial']
                    self.core_dict['Reflectors']['Bottom'] = refl_it
                    continue
        self.core_dict['Blanket'] = None #set default value
        if 'blanket' in data["genome"]["chromosomes"]:
            self.core_dict['Blanket'] = {}
            for key,value in data["genome"]["chromosomes"]["blanket"].items():
                self.core_dict['Blanket']['Cross_Section'] = value['serial']
        if 'inventory' in list(data["genome"].keys()):
            self.core_dict['Inventory'] = data["genome"]['inventory']
        else:
            inventory = {}
            for key,value in data["genome"]["chromosomes"]["fuel"].items():
                inv_it = {}
                inv_it['Max_Limit']=np.inf
                inv_it['In_Design']=0
                inv_it['Tag']=str(value['type'])
                inv_it['Cross_Section']=value['serial']
                if 'Cost' in value:
                    inv_it['Cost']=float(value['Cost'])
                else:
                    inv_it['Cost']=250 #!units?
                inventory[key]=inv_it
            self.core_dict['Inventory'] = inventory
        ## Inventroy_Groups used by 'simple' loading pattern calc
            inventory_groups = {}
            for key,value in self.core_dict['Inventory'].items():
                inventory_groups[key]={'Values': [key],
                                       'Limit': 'Max',
                                       'Limit_Value':np.inf}
            self.core_dict['Inventory_Groups'] = inventory_groups
    
    ## Parse PARCS Data ##
        info = data['genome']['parcs_data']
        if 'xs_library' in info:
            self.library = info['xs_library']
        if 'power' in info:
            self.power = float(info['power'])
        if 'flow' in info:
            self.flow = float(info['flow'])
        if 'inlet_temperature' in info:
            self.inlet_temperature = float(info['inlet_temperature'])
        self.th_fdbk = True #set default value
        if 'th_fdbk' in info:
            self.th_fdbk = str(info['th_fdbk'])
            if self.th_fdbk.upper() in ['N','NO','F','FALSE','OFF']:
                self.th_fdbk = False
            elif self.th_fdbk.upper() in ['Y','YES','T','TRUE','ON']:
                self.th_fdbk = True
            else: #!TODO: warn user option wasn't recognized or error out?
                self.th_fdbk = True
        self.pin_power_recon = True #set default value
        if 'pin_power_recon' in info:
            self.pin_power_recon = str(info['pin_power_recon'])
            if self.pin_power_recon.upper() in ['N','NO','F','FALSE','OFF']:
                self.pin_power_recon = False
            elif self.pin_power_recon.upper() in ['Y','YES','T','TRUE','ON']:
                self.pin_power_recon = True
            else: #!TODO: warn user option wasn't recognized or error out?
                self.pin_power_recon = True
        if 'map_size' in info:
            self.map_size = info['map_size']
        if 'num_axial_nodes' in info:
            self.number_axial = int(info['num_axial_nodes'])
        if 'axial_nodes' in info:
            self.axial_nodes = [a.strip(' ') for a in info['axial_nodes'].split(',')]
        if 'number_assemblies' in info:
            self.number_assemblies = int(info['number_assemblies'])

        if 'fixed_problem' == data['optimization']['reproducer']:
            self.fixed_genome = True
        elif 'unique_genes' == data['optimization']['reproducer']:
            self.fixed_genome = True
        
    ## Mcycle Calculation Types
        if self.calculation_type in ["mcycle_loading_pattern_parcs332",
                                     "mcycle_grouped_loading_pattern_parcs332"]:
            if 'ncycles' in info:
                self.ncycles = int(info['ncycles'])
        elif self.calculation_type == "mcycle_inventory_loading_pattern_parcs332":
            if 'ncycles' in info:
                self.ncycles = int(info['ncycles'])
            if 'design_limits' in data['genome']:
                if 'optimize' in data['genome']['design_limits']:
                    if data['genome']['design_limits']['optimize']=='single_cycle':
                        self.ncycles=1
                self.reload_BU = float(data['genome']['design_limits']['reload_burnup'])
                self.num_feedbatch = int(data['genome']['design_limits']['fresh_feed'])
        
        self.settings = data #!TODO: check if this redundancy is actually used anywhere, then remove it.
        del self.input
        return