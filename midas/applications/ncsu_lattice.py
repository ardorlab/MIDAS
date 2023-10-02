import os
import sys
import pickle
import yaml
from multiprocessing import Pool
from midas.utils.solution_types import Solution, compare_substrings, return_triangular_string
from midas.utils.metrics import Plotter

class Simulate_Lattice(Solution):
    """
    Class specific to NCSU lattice simulator functions. Contains the information to generate an 
    initial solution and evaluate it.
    
    Parameters: None

    Written by Brian Andersen. 1/7/2020
    """ 
    def __init__(self):
        Solution.__init__(self)

        self.plot = Plotter.lattice
        self.boron = None
        self.power = None
        self.title = None
        self.pin_list = None
        self.depletion = None
        self.fuel_temperature = None
        self.control_rod_line = None
        self.spacer_grid_line = None
        self.moderator_temperature = None
        self.reactor = None
        self.void = None
        self.s3c = False

    def add_additional_information(self,settings):
        """
        Adds information on assembly type, the number of pins in the assembly,
        and the model type to the class instance.

        Parameters: 
            settings: Dictionary
                Contains the optimization settings. Used to fill in the additional
                information.

        Written by Brian Andersen. 1/7/2020
        """
        info = settings['genome']['additional']
        if 'boron' in info:
            self.boron = info['boron']
        if 'power' in info:
            self.power = info['power']
        if 'title' in info:
            self.model = info['title']
        if 'pin_list' in info:
            self.pin_list = info['pin_list']
        if 'depletion' in info:
            self.depletion = info['depletion']
        if 'spacer_grid' in info:
            self.spacer_grid_line = info['spacer_grid']
        if 'fuel_temperature' in info:
            self.fuel_temperature = info['fuel_temperature']
        if 'moderator_temperature' in info:
            self.moderator_temperature = info['moderator_temperature']
        if 'reactor' in info:
            self.reactor = info['reactor']
        if 'void' in info:
            self.void = info['void']
        if 's3c' in info:
            self.s3c = info['s3c']
        if 'control_rod' in info:
            self.control_rod_line = info['control_rod']

    def evaluate(self):
        """
        Evaluates the given NCSU lattice solution.

        Parameters: None

        Written by Brian Andersen. 1/7/2020
        """
        infile = open("genome_key",'rb')
        genome_key = pickle.load(infile)
        infile.close()

        fuel_list = []
        pin_list = []           #Fuel pin information is designated with the
        fuel_card_string = ''   #genome in the input file, hence why this info
        for gene in genome_key: #is read here.
            fuel_card_string += genome_key[gene]['card'] + '\n'

        for gene in self.genome:                   #Decode the solution
            fuel_list.append(genome_key[gene]['fuel']) #genome into the info
            pin_list.append(genome_key[gene]['pin'])   #needed by NCSU lattice

        fuel_str = return_triangular_string(fuel_list)
        pin_str  = return_triangular_string(pin_list)

        self.write_input_file(fuel_card_string,fuel_str,pin_str)
        outline = "cd {} ;".format(self.name) 
        outline += "casmo4e -v u1.00.02 {}.inp".format(self.name) 
        os.system(outline)

        file_ = open("{}/{}.out".format(self.name,self.name),'r')
        file_lines = file_.readlines()
        file_.close()
        kinf_list = None

        if 'max_kinf' in self.parameters:
            if 'max_kinf' in self.neural_network:
                pass
            else:
                kinf_list = Extractor.kinf(file_lines)
                self.parameters['max_kinf']['value'] = max(kinf_list)
        #Because peak versus max kinf trips me up enough, they are both now
        #valid keys to use in calculating the highest kinf value in the NCSU lattice 
        #analysis, but only one or the other may be used.
        elif 'peak_kinf' in self.parameters:
            if 'peak_kinf' in self.neural_network:
                pass
            else:
                kinf_list = Extractor.kinf(file_lines)
                self.parameters['peak_kinf']['value'] = max(kinf_list)
        if 'eoc_kinf' in self.parameters:
            if 'eoc_kinf' in self.neural_network:
                pass
            else:
                if not kinf_list:
                    kinf_list = Extractor.kinf(file_lines)
                self.parameters['eoc_kinf']['value'] = kinf_list[-1]
        if 'peak_pin_power' in self.parameters:
            if 'peak_pin_power' in self.neural_network:
                pass
            else:
                max_power_list = Extractor.max_pin_powers(file_lines)
                self.parameters['peak_pin_power']['value'] = max(max_power_list)
        if 'boc_kinf' in self.parameters:
            if not kinf_list:
                kinf_list = Extractor.kinf(file_lines)
            self.paramaters['boc_kinf']['value'] = kinf_list[0]
        if 'enrichment' in self.parameters:
            if 'enrichment' in self.neural_network:
                pass
            else:
                self.parameters['enrichment']['value'] = Extractor.enrichment(file_lines)
        if 'input_enrichment' in self.parameters:
            infile = open("genome_key",'rb')
            genome_key = pickle.load(infile)
            infile.close()

            enrichment_sum = 0.
            for gene in self.genome:
                enrichment_sum += genome_key[gene]['enrichment']
            self.parameters['input_enrichment']['value'] = enrichment_sum

    def write_input_file(self,fuel_card_string,fuel_string,pin_string):
        """
        Writes all of the information stored in the NCSU lattice class instance 
        into a valid NCSU lattice Input File.

        Parameters:
            fuel_card_string: str
                The NCSU lattice fuel cards, in triangular format.
            pin_card_string: str
                The NCSU lattice pin cards, in triangular format.

        Written by Brian Andersen 1/8/2020
        """
        os.system("mkdir {}".format(self.name))
        lattice_file = open("{}/{}.inp".format(self.name,self.name),'w')
        if not self.title:                            #Write Title to input  
            lattice_file.write("TTL * NCSU lattice Analysis\n")      #File
        else:
            lattice_file.write("TTL {}\n".format(self.title))

        lattice_file.write(fuel_card_string)           #Fuel Types

        for card in self.pin_list:                   #Pin Types
            lattice_file.write(card + '\n')

        lattice_file.write('LFU\n')
        lattice_file.write(fuel_string)                #Write Fuel type Map

        lattice_file.write('LPI\n')
        lattice_file.write(pin_string)                 #Write Pin Type Map

        lattice_file.write("PDE {} \n".format(self.power))   #Assembly Power

        if not self.fuel_temperature: #Assembly
            pass                          #Temperature
        else:
            lattice_file.write("TFU = {}\n".format(self.fuel_temperature))  
                                                              
        if not self.moderator_temperature: #write
            pass                               #coolant temperature
        else:
            lattice_file.write("TMO = {}\n".format(self.moderator_temperature)) 

        if not self.boron:
            pass
        else:
            lattice_file.write("BOR = {}\n".format(self.boron))

        if not self.void:
            pass
        else:
            lattice_file.write("VOI = {}\n".format(self.void))

        if not self.reactor:
            pass
        else:
            lattice_file.write("{}\n".format(self.reactor)) #initialize reactor type

        if not self.spacer_grid_line:
            pass
        else:                      #Add spacer grid 
            lattice_file.write("{}\n".format(self.spacer_grid_line)) #info to NCSU lattice

        if not self.control_rod_line:  
            pass
        else:                    #Add spacer grid
            lattice_file.write("{}\n".format(self.control_rod_line)) #info to NCSU lattice

        if not self.s3c:    
            pass
        else:                  #Add S3C card to NCSU lattice
            lattice_file.write("s3c \n") #file

        if not self.depletion:
            pass
        else:                      #Adds specific 
            lattice_file.write("DEP {}\n".format(self.depletion)) #depletion to NCSU lattice

        lattice_file.write("STA"+"\n")
        lattice_file.write("END"+"\n")
        lattice_file.close()


class Extractor(object):
    """
    Class for extracting the values from the output files of NCSU lattice files. 

    Parameters: None

    Written by Brian Andersen. 1/7/2020
    """
    @staticmethod
    def max_pin_powers(file_lines):
        """
        Accepts NCSU lattice file lines as the function argument and returns a list of
        the maximum pin powers over the course of the NCSU lattice depletion.

        Parameters:
            file_lines: str_list
                The file lines of the output file being examined.

        WRitten by Brian Andersen. 1/7/2020
        """
        max_power_list = []
        search_line = "* POWER DISTRIBUTION    PEAK: LL ="
        desired_value = len(search_line.strip().split())
        for line in file_lines:
            elems = line.strip().split()
            if not elems:
                pass
            else:
                if search_line in line:
                    max_power_list.append(float(elems[desired_value]))
        
        return max_power_list

    @staticmethod
    def burnup(file_lines):
        """
        Accepts NCSU lattice file lines as the function argument and returns the
        max kinf value over the course of the depletion
        
        Parameters:
            file_lines: str_list
                The file lines of the output file being examined.

        WRitten by Brian Andersen. 1/7/2020
        """
        burnup_list = []
        search_line = "BURNUP =    MWD/KG   K-INF =     M2 =     B2 = "
        length_string = " BURNUP = "
        desired_value = len(length_string.strip().split())
        for line in file_lines:
            elems = line.strip().split()
            if len(elems) > 0:
                if compare_substrings(search_line,line):
                    burnup_list.append(float(elems[desired_value]))
        return burnup_list   

    @staticmethod
    def kinf(file_lines):
        """
        Accepts NCSU lattice file lines as the function argument and returns the
        max kinf value over the course of the depletion
        
        Parameters:
            file_lines: str_list
                The file lines of the output file being examined.

        WRitten by Brian Andersen. 1/7/2020
        """
        kinf_list = []
        search_line = "BURNUP =    MWD/KG   K-INF =     M2 =     B2 = "
        length_string = " BURNUP =    0.000 MWD/KG   K-INF ="
        desired_value = len(length_string.strip().split())
        for line in file_lines:
            elems = line.strip().split()
            if len(elems) > 0:
                if compare_substrings(search_line,line):
                    kinf_list.append(float(elems[desired_value]))
        return kinf_list   

    @staticmethod
    def enrichment(file_lines):
        """
        Returns the enrichment of the NCSU lattice lattice based on the enrichment
        of the first depletion state.
        
        Parameters:
            file_lines: str_list
                The file lines of the output file being examined.

        WRitten by Brian Andersen. 1/7/2020
        """
        line_count = 1000
        match_found = False
        for line in file_lines:
            if "WEIGHT PERCENT OF  U-235 X 1.E+00" in line:
                match_found = True
                line_count = 0
                string_enrichment_list = []
            elif "MAX VALUE:" in line:
                match_found = False
                break
            if match_found == True and line_count > 0:
                elems = line.strip().split()
                string_enrichment_list.extend(elems[line_count:])
            
            line_count +=1

        float_enrichment_list = [float(enrichment.replace("*","")) for enrichment in string_enrichment_list]

        return sum(float_enrichment_list)



