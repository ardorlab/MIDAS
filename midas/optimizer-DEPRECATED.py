import yaml
from copy import deepcopy
import pickle
## Import MIDAS classes ##
from midas.algorithms import genetic_algorithm as GA
from midas.algorithms import simulated_annealing as SA
from midas.algorithms import random_solutions as RS
from midas.applications.ncsu_lattice import Simulate_Lattice
from midas.applications import ncsu_core
from midas.utils import fitness
from midas.applications import parcs_332
"""
Classes for preparing the optimizer object and calling the necessary optimizer
algorithm and calculation type.

Created by Brian Andersen. 1/7/2019
Updated by Nicholas Rollins. 09/11/2024
"""

def calculate_number_gene_combinations(genome_map):
    """
    Calculates number of possible gene combinations for inputs.

    parameters:
        genome_map: From the input file data dictionary would be the two
        keys: [genome][chromosomes]. 
    
    Written by Brian Andersen. 1/7/2019
    """
    number_changes = 0
    for chromosome in genome_map:
        if chromosome == 'symmetry_list':
            pass
        else:
            number_changes += len(genome_map[chromosome]['map'])  

    return number_changes


class No_Solution_Front(object):
    """
    Class that can be implemented when no solution front
    is desired.
    """
    @staticmethod
    def calculate(solution_list):
        return None


class Optimizer_Factory(object):
    """
    Class for assembling the different optimizations.

    Parameters:
        num_procs: int
            Number of parallel processes that will be used in the optimization
        yaml_file:
            Name of the yaml settings file.
    
    Written by Brian Andersen. 1/7/2019
    Updated by Nicholas Rollins. 09/11/2024
    """
    def __init__(self,num_procs,yaml_lines):
        self.num_procs = num_procs
        ## Inherit input data parameters
        self.input = yaml_lines
        self.file_settings = yaml_lines.file_settings #!temporary redundancy while some input parsing is still being handled elsewhere.
        return

    def build_optimizer(self):
        """
        Assembles all the pieces of the optimization.
        
        Parameters: None
        
        Written by Brian Andersen. 1/7/2019
        """
        methodology = self.input.opt_method
        solution_type_ = self.build_solution() #solution type read from '.yaml' input file: ['optimization']['data_type'] #!still contains file_settings
        population_    = self.build_population()
        if methodology == 'genetic_algorithm':
            generation_    = self.build_generation()
            reproduction_  = self.build_reproduction(generation_.total) #!still contains file_settings
            selection_     = self.build_selection()
            self.optimization = GA.Genetic_Algorithm(solution=solution_type_,
                                                 population=population_,
                                                 generation=generation_,
                                                 reproduction=reproduction_,
                                                 selection=selection_,
                                                 num_procs=self.num_procs,
                                                 inp_obj=self.input,
                                                 file_settings=self.file_settings) #!remove file_settings
        elif methodology == 'simulated_annealing': #!TODO: update with self.input
            generation_    = self.build_generation()
            mutation_      = self.build_mutation(generation_.total) #!still contains file_settings
            cooling_schedule_   = self.build_cooling_schedule()
            fitness_ = self.build_fitness('simulated_annealing')
            self.optimization = SA.Simulated_Annealing(solution=solution_type_,
                                                 population=population_,
                                                 generation=generation_,
                                                 mutation=mutation_,
                                                 num_procs = self.num_procs,
                                                 cooling_schedule= cooling_schedule_,
                                                 fitness=fitness_,
                                                 file_settings=self.file_settings)
        elif methodology == 'reinforcement_learning': #!TODO: update with self.input
            from midas.algorithms import reinforcement_learning as RL
            generation_    = self.build_generation()
            fitness_ = self.build_fitness('genetic_algorithm')
            self.optimization = RL.Reinforcement_Learning(solution=solution_type_,
                                                 population=population_,
                                                 generation=generation_,  
                                                 fitness=fitness_,                     
                                                 num_procs = self.num_procs,
                                                 file_settings=self.file_settings)
        elif methodology == 'random_solutions': #!TODO: update with self.input
            fitness_ = self.build_fitness('genetic_algorithm')
            self.optimization = RS.Random_Solution(solution=solution_type_,
                                                 population=population_,
                                                 fitness=fitness_,
                                                 num_procs= self.num_procs,
                                                 file_settings=self.file_settings)
        elif methodology == 'lava': #!TODO: this optimizer doesn't exist?
            fitness_ = self.build_fitness("lava")
            self.optimization = lava.Volcano(solution=solution_type_,
                                                 fitness=fitness_,
                                                 settings=self.file_settings)
        else:
            raise ValueError("Optimization Type Not Supported")

        return self.optimization

    def build_solution(self):
        """
        Builds the solution type for the optimization.
        
        Parameters: None

        Written by Brian Andersen. 1/7/2019
        """
        data_type_string = self.input.calculation_type

        #!TODO: replace all casmo/simulate calc types with one assignment to ncsu_core
        if data_type_string.lower() == 'pin_lattice':
            solution_type = Simulate_Lattice
            lattice_data = self.file_settings['genome']['additional']
            infile = open("casmo_file_data",'wb') #Information not carried by optimization 
            pickle.dump(lattice_data,infile)        #solutions that is needed to fully write a 
            infile.close()                        #NCSU core simulator input file.
        elif data_type_string.lower() == 'single_assembly_simulate':
            solution_type = ncsu_core.Simulate_Assembly_Solution
            simulate_data = self.file_settings['genome']['additional']

            infile = open("simulate_data",'wb') #Information not carried by optimization 
            pickle.dump(simulate_data,infile)      #solutions that is needed to fully write a 
            infile.close()                      #NCSU core simulator input file.
        elif data_type_string.lower() == "loading_pattern":
            solution_type = ncsu_core.Simulate_Loading_Pattern_Solution
        elif data_type_string.lower() == "loading_pattern_cnn":
            solution_type = ncsu_core.CNN_Loading_Pattern_Solution
            global KERAS 
            KERAS = True
        elif data_type_string.lower() == "fixed_loading_pattern":
            solution_type = ncsu_core.Unique_Assembly_Loading_Pattern_Solution
        elif self.input.ext_model == "PARCS332":
            solution_type = parcs_332.PARCS_Interface #assigned but NOT called (no object created)
        
        genome_key = self.input.chromosomes

        #! Is there any good reason why this needs to be a binary? Why use pickle here?
        infile = open("genome_key",'wb') #The genome key is how the genome used by the
        pickle.dump(genome_key,infile)   #optimization algorithm is decoded into the 
        infile.close()                   #true solution.

        return solution_type

    def build_population(self):
        """
        Builds the population class for the optimization. Intended for the 
        genetic algorithm. Not sure if used in other optimization methodologies

        Parameters: None

        Written by Brian Andersen. 1/7/2019
        """
        population_ = GA.Population()
        population_setting = self.input.population_size
        if population_setting == 'calculate_from_genes':
            number_genes = calculate_number_gene_combinations(self.input.chromosomes)
            population_.calculate_size(number_genes)
        else:
            population_.size = population_setting

        return population_

    def build_generation(self):
        """
        Builds the generation class for the optimization. Intended for genetic
        algorithm, but may be rewritten to incorporate other optimization types.
        There's no reason to have multiple length of optimization classes.

        Parameters: None

        Written by Brian Andersen. 1/7/2019
        """
        generation_ = GA.Generation()
        if self.input.num_generations == 'calculate_from_genes':
            number_genes = calculate_number_gene_combinations(self.input.chromosomes)
            generation_.calculate_total_generations(number_genes)
        else:
            generation_.total = self.input.num_generations

        return generation_
    
    def build_cooling_schedule(self):
        """
        Builds the cooling schedule for the simulated annealing algorithm.
        
        I am unsure how to trigger more than one type of cooling function. This
        is not my biggest priority at this stage of development but will be 
        incorporated in future revisions.
        
        Written by Johnny Klemes. 3/24/2020

        Update Log:
        5/13/2020: Brian Andersen.
            Replaced Johnny Klemes piecewise cooling schedule with the simple exponential decreasing 
            cooling schedule. At time of pull, his code was broken, and there wasn't sufficient documentation
            on how his code was derived to warrant fixing it. 
        """
        temp = self.input.SA_temperature
        alpha = self.input.SA_alpha
        cooling_schedule_= SA.Exponential_Decreasing_Cooling_Schedule(alpha,temp)

        return cooling_schedule_

    def build_mutation(self,number_generations): #!TODO: finish moving this parsing to input_parser
        """
        Builds the mutation class for the optimization. My thought is that new
        optimization methods, such as simulate annealing or tabu search, will use
        a mutation method for generating new solutions. They may be rewritten 
        later to better reflect this. But currently they are still named for 
        mutation methods.

        Available mutation cards:
            mutate_by_type: Only mutates genes into other common genes.
            mutate_by_genome: Randomly mutates based on the provided gene map. Basically the 
                standard mutation method.
            mutate_by_common: Genes are only mutated into genes that already exist in the solution.
                Mostly intended for use with pin optimization. 
            mutate_fixed: Mutations are only performed as translations of genes in the genome.
                This preserves the overall genome, and should be used if every genome is unique.
                This encompasses problems such as the fuel loading problem, where the burnt fuel
                inventory must be preserved.
            double_mutation (specified by list of two mutation methods): Allows the optimization
                to be performed using the two specified mutation methods. 
        
        Parameters: None

        Written by Brian Andersen. 1/7/2019
        """
        MUTE_SETTING = self.input.mutation_settings

        temp1 = MUTE_SETTING['initial_rate'] #Initial mutation rate
        temp2 = MUTE_SETTING['final_rate']   #Final mutation rate
        temp3 = 1                            #Number of mutations per instance
        
        #This sort of works, from what I remember, but I wouldn't mess with it
        #right now. 
        if 'multiple_genes_per_chromosome' in self.file_settings['genome']:
            foo = self.file_settings['genome']['multiple_genes_per_chromosome']
            if foo:
                temp4 = deepcopy(MUTE_SETTING['dictionary'])
                for key in MUTE_SETTING['dictionary']:
                    if 'mutate_by_type' == MUTE_SETTING['dictionary'][key]:
                        if MUTE_SETTING['common_chromosomes'] == "determine_from_maps":
                            temp4['common_gene_list'] = None
                        else:
                            temp4['common_gene_list'] = MUTE_SETTING['common_chromosomes']
                    elif 'mutate_by_genome' == MUTE_SETTING['dictionary'][key]:
                        map_ = self.input.chromosomes
                        temp4['map_genome'] = map_
                    elif 'mutate_by_common' == MUTE_SETTING['dictionary'][key]:
                        map_ = self.input.chromosomes
                        temp4['map_genome'] = map_
                mutator_ = GA.Submutation_Mutator(temp1,temp2,temp3,temp4)
        
        #Lists are used to specify multiple types of mutation. Currently only two different
        #mutation types may be utilized at once. This is hard coded. If I were better at
        #programming it would be easier to make multiple types easily supported.
        elif type(self.file_settings['optimization']['mutation']['method']) == list:
            type0 = self.file_settings['optimization']['mutation']['method'][0]
            type1 = self.file_settings['optimization']['mutation']['method'][1]
            if type0 == 'mutate_by_type':
                type0 = GA.Mutate_By_Type
            elif type0 == 'mutate_by_genome':
                type0 = GA.Mutate_By_Genome
            elif type0 == 'mutate_by_common':
                type0 = GA.Mutate_By_Common
            elif type0 == 'mutate_fixed':
                type0 = GA.Fixed_Genome_Mutator
            if type1 == 'mutate_by_type':
                type1 = GA.Mutate_By_Type
            elif type1 == 'mutate_by_genome':
                type1 = GA.Mutate_By_Genome
            elif type1 == 'mutate_by_common':
                type1 = GA.Mutate_By_Common
            elif type1 == 'mutate_fixed':
                type1 = GA.Fixed_Genome_Mutator
            mutator_ = GA.Double_Mutator(temp1,temp2,temp3,type0,type1,self.file_settings)

        elif MUTE_SETTING['method'].lower() == 'mutate_by_type':
            mutator_ = GA.Mutate_By_Type(temp1,temp2,temp3,self.file_settings)
        elif MUTE_SETTING['method'].lower() == 'mutate_by_genome':
            mutator_ = GA.Mutate_By_Genome(temp1,temp2,temp3,self.file_settings)
        elif MUTE_SETTING['method'].lower() == 'mutate_by_unique':
            mutator_ = GA.Mutate_By_Unique(temp1,temp2,temp3,self.file_settings)
        elif MUTE_SETTING['method'].lower() == 'mutate_by_common':
            print("Original file settings")  
            print(self.file_settings)
            mutator_ = GA.Mutate_By_Common(temp1,temp2,temp3,settings=self.file_settings)

        elif MUTE_SETTING['method'].lower() == 'mutate_fixed':
            mutator_ = GA.Fixed_Genome_Mutator(temp1,temp2,temp3,self.file_settings)

        else:
            raise ValueError("Invalid Mutation Type")
        if mutator_.constant_mutation == True:
            pass
        else:
            mutator_.calculate_rate_increase(number_generations)
        
        return mutator_

    def build_reproduction(self,number_generations):
        """
        Builds the reproduction class for the optimization. Intended for the
        Genetic Algorithm optmization methodology, but may have uses elsewhere.
        Don't really know. Might have to be rethought somewhat to accomodate
        other optimization methodologies.

        Parameters: 
            number_generations: int
                The number of generations over which the optimization occurs.

        There are three types of reproduction specified within the code, but currently
        only two of the three methods work well. The third method seems that it has some
        bugs in it, but I have been having trouble isolating and reproducing the bugs.

        List of implemented reproduction types. #!TODO: this list appears to be outdated.

        multiple_genes_per_chromosome (The one that has unknown bug): The first way I thought to
            accomodate for burned fuel assemblies in the optimization. Burned Fuel Assemblies have
            a sub-genome. But was found to violate conservation of used fuel assemblies.
        fixed_genome_problem: The second way I came up with to preserve burned fuel assemblies.
            Every single genome in the optimization is unique and may only be expressed once. 
        reproduction: Automatically loaded if the other two reproduction methods are not loaded.

        Written by Brian Andersen. 1/7/2019
        """
        mute = self.build_mutation(number_generations)
        method = self.input.reproducer
        if 'multiple_genes_per_chromosome' == method:
            reproducer_ = GA.Dictionary_Genome_Reproducer(mute)
        elif 'fixed_problem' == method:
            reproducer_ = GA.Fixed_Gene_Reproducer(mute)    
        elif 'unique_genes' == method: 
            reproducer_ = GA.Unique_Genome_Reproducer(mute,self.file_settings)
        elif 'mcycle' == method: 
            reproducer_ = GA.Mcycle_Reproducer(mute,self.file_settings)
        else:
            reproducer_ = GA.Reproduction(mute)

        return reproducer_

    def build_fitness(self,method):
        """
        Builds the fitness class for the optimization. All fitness functions are 
        currently written for the genetic algorithm. This and the fitness functions
        will need to be generalized for other optimization methodologies.

        Parameters: None

        Written by Brian Andersen. 1/7/2019
        Updated by Nicholas Rollins. 09/11/2024
        """
        #Fitness card is allowed to be specified either directly under the optimization key,
        #or under the selection key. Under the selection key makes sense for genetic algorithms
        #but it doesn't necessarily make sense for methods such as simulated annealing. So it
        #is allowed in both locations for now. 
        select_line = self.input.selection_settings
        fit_method = self.input.fitness_method

        if fit_method  == 'ranked':
            fmethod_ = fitness.Ranked_Fitness()
        elif fit_method == 'binned':
            fmethod_ = fitness.Binned_Fitness()
        elif fit_method.upper() == 'ACDPF':
            fmethod_ = fitness.ACDPF()
        elif fit_method == 'weighted_positive':
            if method =='genetic_algorithm':
                return fitness.Genetic_Algorithm_Weighted_Positive()
            else:
                return fitness.Fitness()
        elif fit_method == 'adaptive':
            if method =='genetic_algorithm':
                return fitness.Genetic_Algorithm_Adaptive()
            else:
                return fitness.Fitness()
        elif fit_method  == 'weighted':
            if method =='genetic_algorithm':
                return fitness.Genetic_Algorithm_Weighted()
            else:
                return fitness.Fitness()
        elif fit_method == 'quantum':
            return fitness.Quantum_Fitness()
        #The Fitness functions are designed so that minimizing the score is the objective.
        #Genetic algorithms work the opposite way. The fitness maximizer corrects this 
        #situation.
        if method == 'genetic_algorithm':
            fitness_ = fitness.Fitness_Maximizer(fmethod_)
            return fitness_
        else:
            return fmethod_

    def build_selection(self):
        """
        Builds the selection class for the optimization. May need to be rewritten
        for other optimization methodologies.
        
        Parameters: None

        Written by Brian Andersen. 1/7/2019
        """
        if type(self.input.selection_settings) == dict:
            select_line = self.input.selection_settings
            fitness_ = self.build_fitness('genetic_algorithm')
            if 'method' in select_line: 
                method_type = select_line['method'].lower()
                if method_type == 'tournament':
                    method_ = GA.GA_Selection.tournament
                elif method_type == 'roulette':
                    method_ = GA.GA_Selection.roulette
                else: 
                    raise NotImplementedError
            else:
                error = "No method for selection is specified. The two available methods are"
                error += " roulette and tournament."
                raise KeyError(error)
            solution_front = No_Solution_Front()
            selection_ = GA.GA_Selection(fitness_,method_,solution_front)
        elif self.input.selection_settings.upper() == 'MOOGLE':        
            selection_ = GA.MOOGLE()
        else:
            raise NotImplementedError
        
        return selection_
