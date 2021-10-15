import os
import sys
import random
import copy
from metrics import Lava_Metric_Toolbox

class Volcano(object):
    """
    The main function for the lava optimization.
    """
    def __init__(self,solution,fitness,settings):
        self.solution = solution #The solution type solved by the optimization.
        self.file_settings = settings #The general settings provided by the yaml input file.
        self.lava_dictionary = {} #Tracks the different lava flows flowing from the volcano.
        self.number_flows = settings['optimization']['number_flows']
        self.fitness = fitness
 
    def initialize_solution(self,name):
        """
        Generates a blank solution to the optimization problem.

        This might get moved to the lava optimization. It is hard to say
        which is the best move to make.
        """
        foo = self.solution()
        foo.name = name
        foo.parameters = copy.deepcopy(self.file_settings['optimization']['objectives'])
        foo.add_additional_information(self.file_settings)
        
        return foo
    def main_in_serial(self):
        """
        Performs the optimization in serial.
        """
        for i in range(self.number_flows):
            foo = self.initialize_solution(f"flow_{i}_initial") #Generate the initial solution for the flow
            if foo.fixed_genome:                                #If the problem is fixed
                foo.generate_initial_fixed(self.file_settings['genome']['chromosomes'],         #Generate fixed optimization problem.
                                           self.file_settings['optimization']['fixed_groups'])  
                foo.evaluate()                        #Since this is done in serial evaluating solutions now.
                self.fitness.calculate([foo])
                flow = Fixed_Lava_Flow(foo,self.file_settings) #Generate the fixed flow class
                self.lava_dictionary[f"Flow_{i}"] = flow       #Create the dictionary entry for the flow       
            else:
                foo.generate_initial(self.file_settings['genome']['chromosomes'])
                foo.evaluate()
                self.fitness.calculate([foo])
                flow = Unfixed_Lava_Flow(foo,self.file_settings)
                self.lava_dictionary[f"Flow_{i}"] = flow
            
        converged = False    #optimization will continue until all solutions have converged.
        while not converged:
            for key_flow in self.lava_dictionary:           
                current_flow = self.lava_dictionary[key_flow] #To keep lines from getting too long
                if not current_flow.flow_convergence:
                    current_flow.branch()  #Create the appropriate branching flows from the current flow
                    for i,flow in enumerate(current_flow.new_flow_list):
                        solution = self.initialize_solution(f'{key_flow}_{current_flow.generation}_{current_flow.solution_count}')
                        current_flow.solution_count += 1
                        solution.genome = flow
                        solution.evaluate()
                        self.fitness.calculate([solution])
                        current_flow.compare_to_best(solution) #Determines if new solution replaces previous solution.
                    if not current_flow.change_in_flow:
                        current_flow.flows_without_change += 1
                        if current_flow.flows_without_change > 2*len(current_flow.current_best_solution.genome):
                            current_flow.flow_convergence = True
                current_flow.generation += 1
                current_flow.solution_count = 0

class Unfixed_Lava_Flow(object):
    """
    The class where the actual individual optimization is performed. Multiple Lava flows may be combined into a single Volcano optimization.
    """
    def __init__(self,base_solution,settings):
        self.flow_number = None
        self.current_best_solution = base_solution #The current best solution in this lava flow.
        self.new_flow_list = []   #List of solutions created by branching from the current best solution.
        self.flow_convergence = False   #Determines if the solution has converged
        self.flow_map = {}
        self.plotter = Lava_Metric_Toolbox()
        self.generation = 0
        self._fill_flow_map(settings)
        self.solution_count = 0
        self.flows_without_change = 0
        self._moving_foward = random.choice([True,False])    #Records which direction the optimization is moving in for advancing the current position.
        self.current_position = random.randint(0,len(self.current_best_solution.genome)-1)    #The current position at which the solution is being branched. 

    def _fill_flow_map(self,settings):
        """
        Fills in the flow map for the lava flow. This is identical to the genome map
        in genetic algorithm mutation.
        """
        if not settings:
            pass
        else:
            chrom_settings = settings['genome']['chromosomes']
            for chrom in chrom_settings:
                if chrom == "symmetry_list":
                    self.symmetry_list = chrom_settings[chrom]
                else:
                    self.flow_map[chrom] = copy.deepcopy(chrom_settings[chrom]['map'])

    def branch(self):
        """
        Creates the new solution branches for the optimization.
        
        In the unfixed lava flow variation, branching constitutes replacing the
        current decision variable with all valid decision variables for that
        position in the solution problem.
        """
        self.new_flow_list = []
        for key in self.flow_map: #For the number of possible branches
            if self.current_best_solution.genome[self.current_best_solution] == key: #No sense in duplicating the current best flow. 
                pass
            elif self.flow_map[key][self.current_position]:                    #Valid decision variables are represented using a binary key, with 1 meaning yes. 
                new_flow = copy.deepcopy(self.current_best_solution.genome)  #1 will validate to true and 0 goes to false. So if if statement activates the current
                new_flow[self.current_position] = key                        #best genome is copied, then the new genome is inserted into the current position. 
                self.new_flow_list.append(new_flow)

        self._advance_gene_position()

    def _advance_gene_position(self):
        """
        Advances the gene position of the flow. 
        """
        if self._moving_foward:
            self.current_position += 1
            if self.current_position == len(self.current_best_solution.genome): #At the point where reversal should take place, the current position should equal 
                self._moving_foward = False  #the length so the direction is reversed. Then the position is moved two positions backward,       
                self.current_position -= 2   # since moving only one position backward should duplicate what was just done.
        else:
            self.current_position -= 1
            if self.current_position == -1:
                self._moving_foward = True
                self.current_position += 2

    def compare_to_best(self,comparison):
        """
        Compares the current best solution to the new solution provided. If the new solution has a better objective value than the current best value, the new solutoin
        becomes the current best solution.

        Optionally allows for temperature functions to be applied as well in the future, should that be deemed necessary.
        """
        if comparison.fitness < self.current_best_solution.fitness:
            self.current_best_solution = copy.deepcopy(comparison)
            self.change_in_flow = True 
            self.flows_without_change = 0
            self.plotter.update_best_values(self.flow_number,
                                            self.generation,
                                            self.solution_count,
                                            self.current_best_solution)

class Fixed_Lava_Flow(Unfixed_Lava_Flow):
    """
    Lava flow for fixed optimization problems. This optimization type requires the branching to be performed in a different manner, as the unfixed
    method lets all viable decision variables be placed into the problem. 
    """
    def __init__(self,base_solution,settings):
        Unfixed_Lava_Flow.__init__(self,base_solution,settings)
        self._fill_common_decision_list(settings)

    def _fill_common_decision_list(self, settings):
        """
        Fills in the parameter common_gene_list. If the common gene list is provided it
        directly takes this value. Otherwise, it determines the common gene list from
        the provided gene maps of the chromosomes.

        Parameters:
            settings: Dictionary
                The input yaml file, as a dictionary.

        WRitten by Brian Andersen. 1/21/2020.
        """
        MUTE_SETTING = settings['optimization']['mutation']['common_chromosomes']
        if MUTE_SETTING == 'determine_from_maps':
            chromosomes = settings['genome']['chromosomes']
            examine_dictionary = {}
            current_map_count = 0
            for gene in chromosomes:
                if not examine_dictionary:
                    examine_dictionary[current_map_count] = {}
                    examine_dictionary[current_map_count]['map'] = chromosomes[gene]['map']
                    examine_dictionary[current_map_count]['gene_list'] = [gene]
                    current_map_count += 1
                else:
                    no_match_found = True
                    for i in examine_dictionary:
                        if examine_dictionary[i]['map'] == chromosomes[gene]['map']:
                            examine_dictionary[i]['gene_list'].append(gene)
                            no_match_found = False
                            break
                    if no_match_found:
                        examine_dictionary[current_map_count] = {}
                        examine_dictionary[current_map_count]['map'] = chromosomes[gene]['map']
                        examine_dictionary[current_map_count]['gene_list'] = [gene]
                        current_map_count += 1

            common_list = []
            for i in examine_dictionary:
                if len(examine_dictionary[i]['gene_list']) > 1:
                    common_list[i] = examine_dictionary[i]['gene_list']
            self.common_decision_list = common_list
        elif MUTE_SETTING == 'determine_from_fixity':
            chromosomes = settings['genome']['chromosomes']
            self.common_decision_list = {}
            for gene in chromosomes:
                group = chromosomes['gene']['gene_group']
                if group in self.common_decision_list:
                    self.common_decision_list[group].append(gene)
                else:
                    self.common_decision_list[group] = [gene]
        else:
            self.common_decision_list = copy.deepcopy(MUTE_SETTING)

    def branch(self):
        """
        Branching for fixed optimization problems takes place in two stages. IN the first stage, common decision variables, 
        such as assemblies with and without bp but of similar enrichment, are branched. This is basically the same as the Unfixed Optimization.
        In the second stage, all decision variables ahead of the current position are swapped with the decision variable in the current position.
        Only swapping the directions ahead is because for decisions behind this solution should already have been analyzed. This is at least true
        if the solution is not advancing. This also means that some sections will involve a large amount of work, and other sections will be relatively
        simple.
        """
        self.new_flow_list = []
        for key in self.common_decision_list:
            if self.current_best_solution.genome[self.current_position] in self.common_decision_list[key]:
                for decisions in self.common_decision_list[key]:
                    if decisions == self.current_best_solution.genome[self.current_position]:
                        pass
                    else:
                        new_flow = copy.deepcopy(self.current_best_solution.genome)
                        new_flow[self.current_position] = decisions
                        self.new_flow_list.append(new_flow)

        if self._moving_foward:
            for i,decision in enumerate(self.current_best_solution.genome[(self.current_position+1):]): #For all decions from current decision +1 to end of problem
                if decision == self.current_best_solution.genome[self.current_position]:
                    pass
                else:
                    if self.flow_map[decision][self.current_position]:
                        if self.flow_map[self.current_best_solution.genome[self.current_position]][i+self.current_position]:
                            new_flow = copy.deepcopy(self.current_best_solution.genome)
                            new_flow[self.current_position] = decision
                            new_flow[self.current_position+i] = self.current_best_solution.genome[self.current_position]
                            self.new_flow_list.append(new_flow)
        else:
            for i,decision in enumerate(self.current_best_solution.genome[:self.current_position]): #for all decisions from beginning to current decision.
                if decision == self.current_best_solution.genome[self.current_position]:
                    pass
                else:
                    if self.flow_map[decision][self.current_position]:
                        if self.flow_map[self.current_best_solution.genome[self.current_position]][i]:
                            new_flow = copy.deepcopy(self.current_best_solution.genome)
                            new_flow[self.current_position] = decision
                            new_flow[i] = self.current_best_solution.genome[self.current_position]
                            self.new_flow_list.append(new_flow)

        self._advance_gene_position()
