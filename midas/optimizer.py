## Import Block ##
import os
import logging
from pathlib import Path
from shutil import rmtree
from copy import deepcopy
from multiprocessing import Pool
from itertools import repeat
import csv
import pickle

from midas.algorithms import genetic_algorithm as GA
from midas.algorithms import bayesian_optimization as BO
from midas.algorithms import simulated_annealing as SA
from midas.algorithms import parallel_simulated_annealing as PSA
from midas.utils import optimizer_tools as optools
from midas.codes import parcs342, parcs343
from midas.utils import LWR_fuelcyclecost
from midas.utils import LWR_averageenrichment
from midas.codes import nuscale_lut
from midas.utils import termination_criteria as TC


## Classes ##
class Optimizer():
    """
    #!TODO: write docstring.
    
    Written by Nicholas Rollins. 09/23/2024
    """
    def __init__(self, inp_lines):
        self.input = inp_lines
        self.termination_criteria = TC.Termination_Criteria()
    
    def build_optimizer(self):
        """
        Assembles all the pieces of the optimization and create the
        algorithm object.
        
        Written by Brian Andersen. 01/07/2019
        SA added by Johnny Klemes. 03/22/2020
        RL added by G. K. Delipe.  03/24/2023
        Updated by Nicholas Rollins. 09/11/2024
        BO added by Cole Howard. 10/21/2024
        SA updated by Jake Mikouchi 04/25/2025
        """
        methodology = self.input.methodology
        num_gene_combos = self.calculate_number_gene_combinations(self.input.genome)
        
        self.population = optools.Population(self.input.population_size, num_gene_combos)
        self.generation = optools.Generation(self.input.num_generations, num_gene_combos)
        self.fitness    = optools.Fitness()
        self.eval_func  = None
        if self.input.code_interface == "parcs342":
            self.eval_func = parcs342.evaluate #assign, don't execute.
        elif self.input.code_interface == "parcs343":
            self.eval_func = parcs343.evaluate
        elif self.input.code_interface == "nuscale_database":
            self.eval_func = nuscale_lut.evaluate
        
        if methodology == 'genetic_algorithm':
            self.algorithm = GA.Genetic_Algorithm(self.input)
        elif methodology == 'bayesian_optimization':
            self.algorithm = BO.Bayesian_Optimization(self.input)
        elif methodology == 'simulated_annealing' and self.input.num_procs == 1:
            self.algorithm = SA.Simulated_Annealing(self.input)
        elif methodology == 'simulated_annealing' and self.input.num_procs > 1:
            self.algorithm = PSA.Parallel_Simulated_Annealing(self.input)
        #!TODO: Add the other algorithms back in.

        return
    
    def calculate_number_gene_combinations(self, genome):
        """
        Calculates number of possible gene combinations for inputs.

        parameters:
            genome_map: From the input file data dictionary under the keys: 
            [decision_variables][parameters]. 
        
        Written by Brian Andersen. 1/7/2019
        Updated by Nicholas Rollins. 09/24/2024
        """
        number_changes = 0
        for gene in genome:
            number_changes += len(genome[gene]['map'])  
        return number_changes
    
    def generate_solution(self,name,chromosome=None):
        """
        Function for creating new solutions/indiviuals/particles in the optimization space.
        
        Written by Nicholas Rollins. 09/27/2024
        """
        #generate blank solution object
        soln = optools.Solution(f"{name}")
        
        #copy objectives and constraints from input
        soln.parameters = deepcopy(self.input.objectives)
        for key in soln.parameters.keys():
            soln.parameters[key]['value'] = None #placeholder to be filled by objective function.
        
        LWR_core_parameters = [self.input.nrow, self.input.ncol, self.input.num_assemblies, self.input.symmetry]
        if chromosome:
            soln.chromosome = chromosome
        else: #generate random chromosome
            soln.chromosome = soln.generate_initial(self.input.calculation_type, LWR_core_parameters,\
                                                    self.input.genome, self.input.batches) #'batches' is None when not applicable.

        return soln

    def main(self, restart=False):
        """
        Primary execution logic for the Optimizer class. The intention is for all 
        optimization jobs, regardless of algorithm, calculation type, or code model, 
        to use this function.
        
        Written by Nicholas Rollins. 09/27/2024
        """
    ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")
        
        if not (self.input.methodology == 'simulated_annealing' and self.input.num_procs > 1):
            """
            All normal algorithms should go through this execution pathway
            """
            if not restart:
        ## Create results directory and/or clear old
                cwd = Path(os.getcwd())
                results_dir = cwd.joinpath(self.input.results_dir_name)
                if os.path.exists(results_dir):
                    logger.debug("Overwriting existing results directory...")
                    rmtree(results_dir, ignore_errors=True)
                    os.mkdir(results_dir)
                else:
                    os.mkdir(results_dir)
            
        ## Delete old results file, if necessary, to avoid confusion
                # this step is useful since the new results file won't be written until after Generation 0 is completed.
                if os.path.exists("optimizer_results.csv"):
                    logger.debug("Removing existing 'optimizer_results.csv' results file...")
                    os.remove("optimizer_results.csv")
                
        ## Initialize beginning population
                logger.info("Generating initial population of %s individuals...", self.input.population_size)
                self.population.current = []
                for i in range(self.population.size):
                    self.population.current.append(self.generate_solution(f'Gen_0_Indv_{i}'))
                
                pool = Pool(processes=self.input.num_procs) #initialize parallel execution
                
        ## Evaluate fitness
                logger.info("Calculating fitness for generation %s...", self.generation.current)
                ## Execute and parse objective/constraint values
                self.population.current = pool.starmap(self.eval_func, zip(self.population.current, repeat(self.input)))
                if 'cost_fuelcycle' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, self.input)
                if 'av_fuelenrichment' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, self.input)
                ## Calculate fitness from objective/constriant values
                for soln in self.population.current:
                    soln.fitness_value = self.fitness.calculate(soln.parameters)
                logger.info("Done!")
        
        ## Archive initial results
                for soln in self.population.current:
                    self.population.archive['solutions'].append(soln.chromosome)
                    self.population.archive['fitnesses'].append(soln.fitness_value)
                    self.population.archive['parameters'].append(soln.parameters)
                
                ## Only initialize the results file the first time.
                archive_header = ["Generation","Individual","Fitness Value"]
                for param in self.input.objectives.keys():
                    archive_header.append(str(param))
                archive_header.append("Chromosome")
                ## write output file
                with open("optimizer_results.csv", 'w') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
                    csvwriter.writerow(archive_header)
                best_soln_index = [s.fitness_value for s in self.population.current].index(max([s.fitness_value for s in self.population.current]))
                for i in range(len(self.population.current)):
                    soln = self.population.current[i]
                    soln_result_list = [str(self.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                    for param in soln.parameters.keys():
                        if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                            soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                        else:
                            try:
                                soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                            except:
                                soln_result_list.append(str(soln.parameters[param]['value']))
                    for gene in soln.chromosome:
                        soln_result_list.append(str(gene))
                    ## write to output file
                    with open("optimizer_results.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_list)
                    if i == best_soln_index:
                        best_soln_string = ",".join(soln_result_list)
            
        ## Archive initial results
                for soln in self.population.current:
                    self.population.archive['solutions'].append(soln.chromosome)
                    self.population.archive['fitnesses'].append(soln.fitness_value)
                    self.population.archive['parameters'].append(soln.parameters)
                
                ## Only initialize the results file the first time.
                archive_header = ["Generation","Individual","Fitness Value"]
                for param in self.input.objectives.keys():
                    archive_header.append(str(param))
                archive_header.append("Chromosome")
                ## write output file
                with open("optimizer_results.csv", 'w') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
                    csvwriter.writerow(archive_header)
                
                best_soln_index = [s.fitness_value for s in self.population.current].index(max([s.fitness_value for s in self.population.current]))
                for i in range(len(self.population.current)):
                    soln = self.population.current[i]
                    soln_result_list = [str(self.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                    for param in soln.parameters.keys():
                        if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                            soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                        else:
                            try:
                                soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                            except:
                                soln_result_list.append(str(soln.parameters[param]['value']))
                    for gene in soln.chromosome:
                        soln_result_list.append(str(gene))
                    ## write to output file
                    with open("optimizer_results.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_list)
                    if i == best_soln_index:
                        best_soln_string = ",".join(soln_result_list)
                
                
                logger.info("Generation %s Best Individual:", self.generation.current)
                logger.info(', '.join(archive_header)+'\n'+best_soln_string+'\n')
                
        ## Create restart file
                logger.debug("Writing restart file %s...",self.input.job_name+".rst")
                with open(self.input.job_name+".rst", "wb") as f: # Open in binary write mode
                    pickle.dump(self, f)

        ## Clear solution files to save disk space
                if self.input.clear_results == "all":
                    logger.info("Clearing solution files for Generation 0...")
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_0_Indv_*')
                    logger.info("Done!\n")
                elif self.input.clear_results == "all_but_best":
                    logger.info("Clearing all but best solution files for Generation 0...")
                    os.system(f'mv ./{self.input.results_dir_name}/Gen_0_Indv_{best_soln_index} ./{self.input.results_dir_name}/safeGen_0_Indv_{best_soln_index}')
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_0_Indv_*')
                    os.system(f'mv ./{self.input.results_dir_name}/safeGen_0_Indv_{best_soln_index} ./{self.input.results_dir_name}/Gen_0_Indv_{best_soln_index}')
                    logger.info("Done!\n")
                elif self.input.clear_results == "output_files":
                    logger.info(f"Clearing all output files for Generation {self.generation.current}...")
                    for i in range(len(self.population.current)):
                        try:
                            os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.inp ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/safe.inp ')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.*')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/boc_exp.dep')
                            os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/safe.inp  ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.inp')
                        except: 
                            continue
                    logger.info("Done!\n")
        
    ## restart previous optimization routine ##
            else:
            ## Initialize optimization variables
                pool = Pool(processes=self.input.num_procs) #initialize parallel execution
                archive_header = ["Generation","Individual","Fitness Value"]
                for param in self.input.objectives.keys():
                    archive_header.append(str(param))
                archive_header.append("Chromosome")
            ## Check that results files exist
                cwd = Path(os.getcwd())
                results_dir = cwd.joinpath(self.input.results_dir_name)
                if not os.path.exists(results_dir):
                    logger.debug("Results directory is missing. Creating results directory...")
                    os.mkdir(results_dir)
                if not os.path.exists(cwd.joinpath("optimizer_results.csv")):
                    logger.debug("'optimizer_results.csv' results file is missing and will be recreated.")
                    ## write output file
                    with open("optimizer_results.csv", 'w') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
                        csvwriter.writerow(archive_header)
            
    ## Iterate Over Generations  ##

            for iter_gens in range(1,self.generation.total):
                self.generation.current += 1
            ## Create new generation
                logger.info("Creating population of %s individuals for generation %s...", self.input.population_size, self.generation.current)
                new_chromosome_list = self.algorithm.reproduction(self.population.current, self.generation.current)
                self.population.current = []
                for i in range(len(new_chromosome_list)):
                    self.population.current.append(self.generate_solution(f'Gen_{self.generation.current}_Indv_{i}', new_chromosome_list[i]))

            ## Evaluate fitness
                ## If chromosome exists in previous generations, skip call to external model.
                inactive_solutions = []
                for soln in self.population.current:
                    try:
                        soln_index = self.population.archive['solutions'].index(soln.chromosome)
                        soln.fitness_value = self.population.archive['fitnesses'][soln_index]
                        soln.parameters = self.population.archive['parameters'][soln_index]
                        inactive_solutions.append(soln)
                        self.population.current.remove(soln)
                        logger.debug(f"Fitness value for solution '{soln.name}' will be taken from archive entry: {soln_index}.")
                    except ValueError:
                        continue #chromosome is unique, do nothing.
                logger.info("Calculating fitness for generation %s...", self.generation.current)
                ## Execute and parse objective/constraint values
                self.population.current = pool.starmap(self.eval_func, zip(self.population.current, repeat(self.input)))
                if 'cost_fuelcycle' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, self.input)
                if 'av_fuelenrichment' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, self.input)
                
                ## Calculate fitness from objective/constriant values
                for soln in self.population.current:
                    soln.fitness_value = self.fitness.calculate(soln.parameters)
                logger.info("Done!")
                
                ## Recombine active and inactive solutions.
                for soln in inactive_solutions:
                    self.population.current.append(soln)

            ## Archive results
                for soln in self.population.current:
                    self.population.archive['solutions'].append(soln.chromosome)
                    self.population.archive['fitnesses'].append(soln.fitness_value)
                    self.population.archive['parameters'].append(soln.parameters)
                
                best_soln_index = [s.fitness_value for s in self.population.current].index(max([s.fitness_value for s in self.population.current]))
                for i in range(len(self.population.current)):
                    soln = self.population.current[i]
                    soln_result_list = [str(self.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                    for param in soln.parameters.keys():
                        if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                            soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                        else:
                            try:
                                soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                            except:
                                soln_result_list.append(str(soln.parameters[param]['value']))
                    for gene in soln.chromosome:
                        soln_result_list.append(str(gene))
                    ## write to output file
                    with open("optimizer_results.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_list)
                    if i == best_soln_index:
                        best_soln_string = ",".join(soln_result_list)
                
                logger.info("Generation %s Best Individual:", self.generation.current)
                logger.info(', '.join(archive_header)+'\n'+best_soln_string+'\n')
                
            ## Create restart file
                logger.debug("Rewriting restart file %s...",self.input.job_name+".rst")
                with open(self.input.job_name+".rst", "wb") as f: # Open in binary write mode
                    pickle.dump(self, f)
                
            ## Clear solution files to save disk space
                if self.input.clear_results == "all":
                    logger.info(f"Clearing solution files for Generation {self.generation.current}...")
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_*')
                    logger.info("Done!\n")
                elif self.input.clear_results == "all_but_best":
                    logger.info(f"Clearing all but best solution files for Generation {self.generation.current}...")
                    os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{best_soln_index} ./{self.input.results_dir_name}/safeGen_{self.generation.current}_Indv_{best_soln_index}')
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_*')
                    os.system(f'mv ./{self.input.results_dir_name}/safeGen_{self.generation.current}_Indv_{best_soln_index} ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{best_soln_index}')
                    logger.info("Done!\n")
                elif self.input.clear_results == "output_files":
                    logger.info(f"Clearing all output files for Generation {self.generation.current}...")
                    for i in range(len(self.population.current)):
                        try:
                            os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.inp ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/safe.inp ')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.*')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/boc_exp.dep')
                            os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/safe.inp  ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{i}/Gen_{self.generation.current}_Indv_{i}.inp')
                        except: 
                            continue
                    logger.info("Done!\n")
            
                terminate = self.termination_criteria.TC_methods(self.population.current, self.input.termination_criteria)
                if terminate == True: 
                    logger.info("--Run terminated due to termination criterion being met--\n")
                    break

            ## Optimization concluded
            # Report best solution information in output file
            optimization_information = optools.Solution_Reporting()
            best_solution_info = optimization_information.best_solution_information("optimizer_results.csv")
            logger.info("Best solution found in optimization: \n")
            for key in best_solution_info:
                logger.info(f'{key}: {best_solution_info[key]}')
            
            statistics = optimization_information.optimization_statistics()
            last_gen_key = max(statistics.keys(), key=lambda x: int(x.split('_')[1]))  # Extracts last generation key
            last_gen_data = statistics[last_gen_key]  # Retrieves the statistics for that generation
            #Print statistics of last generation to output file
            logger.info("\nStatistics for last generation: \n")
            for key in last_gen_data:
                logger.info(f'{key}: {last_gen_data[key]}')
            
            #Create output statistics file
            with open('optimization_statistics.csv','w') as file:
                file.write('Generation, Average Fitness, Maximum Fitness, Standard Deviation of Fitness\n')
                for gen in statistics:
                    avg_fit = statistics[gen]['Average_Fitness']
                    max_fit = statistics[gen]['Average_Fitness']
                    std_fit = statistics[gen]['Std_Fitness']
                    file.write(f'{int(gen.split('_')[1])},{avg_fit},{max_fit},{std_fit}\n')
            #Plot statistics info if user turned on plot
            if self.input.statistics_plots:
                optimization_information.plot_optimization_statistics()
            #Plot convergence if user turned on convergence
            if self.input.convergence_plot:
                optimization_information.plot_optimization_convergence()
        
            return
        

        elif self.input.methodology == 'simulated_annealing' and self.input.num_procs > 1:
            """
            This functions as the main execution logic for PSA. PSA requires a unique function because the algorithm structure
            is very different from the other evoloutionary algorithms. Creating its own execution logic was the best way to keep everything clean.
            Future developers should do everything they can to integrate their algorithms into the original execution logic.

            TODO add restart ability with PSA 

            Written by Jake Mikouchi. 04/25/2025
            """
        ## Initialize logging for the present file
            logger = logging.getLogger("MIDAS_logger")
            
            if not restart:
        ## Create results directory and/or clear old
                cwd = Path(os.getcwd())
                results_dir = cwd.joinpath(self.input.results_dir_name)
                if os.path.exists(results_dir):
                    logger.debug("Overwriting existing results directory...")
                    rmtree(results_dir, ignore_errors=True)
                    os.mkdir(results_dir)
                else:
                    os.mkdir(results_dir)
            
        ## Delete old results file, if necessary, to avoid confusion
                # this step is useful since the new results file won't be written until after Generation 0 is completed.
                if os.path.exists("optimizer_results.csv"):
                    logger.debug("Removing existing 'optimizer_results.csv' results file...")
                    os.remove("optimizer_results.csv")
                
        ## Initialize beginning population
                logger.info("Generating Buffer of %s individuals...", self.input.buffer_size)
                self.population.current = []
                for i in range(self.input.buffer_size):
                    self.population.current.append(self.generate_solution(f'Buffer_Indv_{i}'))
                
                pool = Pool(processes=self.input.num_procs) #initialize parallel execution
                
        ## Evaluate fitness
                logger.info("Calculating fitness of buffer solutions")
                ## Execute and parse objective/constraint values
                self.population.current = pool.starmap(self.eval_func, zip(self.population.current, repeat(self.input)))
                if 'cost_fuelcycle' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, self.input)
                if 'av_fuelenrichment' in self.input.objectives.keys():
                    for soln in self.population.current:
                        soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, self.input)
                ## Calculate fitness from objective/constriant values
                for soln in self.population.current:
                    soln.fitness_value = self.fitness.calculate(soln.parameters)
                logger.info("Done!")
                    
        ## Archive initial results
                for soln in self.population.current:
                    self.population.archive['solutions'].append(soln.chromosome)
                    self.population.archive['fitnesses'].append(soln.fitness_value)
                    self.population.archive['parameters'].append(soln.parameters)
                
                ## Only initialize the results file the first time.
                archive_header = ["Generation","processor","Fitness Value"]
                for param in self.input.objectives.keys():
                    archive_header.append(str(param))
                archive_header.append("Chromosome")
                ## write output file
                with open("optimizer_results.csv", 'w') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
                    csvwriter.writerow(archive_header)
                
                best_soln_index = [s.fitness_value for s in self.population.current].index(max([s.fitness_value for s in self.population.current]))
                for i in range(len(self.population.current)):
                    soln = self.population.current[i]
                    soln_result_list = [str(self.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                    for param in soln.parameters.keys():
                        if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                            soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                        else:
                            try:
                                soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                            except:
                                soln_result_list.append(str(soln.parameters[param]['value']))
                    for gene in soln.chromosome:
                        soln_result_list.append(str(gene))
                    ## write to output file
                    with open("optimizer_results.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_list)
                    if i == best_soln_index:
                        best_soln_string = ",".join(soln_result_list)
                
        ## Create restart file
                logger.debug("Writing restart file %s...",self.input.job_name+".rst")
                with open(self.input.job_name+".rst", "wb") as f: # Open in binary write mode
                    pickle.dump(self, f)

        ## Clear solution files to save disk space
                if self.input.clear_results == "all":
                    logger.info("Clearing solution files for Buffer...")
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Buffer_Indv_*')
                    logger.info("Done!\n")
                elif self.input.clear_results == "all_but_best":
                    logger.info("Clearing all but best solution files for Buffer...")
                    os.system(f'mv ./{self.input.results_dir_name}/Buffer_Indv_{best_soln_index} ./{self.input.results_dir_name}/safeGen_0_Indv_{best_soln_index}')
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Buffer_Indv_*')
                    os.system(f'mv ./{self.input.results_dir_name}/safeGen_0_Indv_{best_soln_index} ./{self.input.results_dir_name}/Buffer_Indv_{best_soln_index}')
                    logger.info("Done!\n")
                elif self.input.clear_results == "output_files":
                    logger.info(f"Clearing all output files for Buffer...")
                    for i in range(self.input.buffer_size):
                        try:
                            os.system(f'mv ./{self.input.results_dir_name}/Buffer_Indv_{i}/Buffer_Indv_{i}.inp ./{self.input.results_dir_name}/Buffer_Indv_{i}/safe.inp ')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Buffer_Indv_{i}/Buffer_Indv_{i}.*')
                            os.system(f'rm -rf ./{self.input.results_dir_name}/Buffer_Indv_{i}/boc_exp.dep')
                            os.system(f'mv ./{self.input.results_dir_name}/Buffer_Indv_{i}/safe.inp  ./{self.input.results_dir_name}/Buffer_Indv_{i}/Buffer_Indv_{i}.inp')
                        except: 
                            continue
                    logger.info("Done!\n")

    ## Iterate Over Generations  ##
            for iter_gens in range(self.generation.total):
                ## create buffer / update buffer
                self.algorithm.update_buffer(self.population.current, self.generation.current)
                self.generation.current += 1
                logger.info("Running generation %s with %s parallel SA algorithms for %s iterations....", self.generation.current, self.input.num_procs, self.input.population_size)
                for population_step in range(self.input.population_size):
                    new_chromosome_list = []
                    for i in range(self.input.num_procs):
                        new_chromosome_list.append(self.algorithm.reproduction(population_step , i))

                    self.population.current = []
                    for i in range(self.input.num_procs):
                        proc_data = []
                        for j in range(len(new_chromosome_list[i])):
                            proc_data.append(self.generate_solution(f'Gen_{self.generation.current}_Indv_{population_step}_proc_{i}', new_chromosome_list[i][j]))
                        self.population.current.append(proc_data)

                ## Evaluate fitness
                    ## If chromosome exists in previous generations, skip call to external model.
                    inactive_solutions = []
                    for i in range(self.input.num_procs):
                        for soln in self.population.current[i]:
                            try:
                                soln_index = self.population.archive['solutions'].index(soln.chromosome)
                                soln.fitness_value = self.population.archive['fitnesses'][soln_index]
                                soln.parameters = self.population.archive['parameters'][soln_index]
                                inactive_solutions.append(soln)
                                self.population.current[i].remove(soln)
                                logger.debug(f"Fitness value for solution '{soln.name}' will be taken from archive entry: {soln_index}.")
                            except ValueError:
                                continue #chromosome is unique, do nothing.

                    ## Execute and parse objective/constraint values
                    soln_index = []
                    calculations = []
                    for i in range(self.input.num_procs):
                        soln_index.append(len(self.population.current[i]))
                        for j in self.population.current[i]:
                            calculations.append(j)

                    calculations = pool.starmap(self.eval_func, zip(calculations, repeat(self.input)))

                    x = 0
                    for i in range(self.input.num_procs):
                        for j in range(soln_index[i]):
                            self.population.current[i][j] = calculations[x]
                            x += 1

                    for i in range(self.input.num_procs):
                        if 'cost_fuelcycle' in self.input.objectives.keys():
                            for soln in self.population.current[i]:
                                soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, self.input)
                        if 'av_fuelenrichment' in self.input.objectives.keys():
                            for soln in self.population.current[i]:
                                soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, self.input)

                    ## Calculate fitness from objective/constriant values
                    for i in range(self.input.num_procs):
                        for soln in self.population.current[i]:
                            soln.fitness_value = self.fitness.calculate(soln.parameters)
                    
                    ## Recombine active and inactive solutions.
                    for soln in inactive_solutions:
                        self.population.current[inactive_solutions.index(soln)].append(soln)

                    self.algorithm.update_active(self.population.current)

                                    ## Archive results
                    for i in range(self.input.num_procs):
                        for soln in self.population.current[i]:
                            self.population.archive['solutions'].append(soln.chromosome)
                            self.population.archive['fitnesses'].append(soln.fitness_value)
                            self.population.archive['parameters'].append(soln.parameters)
                logger.info("Done!")
                
                best_indv_fitness = -10000000
                for i in range(self.input.num_procs):  
                    if max([s.fitness_value for s in self.population.current[i]]) > best_indv_fitness:
                        best_indv_fitness = max([s.fitness_value for s in self.population.current[i]])
                        best_soln_index = (i, [s.fitness_value for s in self.population.current[i]].index(max([s.fitness_value for s in self.population.current[i]])))

                for i in range(self.input.num_procs):  
                    for j in range(len(self.population.current[i])):
                        soln = self.population.current[i][j]
                        soln_result_list = [str(self.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                        for param in soln.parameters.keys():
                            if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                                soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                            else:
                                try:
                                    soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                                except:
                                    soln_result_list.append(str(soln.parameters[param]['value']))
                        for gene in soln.chromosome:
                            soln_result_list.append(str(gene))
                        ## write to output file
                        with open("optimizer_results.csv", 'a') as csvfile:
                            csvwriter = csv.writer(csvfile, delimiter=',')
                            csvwriter.writerow(soln_result_list)
                        if i == best_soln_index[0] and j == best_soln_index[1]:
                            best_soln_string = ",".join(soln_result_list)
                    
                logger.info("Generation %s Best Individual:", self.generation.current)
                logger.info(', '.join(archive_header)+'\n'+best_soln_string+'\n')

            ## Create restart file
                logger.debug("Rewriting restart file %s...",self.input.job_name+".rst")
                with open(self.input.job_name+".rst", "wb") as f: # Open in binary write mode
                    pickle.dump(self, f)
                
            ## Clear solution files to save disk space
                if self.input.clear_results == "all":
                    logger.info(f"Clearing solution files for Generation {self.generation.current}...")
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}*')
                    logger.info("Done!\n")
                elif self.input.clear_results == "all_but_best":
                    logger.info(f"Clearing all but best solution files for Generation {self.generation.current}...")
                    os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{best_soln_index[1]}_proc{best_soln_index[0]} ./{self.input.results_dir_name}/safeGen_{self.generation.current}_Indv_{best_soln_index[1]}_proc{best_soln_index[0]}')
                    os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_*')
                    os.system(f'mv ./{self.input.results_dir_name}/safeGen_{self.generation.current}_Indv_{best_soln_index[1]}_proc{best_soln_index[0]} ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{best_soln_index[1]}_proc{best_soln_index[0]}')
                    logger.info("Done!\n")
                elif self.input.clear_results == "output_files":
                    logger.info(f"Clearing all output files for Generation {self.generation.current}...")
                    for i in range(self.input.num_procs):
                        for population_step in range(self.input.population_size):
                            try:
                                os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}.inp ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/safe.inp ')
                                os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}.*')
                                os.system(f'rm -rf ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/boc_exp.dep')
                                os.system(f'mv ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/safe.inp  ./{self.input.results_dir_name}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}/Gen_{self.generation.current}_Indv_{population_step}_proc_{i}.inp')
                            except: 
                                continue
                    logger.info("Done!\n")
            
                terminate = self.termination_criteria.TC_methods(self.population.current, self.input.termination_criteria)
                if terminate == True: 
                    logger.info("--Run terminated due to termination criterion being met--\n")
                    break

            ## Optimization concluded
            # Report best solution information in output file
            optimization_information = optools.Solution_Reporting()
            best_solution_info = optimization_information.best_solution_information("optimizer_results.csv")
            logger.info("Best solution found in optimization: \n")
            for key in best_solution_info:
                logger.info(f'{key}: {best_solution_info[key]}')
            
            statistics = optimization_information.optimization_statistics()
            last_gen_key = max(statistics.keys(), key=lambda x: int(x.split('_')[1]))  # Extracts last generation key
            last_gen_data = statistics[last_gen_key]  # Retrieves the statistics for that generation
            #Print statistics of last generation to output file
            logger.info("\nStatistics for last generation: \n")
            for key in last_gen_data:
                logger.info(f'{key}: {last_gen_data[key]}')
            
            #Create output statistics file
            with open('optimization_statistics.csv','w') as file:
                file.write('Generation, Average Fitness, Maximum Fitness, Standard Deviation of Fitness\n')
                for gen in statistics:
                    avg_fit = statistics[gen]['Average_Fitness']
                    max_fit = statistics[gen]['Average_Fitness']
                    std_fit = statistics[gen]['Std_Fitness']
                    file.write(f'{int(gen.split('_')[1])},{avg_fit},{max_fit},{std_fit}\n')
            #Plot statistics info if user turned on plot
            if self.input.statistics_plots:
                optimization_information.plot_optimization_statistics()
            #Plot convergence if user turned on convergence
            if self.input.convergence_plot:
                optimization_information.plot_optimization_convergence()
        
            return

        