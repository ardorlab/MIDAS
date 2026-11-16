from copy import deepcopy
from zipfile import Path

from midas.codes import surrogatemodel_BK
from midas.utils import optimizer_tools as optools
from midas.utils import problem_preparation
import midas_data
import csv
import pickle
from shutil import rmtree
from midas.codes import parcs343
from midas.utils import LWR_fuelcyclecost
from midas.utils import LWR_averageenrichment
from surmodel.trainmodelRPF import trainmodelrpf
from surmodel.trainmodelcore import traincoremodel_update
import joblib
from itertools import repeat

import os
import logging
import time as TT

## Initialize logging for the present file
logger = logging.getLogger("MIDAS_logger")

def evaluate(solution, mc_input, xsdict):
    """
    interface for evaluating a multicycle solution by running each cycle through the optimizer and collecting the results.

    C.King 2026-11
    """

    # initialize parameter result collector
    cycle_param_results = {}
    for obj_key in mc_input.objectives.keys():
        cycle_param_results[obj_key] = {}

    print("MULTICYCLE INPUT:")
    print("Multicycle Chromosome:",solution.chromosome)
    print("Multicycle Output directory:",mc_input.results_dir_name)
    # print("Num cycles:", mc_input.num_cycles)
    # print("Optimizer:",mc_input.methodology)
    # print("Pop Size:",mc_input.population_size)
    # print("Num Gens:",mc_input.num_generations)
    # print("Objectives:",mc_input.objectives)
    # print("Perturbs:",mc_input.perturbation_type)   
    print("Multicycle Assembly Options:",mc_input.fa_options)
    print("Multicycle Assembly Params/Genome:",mc_input.genome)

    for cycle_num in range(mc_input.num_cycles):
        # include a more interpretable printout of each 
        logger.info(f"Evaluating CYCLE {cycle_num+1} of {mc_input.num_cycles+1}")
        cycle_input = deepcopy(mc_input) # copy all settings from the input object to avoid modifying the original multicycle problem
        
        # Change settings to be for single-cycle. Not all settings need to be changed.
        # the current multicycle problem is not affected as this a copy

        # cycle_input.calculation_type = "single_cycle" # switch to single-cycle for evaluating the single-cycle problem
        # we don't have to switch to single_cycle since 

        # Some of the inputs for multicycle can be copied over to single-cycle. Others need to be adjusted
        # Since the inputs have been copied, we do not have to worry about breaking the original multicycle problem

        # GENERAL BLOCK
        cycle_input.calculation_type = "single_cycle" # set to single cycle
        cycle_input.methodology = "simulated_annealing"
        # set initial_temperature and other PSA defaults
        # code_type: surrogate (maybe?)
        # results_dreictory_name: smoething not output_files
        cycle_input.results_dir_name = "output_files" + "_cycle_" + str(cycle_num + 1) # might need to adjust this for each cycle that is made
        # BIG ISSUE: might blow up file size
        
        # OPTIMIZER BLOCK
        cycle_input.population_size = 1 # pop size is 1 for simulated annealing
        cycle_input.num_generations = 100 # set a reasonable number of generations for simulated annealing
        cycle_input.objectives = {'max_boron': {'goal': 'less_than_target', 'target': 1300, 'weight': 1.0},
                                  'PinPowerPeaking': {'goal': 'less_than_target', 'target': 5000, 'weight': 2.1},
                                  'FDeltaH': {'goal': 'less_than_target', 'target': 1.48, 'weight': 6000},
                                  'cycle_length': {'goal': 'maximize', 'weight': 1.0}                                  
        }

        # ALGORITHM BLOCK
        cycle_input.perturbation_type = {"method": "perturb_by_swap", "num_permutations": 4}

        cycle_input.initial_temperature = 100
        if cycle_input.num_procs > 1:
            cycle_input.cooling_schedule = 'lam' 
        if cycle_input.num_procs <= 1:
            cycle_input.cooling_schedule = 'exponential_decrease'
            if cycle_input.cooling_schedule == "lam":
                raise ValueError(f"Cooling schedule '{cycle_input.cooling_schedule}' only available for parallel simulated annealing.")
        cycle_input.secondary_cooling_schedule = 'exponential_decrease'

        print(f"CYCLE {cycle_num+1} INPUT:")
        cycle_input.genome = deepcopy(solution.single_cycle_from_multicycle(cycle_num,cycle_input.fa_options)) # single_cycle_from_mutlicycle modifies asssembly_optionsdirectly
        # add burn assemblies will have to be in the process above!
        # parcs_data can remain the same
        
        print("Output directory:",cycle_input.results_dir_name)
        # print("Num cycles:", cycle_input.num_cycles) # this may be an issue for anything that uses the number of cycles as a check
        # print("Optimizer:",cycle_input.methodology)
        # print("Pop Size:",cycle_input.population_size)
        # print("Num Gens:",cycle_input.num_generations)
        # print("Objectives:",cycle_input.objectives)
        # print("Perturbs:",cycle_input.perturbation_type)
        print("Assembly Options:",cycle_input.fa_options)


        print("Running LP map loader")
        optools.LP_map_loader(cycle_input)

        print("Assembly Params/Genome:",cycle_input.genome)

        print("Calc Type:",cycle_input.calculation_type)


        cycle_input = problem_preparation.Prepare_Problem_Values.prepare_calc(cycle_input) # prepare calculation values
        run_single_cycle(cycle_input) # run the single-cycle optimization for this cycle
        # extract info later

    print("MULTICYCLE INPUT MODIFICATION CHECK")
    print("Chromosome:",solution.chromosome)
    print("Output directory:",mc_input.results_dir_name)
    # print("Num cycles:",mc_input.num_cycles)
    # print("Optimizer:",mc_input.methodology)
    # print("Pop Size:",mc_input.population_size)
    # print("Num Gens:",mc_input.num_generations)
    # print("Objectives:",mc_input.objectives)
    # print("Perturbs:",mc_input.perturbation_type)
    print("Assembly Options:",mc_input.fa_options)
    print("Assembly Params/Genome:",mc_input.genome)

    # Sum up results of each cycle
    for item in cycle_param_results:
        cycle_param_results[item] = 5
        solution.parameters[item]['value'] = cycle_param_results[item]

    # process cycle_param_results to get final numbers for solution fitness function params

    return solution # return solution object and return to standard MIDAS execution

def run_single_cycle(cycle_input):
    """
    Runs a single-cycle optimization using PSA with the MIDAS optimization technique.
    Primarily copied from Optimizer.py() - parallelization cannot be nested so the pool technique leads to an error

    """ 
    from midas.optimizer import Optimizer # import in function to prevent circular input

    print("Calc Type:",cycle_input.calculation_type)
    print(f"Optimizer is activated! Beware!")
    optimizer2 = Optimizer(cycle_input) # Run the optimizer with the cycle_input parameters
    
    print("Calc type after optimizer initialized:",cycle_input.calculation_type)
    optimizer2.build_optimizer()

    logger.info(f"Completed secondary Optimizer assembly. Running secondary optimizer...\n")
    ts = TT.time()

    # equivalent of optimizer2.main()

    restart = False
    ## Initialize logging for the present file
    logger = logging.getLogger("MIDAS_logger")
    if not restart:
## Create results directory and/or clear old
        cwd = Path(os.getcwd())
        results_dir = cwd.joinpath(optimizer2.input.results_dir_name)
        if os.path.exists(results_dir):
            logger.debug("Overwriting existing results directory...")
            rmtree(results_dir)
            os.makedirs(results_dir)
        else:
            os.makedirs(results_dir)
        ## cleaning retrain data if any
        if os.path.exists(midas_data.__path_to_store_retrain_data__):
            rmtree(midas_data.__path_to_store_retrain_data__)
    
## Delete old results file, if necessary, to avoid confusion
        # this step is useful since the new results file won't be written until after Generation 0 is completed.
        if os.path.exists(midas_data.__odir__ + "/optimizer_results.csv"):
            logger.debug("Removing existing 'optimizer_results.csv' results file...")
            os.remove(midas_data.__odir__ + "/optimizer_results.csv")
        
## Initialize beginning population
        optimizer2.population.current = []
        namedir = []
        if optimizer2.input.methodology == 'simulated_annealing' and optimizer2.input.num_procs > 1:
            # for parallel simulated annealing the initial population is the size of the buffer
            logger.info("Generating initial population of %s individuals...", optimizer2.input.buffer_size)
            for i in range(optimizer2.input.buffer_size):
                chromosome = optimizer2.get_initial_population(i) 
                optimizer2.population.current.append(optimizer2.generate_solution(f'Gen_0_Indv_{i}', chromosome))
        else:
            logger.info("Generating initial population of %s individuals...", optimizer2.input.population_size)
            for i in range(optimizer2.population.size):
                chromosome = optimizer2.get_initial_population(i) 
                optimizer2.population.current.append(optimizer2.generate_solution(f'Gen_0_Indv_{i}', chromosome))
                namedir.append(f'Gen_0_Indv_{i}')

        
        #if optimizer2.input.code_interface == "surrogate":
        #    import multiprocessing as mp
        #    mp.set_start_method('spawn', force=True) ### NOTE: do not use this for debug mode!!!!!
        #    from surmodel.paralel_MLmodel import init_worker, init_worker_updated
        #    chunksize = max(1, optimizer2.population.size // (optimizer2.input.num_procs * 4))
        xsdict = joblib.load(midas_data.__path_xs_pickle__)#,mmap_mode='r')
        #    ### set pool once 
        #    pool = Pool(processes=optimizer2.input.num_procs, initializer=init_worker)
        #else:
        #    pool = Pool(processes=optimizer2.input.num_procs) #initialize parallel execution

        if optimizer2.input.methodology == 'simulated_annealing'and optimizer2.input.num_procs == 1:
            logger.info(f"Initial Temperature: {optimizer2.input.initial_temperature}")
        
## Evaluate fitness
        logger.info("Calculating fitness for generation %s...", optimizer2.generation.current)
        ## Execute and parse objective/constraint values
        ### initial with PARCRs running first 
        if optimizer2.generation.current in optimizer2.input.parcs_steps:
            optimizer2.eval_func = parcs343.evaluate
            print('Run parcs initially ', optimizer2.generation.current )
            # if 'pool' in locals():
            #        pool.close()
            #        pool.join()
            #        del pool
            #        gc.collect()
            #        pool = Pool(processes=optimizer2.input.num_procs)
            optimizer2.population.current= optimizer2.eval_func(optimizer2.population.current, optimizer2.input)
            temppop = optimizer2.population.current
            surrogatemodel_BK.extract_data(namedir,optimizer2.input, xsdict) ## extract data after PARCs done
            # # temparcs = optimizer2.population.current
            if optimizer2.generation.current in optimizer2.input.retrain_steps:
                trainmodelrpf(optimizer2.input.if_calibrate)
                traincoremodel_update(optimizer2.input.if_calibrate)
                if optimizer2.input.keep_retraindata==False: ## delete retrain data after retrainning
                    rmtree(midas_data.__path_to_store_retrain_data__)
                ## update pool 
                # if 'pool' in locals():
                #    pool.close()
                #    pool.join()
                #    del pool
                #    gc.collect()
                #    pool = Pool(processes=optimizer2.input.num_procs, initializer=init_worker_updated)
                # ## for printing out fitness value at this step 
                optimizer2.eval_func = surrogatemodel_BK.evaluate
                surrogate_sol = optimizer2.eval_func(temppop,optimizer2.input,xsdict)

                for soln in surrogate_sol:
                    soln.fitness_value = optimizer2.fitness.calculate(soln.parameters)
                for i in range(len(temppop)):
                    soln = surrogate_sol[i]
                    soln_result_list = [str(optimizer2.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                    for param in soln.parameters.keys():
                        if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                            soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                        else:
                            soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                    for gene in soln.chromosome:
                        soln_result_list.append(str(gene))
                    ## write to output file
                    with open("optimizer_results_at_retrain.csv", 'a') as csvfile:
                        csvwriter = csv.writer(csvfile, delimiter=',')
                        csvwriter.writerow(soln_result_list)
            optimizer2.eval_func = surrogatemodel_BK.evaluate ## reset
       
            # for i in range(optimizer2.population.size):
            #     optimizer2.population.current[i] = optimizer2.eval_func(optimizer2.population.current[i], optimizer2.input) ## serial test
    
        else:
            optimizer2.population.current = optimizer2.eval_func(optimizer2.population.current, optimizer2.input, xsdict)
        
        # if 'cost_fuelcycle' in optimizer2.input.objectives.keys():
        #     for soln in optimizer2.population.current:
        #         soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, optimizer2.input)
        # if 'av_fuelenrichment' in optimizer2.input.objectives.keys():
        #     for soln in optimizer2.population.current:
        #         soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, optimizer2.input)
        ## Calculate fitness from objective/constriant values
        for soln in optimizer2.population.current:
            soln.fitness_value = optimizer2.fitness.calculate(soln.parameters)
        logger.info("Done!")
    
## Archive initial results
        for soln in optimizer2.population.current:
            optimizer2.population.archive['solutions'].append(soln.chromosome)
            optimizer2.population.archive['fitnesses'].append(soln.fitness_value)
            optimizer2.population.archive['parameters'].append(soln.parameters)
        
        ## Only initialize the results file the first time.
        archive_header = ["Generation","Individual","Fitness Value"]
        for param in optimizer2.input.objectives.keys():
            archive_header.append(str(param))
        archive_header.append("Chromosome")
        ## write output file
        with open(midas_data.__odir__ + "/optimizer_results.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
            csvwriter.writerow(archive_header)
        
        best_soln_index = [s.fitness_value for s in optimizer2.population.current].index(max([s.fitness_value for s in optimizer2.population.current]))
        for i in range(len(optimizer2.population.current)):
            soln = optimizer2.population.current[i]
            soln_result_list = [str(optimizer2.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
            for param in soln.parameters.keys():
                if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                    soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                else:
                    soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
            for gene in soln.chromosome:
                soln_result_list.append(str(gene))
            ## write to output file
            with open(midas_data.__odir__ + "/optimizer_results.csv", 'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',')
                csvwriter.writerow(soln_result_list)
            if i == best_soln_index:
                best_soln_string = ",".join(soln_result_list)
        
        logger.info("Generation %s Best Individual:", optimizer2.generation.current)
        logger.info(', '.join(archive_header)+'\n'+best_soln_string+'\n')
        
## Create restart file
        logger.debug("Writing restart file %s...",optimizer2.input.job_name+".rst")
        with open(optimizer2.input.job_name+".rst", "wb") as f: # Open in binary write mode
            pickle.dump(optimizer2, f)
    
## Clear solution files to save disk space
        if optimizer2.input.clear_results == "all":
            logger.info("Clearing solution files for Generation 0...")
            os.system(f'rm -rf ./{optimizer2.input.results_dir_name}/Gen_0_Indv_*')
            logger.info("Done!\n")
        elif optimizer2.input.clear_results == "all_but_best":
            logger.info("Clearing all but best solution files for Generation 0...")
            os.system(f'mv ./{optimizer2.input.results_dir_name}/Gen_0_Indv_{best_soln_index} ./{optimizer2.input.results_dir_name}/safeGen_0_Indv_{best_soln_index}')
            os.system(f'rm -rf ./{optimizer2.input.results_dir_name}/Gen_0_Indv_*')
            os.system(f'mv ./{optimizer2.input.results_dir_name}/safeGen_0_Indv_{best_soln_index} ./{optimizer2.input.results_dir_name}/Gen_0_Indv_{best_soln_index}')
            logger.info("Done!\n")
        
## Iterate Over Generations  ##

    for iter_gens in range(1,optimizer2.generation.total):
        optimizer2.generation.current += 1
    ## Create new generation
        logger.info("Creating population of %s individuals for generation %s...", optimizer2.input.population_size, optimizer2.generation.current)
        new_chromosome_list = optimizer2.algorithm.reproduction(optimizer2.population.current, optimizer2.generation)

        if optimizer2.input.methodology == 'simulated_annealing' and optimizer2.input.num_procs > 1: # if parallel simulated annealing skip evaluation
            optimizer2.population.current = []
            for i in range(len(new_chromosome_list)):
                if os.path.exists(f"{optimizer2.input.results_dir_name}/{new_chromosome_list[i].name}"): 
                    os.rename(f"{optimizer2.input.results_dir_name}/{new_chromosome_list[i].name}", f'{optimizer2.input.results_dir_name}/Gen_{optimizer2.generation.current}_Indv_{i}')
                new_chromosome_list[i].name = f'Gen_{optimizer2.generation.current}_Indv_{i}'
                optimizer2.population.current.append(new_chromosome_list[i])
            logger.info("Calculating fitness for generation %s...", optimizer2.generation.current)
            logger.info("Done!")

        else: #every other algorithm
            optimizer2.population.current = []
            namedir = []
            for i in range(len(new_chromosome_list)):
                optimizer2.population.current.append(optimizer2.generate_solution(f'Gen_{optimizer2.generation.current}_Indv_{i}', new_chromosome_list[i]))
                namedir.append(f'Gen_{optimizer2.generation.current}_Indv_{i}')
        
        ## Evaluate fitness
            ## If chromosome exists in previous generations, skip call to external code.
            ### testtttttttttt
            # inactive_solutions = []
            # for soln in optimizer2.population.current:
            #     try:
            #         soln_index = optimizer2.population.archive['solutions'].index(soln.chromosome)
            #         soln.fitness_value = optimizer2.population.archive['fitnesses'][soln_index]
            #         soln.parameters = optimizer2.population.archive['parameters'][soln_index]
            #         inactive_solutions.append(soln)
            #         ## skip the removal for test??
            #         optimizer2.population.current.remove(soln)
            #         ## take out some in active solution 
            #         namedir.remove(soln.name)
            #         logger.debug(f"Fitness value for solution '{soln.name}' will be taken from archive entry: {soln_index}.")
            #     except ValueError:
            #         continue #chromosome is unique, do nothing.
            
            logger.info("Calculating fitness for generation %s...", optimizer2.generation.current)
            ## Execute and parse objective/constraint values
            if optimizer2.generation.current in optimizer2.input.parcs_steps:
                optimizer2.eval_func = parcs343.evaluate
                print('Run parcs initially ', optimizer2.generation.current )
                # if 'pool' in locals():
                #    pool.close()
                #    pool.join()
                #    del pool
                #    gc.collect()
                #    pool = Pool(processes=optimizer2.input.num_procs)
                ## new pool 
                optimizer2.population.current= optimizer2.eval_func(optimizer2.population.current, optimizer2.input)
                temppop = optimizer2.population.current 
                # ## extract data after PARCS done 
                surrogatemodel_BK.extract_data(namedir, optimizer2.input, xsdict)
                if optimizer2.generation.current in optimizer2.input.retrain_steps:
                    trainmodelrpf(optimizer2.input.if_calibrate)
                    traincoremodel_update(optimizer2.input.if_calibrate)
                    if optimizer2.input.keep_retraindata==False: ## delete retrain data after retrainning
                        rmtree(midas_data.__path_to_store_retrain_data__)
                    # update pool
                    #if 'pool' in locals():
                    #    pool.close()
                    #    pool.join()
                    #    del pool
                    #    gc.collect()
                    #    pool = Pool(processes=optimizer2.input.num_procs, initializer=init_worker_updated)
                    ## for printing out fitness value at this step 
                    optimizer2.eval_func = surrogatemodel_BK.evaluate
                    surrogate_sol = optimizer2.eval_func(temppop, optimizer2.input, xsdict)
                    
                    for soln in surrogate_sol:
                        soln.fitness_value = optimizer2.fitness.calculate(soln.parameters)
                    for i in range(len(temppop)):
                        soln = surrogate_sol[i]
                        soln_result_list = [str(optimizer2.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
                        for param in soln.parameters.keys():
                            if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                                soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                            else:
                                soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
                        for gene in soln.chromosome:
                            soln_result_list.append(str(gene))
                        ## write to output file
                        with open("optimizer_results_at_retrain.csv", 'a') as csvfile:
                            csvwriter = csv.writer(csvfile, delimiter=',')
                            csvwriter.writerow(soln_result_list)

                optimizer2.eval_func = surrogatemodel_BK.evaluate ## reset
            else:
                optimizer2.population.current = optimizer2.eval_func(optimizer2.population.current, optimizer2.input, xsdict)

            if 'cost_fuelcycle' in optimizer2.input.objectives.keys():
                for soln in optimizer2.population.current:
                    soln.parameters = LWR_fuelcyclecost.get_fuelcycle_cost(soln, optimizer2.input)
            if 'av_fuelenrichment' in optimizer2.input.objectives.keys():
                for soln in optimizer2.population.current:
                    soln.parameters = LWR_averageenrichment.get_avfuelenrichment(soln, optimizer2.input)
            
            ## Calculate fitness from objective/constriant values
            for soln in optimizer2.population.current:
                soln.fitness_value = optimizer2.fitness.calculate(soln.parameters)
            logger.info("Done!")
            
            ## Recombine active and inactive solutions.
            # for soln in inactive_solutions:
            #     optimizer2.population.current.append(soln)
    
    ## Archive results
        for soln in optimizer2.population.current:
            optimizer2.population.archive['solutions'].append(soln.chromosome)
            optimizer2.population.archive['fitnesses'].append(soln.fitness_value)
            optimizer2.population.archive['parameters'].append(soln.parameters)
        
        best_soln_index = [s.fitness_value for s in optimizer2.population.current].index(max([s.fitness_value for s in optimizer2.population.current]))
        for i in range(len(optimizer2.population.current)):
            soln = optimizer2.population.current[i]
            soln_result_list = [str(optimizer2.generation.current),str(i),'{0:.3f}'.format(soln.fitness_value)]
            for param in soln.parameters.keys():
                if param == 'av_fuelenrichment': #reformat this parameter prior to printing
                    soln_result_list.append('{0:.3f}'.format(100*soln.parameters[param]['value'])) #convert w.t. to wo%
                else:
                    soln_result_list.append('{0:.3f}'.format(soln.parameters[param]['value']))
            for gene in soln.chromosome:
                soln_result_list.append(str(gene))
            ## write to output file
            with open(midas_data.__odir__ + "/optimizer_results.csv", 'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',')
                csvwriter.writerow(soln_result_list)
            if i == best_soln_index:
                best_soln_string = ",".join(soln_result_list)
        
        logger.info("Generation %s Best Individual:", optimizer2.generation.current)
        logger.info(', '.join(archive_header)+'\n'+best_soln_string+'\n')
        
    ## Create restart file
        logger.debug("Rewriting restart file %s...",optimizer2.input.job_name+".rst")
        with open(optimizer2.input.job_name+".rst", "wb") as f: # Open in binary write mode
            pickle.dump(optimizer2, f)
        
    ## Clear solution files to save disk space
        if optimizer2.input.clear_results == "all":
            logger.info(f"Clearing solution files for Generation {optimizer2.generation.current}...")
            os.system(f'rm -rf ./{optimizer2.input.results_dir_name}/Gen_{optimizer2.generation.current}_Indv_*')
            logger.info("Done!\n")
        elif optimizer2.input.clear_results == "all_but_best":
            logger.info(f"Clearing all but best solution files for Generation {optimizer2.generation.current}...")
            os.system(f'mv ./{optimizer2.input.results_dir_name}/Gen_{optimizer2.generation.current}_Indv_{best_soln_index} ./{optimizer2.input.results_dir_name}/safeGen_{optimizer2.generation.current}_Indv_{best_soln_index}')
            os.system(f'rm -rf ./{optimizer2.input.results_dir_name}/Gen_{optimizer2.generation.current}_Indv_*')
            os.system(f'mv ./{optimizer2.input.results_dir_name}/safeGen_{optimizer2.generation.current}_Indv_{best_soln_index} ./{optimizer2.input.results_dir_name}/Gen_{optimizer2.generation.current}_Indv_{best_soln_index}')
            logger.info("Done!\n")
        terminate = optimizer2.termination_criteria.TC_methods(optimizer2.population.current, optimizer2.input.termination_criteria)
        if terminate == True: 
            logger.info("--Run terminated due to termination criterion being met--\n")
            break

    ## Optimization concluded
    ## cleanning
    #if 'pool' in locals():
    #    pool.close()
    #    pool.join()
    #    gc.collect()
    # Report best solution information in output file
    optimization_information = optools.Solution_Reporting()
    best_solution_info = optimization_information.best_solution_information(midas_data.__odir__ + "/optimizer_results.csv")
    logger.info("Best solution found in optimization: \n")
    for key in best_solution_info:
        logger.info(f'{key}: {best_solution_info[key]}')
        
    with open(midas_data.__odir__ + "/best_solution_chromosome.csv", 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
        csvwriter.writerow(best_solution_info['Chromosome'])
    
    statistics = optimization_information.optimization_statistics()
    last_gen_key = max(statistics.keys(), key=lambda x: int(x.split('_')[1]))  # Extracts last generation key
    last_gen_data = statistics[last_gen_key]  # Retrieves the statistics for that generation
    #Print statistics of last generation to output file
    logger.info("\nStatistics for last generation: \n")
    for key in last_gen_data:
        logger.info(f'{key}: {last_gen_data[key]}')
    ## cleaning and close pool 
    #gc.collect() 
    #Create output statistics file
    with open(midas_data.__odir__ + "/optimization_statistics.csv", 'w') as file:
        file.write('Generation, Average Fitness, Maximum Fitness, Standard Deviation of Fitness\n')
        for gen in statistics:
            avg_fit = statistics[gen]['Average_Fitness']
            max_fit = statistics[gen]['Average_Fitness']
            std_fit = statistics[gen]['Std_Fitness']
            file.write(f"{int(gen.split('_')[1])},{avg_fit},{max_fit},{std_fit}\n")
    #Plot statistics info if user turned on plot
    if optimizer2.input.statistics_plots:
        optimization_information.plot_optimization_statistics()
    #Plot convergence if user turned on convergence
    if optimizer2.input.convergence_plots:
        optimization_information.plot_optimization_convergence()



    # nevertheless, this can be resolved later
    logger.info(f'Execution completed after {TT.time()-ts:.5f} seconds')

    
    return