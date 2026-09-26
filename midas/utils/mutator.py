## Import Block ##
import logging
from copy import deepcopy
import random
from midas.utils import optimizer_tools as optools

class Mutator():
    """
    This class contains all mutation methodologies available in MIDAS.
    Mutators are methods for perturbing solutions within optimization processes. 
    Many Mutators are common between optimizers (even if the name is not the same)
    this class contains all mutators which can be shared between algorithms. 

    Mutators which are algorithm specific should be coded into the specific algorithm class. 

    Written by Jake Mikouchi 09/24/26
    """

    def mutator_methods(input_obj, chromosome):
        """
        Method for distributing to the requested mutation method.
        
        Written by Jake Mikouchi. 09/24/26
        """
        if input_obj.mutation_type['method'] == 'mutate_by_gene':
            child_chromosome = Mutator.mutate_by_gene(input_obj, chromosome)
        elif input_obj.mutation_type['method'] == 'mutate_by_swap':
            child_chromosome = Mutator.mutate_by_swap(input_obj, chromosome)
        else:
            raise ValueError("Requested mutation type not recognized.")

        return child_chromosome
    
    def mutate_by_gene(input_obj, chromosome):
        """
        Generates a new solution by randomly assigning a new value to a selected gene.
        Users can specifiy the number of genes to undergo assignment.
        
        Updated by Nicholas Rollins. 09/27/2024
        Updated by Jake Mikouchi. 09/24/2026
        """
       ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")
        
        core_parameters = [input_obj.nrow, input_obj.ncol, input_obj.num_assemblies,
                            input_obj.symmetry, input_obj.calculation_type]
        
        if input_obj.calculation_type in ["eq_cycle"]:
            zone_chromosome = [loc[0] for loc in chromosome]
            child_zone_chromosome = deepcopy(zone_chromosome)
            old_soln = zone_chromosome
            new_soln = child_zone_chromosome
            all_gene_options = input_obj.batches
            all_genes_list = list(input_obj.batches.keys())
        else:
            child_chromosome = deepcopy(chromosome)
            old_soln = chromosome
            new_soln = child_chromosome
            all_gene_options = input_obj.genome
            all_genes_list = list(input_obj.genome.keys())

        num_mutations = input_obj.mutation_type['num_mutations'] 
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            new_soln = deepcopy(old_soln) #in the case of abortion, start from scratch.
            while new_soln == old_soln:
                for i in range(num_mutations):
                    loc_to_mutate = random.randint(0, len(new_soln)-1) #choose a random gene
                    old_gene = new_soln[loc_to_mutate]
                    gene_options = optools.Gene_Validity_check.contraceptive_check(input_obj, all_genes_list, all_gene_options,
                                                                                    core_parameters, old_soln, [], loc_to_mutate)
                    if gene_options == [0,1]:
                        new_gene = random.uniform(0,1)
                    else:
                        new_gene = random.choice(gene_options)
                    if new_gene != old_gene:
                        if input_obj.calculation_type in ["single_cycle","eq_cycle", "lattice_physics"] and all_gene_options[new_gene]['map'][loc_to_mutate] == 1:
                            new_soln[loc_to_mutate] = new_gene
                        else:
                            new_soln[loc_to_mutate] = new_gene
            chromosome_is_valid = optools.Gene_Validity_check.abortive_check(input_obj, all_genes_list,all_gene_options,\
                                                                            core_parameters,new_soln)
            if not chromosome_is_valid:
                attempts += 1
                if attempts > 100000:
                    logger.error("Mutate-by-Gene has failed after 100,000 attempts; the Individual will be restored. Consider relaxing the constraints on the input space.")
                    return chromosome

        if input_obj.calculation_type in ["eq_cycle"]:
            #recreate child_chromosome
            child_chromosome = []
            for i in range(len(new_soln)):
                if new_soln[i] == chromosome[i][0]:
                    child_chromosome.append(chromosome[i])
                else:
                    child_chromosome.append((new_soln[i],None))
            child_chromosome = optools.Solution.EQ_reload_fuel(input_obj.genome,core_parameters,child_chromosome)

        else:
            child_chromosome = new_soln
            
        return child_chromosome

    def mutate_by_swap(input_obj, chromosome):
        """
        Generates a new solution by randomly selecting two genes and swapping their values.
        Users can specifiy the number of swaps are performed.
        
        Written by Jake Mikouchi. 09/24/2026
        """
       ## Initialize logging for the present file
        logger = logging.getLogger("MIDAS_logger")
        
        core_parameters = [input_obj.nrow, input_obj.ncol, input_obj.num_assemblies,
                            input_obj.symmetry, input_obj.calculation_type]
        
        if input_obj.calculation_type in ["eq_cycle"]:
            zone_chromosome = [loc[0] for loc in chromosome]
            child_zone_chromosome = deepcopy(zone_chromosome)
            old_soln = zone_chromosome
            new_soln = child_zone_chromosome
            all_gene_options = input_obj.batches
            all_genes_list = list(input_obj.batches.keys())
        else:
            child_chromosome = deepcopy(chromosome)
            old_soln = chromosome
            new_soln = child_chromosome
            all_gene_options = input_obj.genome
            all_genes_list = list(input_obj.genome.keys())

        num_mutations = input_obj.mutation_type['num_mutations'] 
        chromosome_is_valid = False
        attempts = 0
        while not chromosome_is_valid:
            new_soln = deepcopy(old_soln) #in the case of abortion, start from scratch.
            for i in range(num_mutations):
                loc_1 = random.randint(0, len(new_soln)-1) #choose a random gene
                loc_2 = random.randint(0, len(new_soln)-1)
                while loc_1 == loc_2:
                    loc_2 = random.randint(0, len(new_soln)-1)

                gene_1 = new_soln[loc_1]
                gene_2 = new_soln[loc_2]
                new_soln[loc_1] = gene_2
                new_soln[loc_2] = gene_1

            chromosome_is_valid = optools.Gene_Validity_check.abortive_check(input_obj, all_genes_list,all_gene_options,\
                                                                            core_parameters,new_soln)
            if not chromosome_is_valid:
                attempts += 1
                if attempts > 100000:
                    logger.error("Mutate-by-Gene has failed after 100,000 attempts; the Individual will be restored. Consider relaxing the constraints on the input space.")
                    return chromosome

        if input_obj.calculation_type in ["eq_cycle"]:
            #recreate child_chromosome
            child_chromosome = []
            for i in range(len(new_soln)):
                if new_soln[i] == chromosome[i][0]:
                    child_chromosome.append(chromosome[i])
                else:
                    child_chromosome.append((new_soln[i],None))
            child_chromosome = optools.Solution.EQ_reload_fuel(input_obj.genome,core_parameters,child_chromosome)

        else:
            child_chromosome = new_soln
        
        # import pdb; pdb.set_trace()
            
        return child_chromosome