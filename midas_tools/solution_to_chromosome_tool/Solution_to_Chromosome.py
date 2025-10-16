"""
Created by Jake Mikouchi 10/2/2025

This is a tool to convert user created solutions into a chromosome so that they can be used as an initial solution in an optimization.
Specific optimization types like LP or control blade sequence optimization will require special handeling to convert but generic 
optimizations like discrete and continuous can hopefully be done generally.

Conversions that are currently available:
Single Cycle:
    BWR:
        printing symmetry: full     optimization symmetry: full
        printing symmetry: quarter  optimization symmetry: quarter
        printing symmetry: quarter  optimization symmetry: octant
    PWR:
        printing symmetry: full     optimization symmetry: full
        printing symmetry: quarter  optimization symmetry: quarter
        printing symmetry: quarter  optimization symmetry: octant
Equilibrium Cycle:
    BWR:
        printing symmetry: quarter  optimization symmetry: quarter
        printing symmetry: quarter  optimization symmetry: octant
    PWR:
        printing symmetry: quarter  optimization symmetry: quarter
Lattice Physics:
    BWR:
        lattice symmetry: full
        lattice symmetry: diagonal
    PWR: 
        lattice symmetry: full
        lattice symmetry: se


HOW TO USE:
    set up txt file:   
    Simply copy the solutions as they would appear in the input file for the respective problem types and paste into a txt file. 
    For lattices copy the pinmap from the polaris input file.
    For single cycles copy the rad_conf from the parcs input file.
    In both lattice and single cycles just the map is sufficient. 
    For equilibrium cycles copy the location and shuf_maps from the parcs input file.
    In the case of equilibrium cycles, the LOCATION and SHUF_MAP flags need to be present in the txt file. 
    It is easiest to copy them both at the same time and paste without any modifications. 
    The location map must appear first in the txt file.

    set up problem info:
    Fill in the problem info dictionary with all of the necessary information. 
    All of the inputs have the same options as midas. 

    set up parameter types: 
    For single cycle and equilibrium cycle optimizations, fill in "assemblies" element with information about chromosome. 
        The first element in the tuple is the gene name and the second is how PARCS represents the assembly (type).
        This must be defined in the exact same way as in the midas .yaml file under 'assembly_options' or the optimization will not work.
    For lattice optimizations, fill in "pins" element with information about the chromosome.
        The format is the same as the "assemblies" element and the information must be the same as "rod_options" in the midas .yaml file.
    For equilibrium optimizations, the "batches" element must be filled in along with "assemblies".
        To do this the name of each batch needs to be given in order from first batch to last. 
        These must be defined in the same way as the batches in the "assembly_parameters" option in the midas .yaml file.

"""
import csv
import copy


problem_info = {
    "calc_type": "single_cycle", # same options as midas
    "file_location": "./PWR_single_cycle.txt",

    # LP optimization requirements, same options as midas
    "LWR_System": "PWR", # needed for lattice too 
    "LP_printing_symmetry": "quarter", 
    "LP_optimization_symmetry": "quarter",

    # LP lattice requirements
    "lattice_symmetry": "full" # same options as midas
}


parameter_types = {
                #[(MIDAS gene, PARCS assembly), ...]
    "assemblies": [("2ENR", 20), ("3ENR", 30), ("4ENR", 40), ("5ENR", 50)],
                #[(MIDAS gene, polaris pin), ...]
    "pins": [("FUEL1", 1), ("UGAD4", 4), ("UGAD8", 8), ("INSTRTUBE", 2)],
                # write batches in order from fresh to last batch
    "batches": ["batch 0", "batch 1", "batch 2"],

    "continuous": [],
    "discrete": []

}

class Solution_to_Chromosome():

    def distribute_conversion(problem_info, parameter_types):
        pt = problem_info["calc_type"]
        if pt == "single_cycle":
            chromosome = Solution_to_Chromosome.convert_single_cycle_LWR(problem_info, parameter_types)
        elif pt == "eq_cycle":
            chromosome = Solution_to_Chromosome.convert_eq_cycle_LWR(problem_info, parameter_types)
        elif pt == "lattice_physics":
            chromosome = Solution_to_Chromosome.convert_lattice(problem_info, parameter_types)
        else:
            raise ValueError(f"Conversion {pt} not currently supported")


        Solution_to_Chromosome.store_chromosome(chromosome, pt)

    def store_chromosome(chromosome, pt):
        if pt != "eq_cycle":
            with open("solution_to_chromosme.csv", 'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONE)
                csvwriter.writerow(chromosome)
        else: 
            with open("solution_to_chromosme.csv", 'a') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter=',')
                csvwriter.writerow(chromosome)

    def convert_single_cycle_LWR(problem_info, parameter_types):
        core_map = []
        chromosome = []
        # retirve chromosome
        with open(problem_info["file_location"], "r") as f:
            file = f.readlines()
            for line in file: 
                core_map.append([ int(assemb) for assemb in line.split()])

        # remove empty assemblies and reflectors
        for i in range(len(core_map)):
            core_map[i] = list(filter(LP_conversion_tools.is_assembly, core_map[i]))

        # remove whitespaces
        core_map= list(filter(None, core_map))

        # create preprocessed chromosome
        if problem_info["LWR_System"].lower() == "bwr":           
            if problem_info["LP_printing_symmetry"].lower() == "full":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    for i in range(len(core_map)):
                        chromosome.extend(core_map[i])
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to quarter core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to octant core not yet available")

            elif problem_info["LP_printing_symmetry"].lower() == "quarter":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to expand out to full core
                    raise ValueError("Expanding from quarter core to full core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    for i in range(len(core_map)):
                        chromosome.extend(core_map[i])
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    for i in range(len(core_map)):
                        if (i+1) < len(core_map[i]):
                            chromosome.extend(core_map[i][:(i+1)])
                        else: 
                            chromosome.extend(core_map[i])
        
        if problem_info["LWR_System"].lower() == "pwr":           
            if problem_info["LP_printing_symmetry"].lower() == "full":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    for i in range(len(core_map)):
                        chromosome.extend(core_map[i])
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to quarter core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to octant core not yet available")

            elif problem_info["LP_printing_symmetry"].lower() == "quarter":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to expand out to full core
                    raise ValueError("Expanding from quarter core to full core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    for i in range(len(core_map)):
                        if i == 0:
                            chromosome.extend([core_map[i][i]])
                        else:
                            chromosome.extend(core_map[i])
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    for i in range(len(core_map)):
                        if (i+1) < len(core_map[i]):
                            chromosome.extend(core_map[i][:(i+1)])
                        else: 
                            chromosome.extend(core_map[i])

        chromosome = LP_conversion_tools.process_single_cycle_chromosome(parameter_types, chromosome)

        return chromosome

    def convert_eq_cycle_LWR(problem_info, parameter_types):
        location_map = []
        shuf_map = []
        chromosome = []
        # retirve chromosome
        location_map_flag = False 
        shuf_map_flag = False
        with open(problem_info["file_location"], "r") as f:
            file = f.readlines()
            for line in file: 
                if "location" in line.lower():
                    location_map_flag = True 
                    shuf_map_flag = False
                if "shuf_map" in line.lower():
                    location_map_flag = False 
                    shuf_map_flag = True        
                if location_map_flag and ("location" not in line.lower() and len(line.strip()) != 0):              
                    location_map.append([ str(assemb) for assemb in line.split()])
                if shuf_map_flag and ("shuf_map" not in line.lower() and len(line.strip()) != 0):              
                    shuf_map.append([ str(assemb) for assemb in line.split()])

        for i in range(len(shuf_map)):
            for j in range(len(shuf_map[i])):
                if shuf_map[i][j] == '0' or shuf_map[i][j] == '00': 
                    shuf_map[i][j] = 0 
                if location_map[i][j] == '0' or location_map[i][j] == '00':
                    location_map[i][j] = 0

        # remove empty assemblies and reflectors
        for i in range(len(shuf_map)):
            shuf_map[i] = list(filter(LP_conversion_tools.is_assembly, shuf_map[i]))
            location_map[i] = list(filter(LP_conversion_tools.is_assembly, location_map[i]))

        # remove whitespaces
        shuf_map= list(filter(None, shuf_map))
        location_map = list(filter(None, location_map))

        solution_holder = copy.deepcopy(shuf_map)
        for i in range(len(shuf_map)):
            for j in range(len(shuf_map[i])):
                solution_holder[i][j] = 0                

        # create preprocessed chromosome
        if problem_info["LWR_System"].lower() == "bwr":           
            if problem_info["LP_printing_symmetry"].lower() == "full":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to read full core
                    raise ValueError("Reading from full core values not currently available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to quarter core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to octant core not yet available")
            
            elif problem_info["LP_printing_symmetry"].lower() == "quarter":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to expand out to full core
                    raise ValueError("Expanding from quarter core to full core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    # create chromosome
                    for i in range(len(shuf_map)):
                        for j in range(len(shuf_map[i])):
                            chromosome.append([0,0])
                        
                    # fill in batch 0 by looking for fresh fuel locations
                    index = 0
                    for i in range(len(shuf_map)):
                        for j in range(len(shuf_map[i])):
                            try: 
                                if isinstance(int(shuf_map[i][j]), int):
                                    chromosome[index][0] = parameter_types["batches"][0]
                                    fresh_fuel_type_index = [assembly[1] for assembly in parameter_types["assemblies"]].index(abs(int(shuf_map[i][j])))
                                    chromosome[index][1] = parameter_types["assemblies"][fresh_fuel_type_index][0]
                                    solution_holder[i][j] = 10
                            except: 
                                pass
                            index += 1

                    batch_holder = 10
                    for batch_index in range(len(parameter_types["batches"][1:])):
                        index = 0
                        for i in range(len(shuf_map)):
                            for j in range(len(shuf_map[i])):
                                # try block to determine if location is not fresh fuel
                                try:
                                    isinstance(int(shuf_map[i][j]), int)
                                except:
                                    # for loop to find row and column index of assembly in location_map
                                    for loc_row_index, row in enumerate(location_map):
                                        try:
                                            loc_col_index = row.index(shuf_map[i][j])
                                            break
                                        except:
                                            pass
                                    # contruct gene which contains batch number and index from which the assembly came from
                                    if solution_holder[loc_row_index][loc_col_index] == batch_holder:
                                        chromosome[index][0] = parameter_types["batches"][batch_index+1]
                                        chromosome[index][1] =  sum(len(row) for row in solution_holder[:loc_row_index]) + loc_col_index 
                                        solution_holder[i][j] = batch_holder + 10 
                                index += 1
                        batch_holder += 10 
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    # create chromosome
                    for i in range(len(shuf_map)):
                        if (i+1) < len(shuf_map[i]):
                            chromosome.extend([[0,0] for j in range(len(shuf_map[i][:(i+1)]))])
                        else: 
                            chromosome.extend([[0,0] for j in range(len(shuf_map[i]))])
                        
                    # fill in batch 0 by looking for fresh fuel locations
                    index = 0
                    for i in range(len(shuf_map)):
                        for j in range(len(shuf_map[i])):
                            if j <= i:
                                try: 
                                    if isinstance(int(shuf_map[i][j]), int):
                                        chromosome[index][0] = parameter_types["batches"][0]
                                        fresh_fuel_type_index = [assembly[1] for assembly in parameter_types["assemblies"]].index(abs(int(shuf_map[i][j])))
                                        chromosome[index][1] = parameter_types["assemblies"][fresh_fuel_type_index][0]
                                        solution_holder[i][j] = 10
                                except: 
                                    pass
                                index += 1

                    batch_holder = 10
                    for batch_index in range(len(parameter_types["batches"][1:])):
                        index = 0
                        for i in range(len(shuf_map)):
                            for j in range(len(shuf_map[i])):
                                if j <= i:
                                    # try block to determine if location is not fresh fuel
                                    try:
                                        isinstance(int(shuf_map[i][j]), int)
                                    except:
                                        # for loop to find row and column index of assembly in location_map
                                        for loc_row_index, row in enumerate(location_map):
                                            try:
                                                loc_col_index = row.index(shuf_map[i][j])
                                                break
                                            except:
                                                pass
                                        # contruct gene which contains batch number and index from which the assembly came from
                                        if solution_holder[loc_row_index][loc_col_index] == batch_holder:
                                            chromosome[index][0] = parameter_types["batches"][batch_index+1]
                                            solution_index_sum = 0 
                                            stop = False
                                            for k in range(len(solution_holder)):
                                                for l in range(len(solution_holder[k])):   
                                                    if l <= k: 
                                                        solution_index_sum += 1
                                                    if k == loc_row_index and l == loc_col_index:
                                                        stop = True
                                                        break
                                                if stop:
                                                    break
                                            chromosome[index][1] =  solution_index_sum -1 #+ loc_col_index 

                                            solution_holder[i][j] = batch_holder + 10 
                                    index += 1
                        batch_holder += 10 

        if problem_info["LWR_System"].lower() == "pwr":           
            if problem_info["LP_printing_symmetry"].lower() == "full":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to read full core
                    raise ValueError("Reading from full core values not currently available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to quarter core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from full core to octant core not yet available")

            elif problem_info["LP_printing_symmetry"].lower() == "quarter":
                if problem_info["LP_optimization_symmetry"].lower() == "full":
                    #TODO add ability to expand out to full core
                    raise ValueError("Expanding from quarter core to full core not yet available")
                if problem_info["LP_optimization_symmetry"].lower() == "quarter":
                    # create chromosome
                    for i in range(len(shuf_map)):
                        for j in range(len(shuf_map[i])):
                            if (i == 0 and j == 0) or i > 0:
                                chromosome.append([0,0])

                    # fill in batch 0 by looking for fresh fuel locations
                    index = 0
                    for i in range(len(shuf_map)):
                        for j in range(len(shuf_map[i])):
                            if (i == 0 and j == 0) or i > 0:
                                try: 
                                    if isinstance(int(shuf_map[i][j]), int):
                                        chromosome[index][0] = parameter_types["batches"][0]
                                        fresh_fuel_type_index = [assembly[1] for assembly in parameter_types["assemblies"]].index(abs(int(shuf_map[i][j])))
                                        chromosome[index][1] = parameter_types["assemblies"][fresh_fuel_type_index][0]
                                        solution_holder[i][j] = 10
                                except: 
                                    pass
                                index += 1

                    batch_holder = 10
                    for batch_index in range(len(parameter_types["batches"][1:])):
                        index = 0
                        for i in range(len(shuf_map)):
                            for j in range(len(shuf_map[i])):
                                if (i == 0 and j == 0) or i > 0:
                                    # try block to determine if location is not fresh fuel
                                    try:
                                        isinstance(int(shuf_map[i][j]), int)
                                    except:
                                        # for loop to find row and column index of assembly in location_map
                                        for loc_row_index, row in enumerate(location_map):
                                            try:
                                                loc_col_index = row.index(shuf_map[i][j])
                                                break
                                            except:
                                                pass
                                        # contruct gene which contains batch number and index from which the assembly came from
                                        if solution_holder[loc_row_index][loc_col_index] == batch_holder:
                                            chromosome[index][0] = parameter_types["batches"][batch_index+1]
                                            chromosome[index][1] =  sum(len(row) for row in solution_holder[:loc_row_index]) + loc_col_index 
                                            solution_holder[i][j] = batch_holder + 10 
                                    index += 1
                        batch_holder += 10 
                if problem_info["LP_optimization_symmetry"].lower() == "octant":
                    #TODO add ability to expand out to full core
                    raise ValueError("Compressing from quarter core to octant core not yet available")

        chromosome = LP_conversion_tools.process_eq_cycle_chromosome(parameter_types, chromosome)
        
        return chromosome
    
    def convert_lattice(prblem_info, parameter_types):
        lattice_map = []
        chromosome = []
        # retirve chromosome
        with open(problem_info["file_location"], "r") as f:
            file = f.readlines()
            for line in file: 
                lattice_map.append([ int(pin) for pin in line.split()])

        # remove whitespaces
        lattice_map = list(filter(None, lattice_map))
        if problem_info["LWR_System"].lower() == "pwr": 
            if problem_info["lattice_symmetry"].lower() == "full":
                for i in range(len(lattice_map)):
                    chromosome.extend(lattice_map[i])
            if problem_info["lattice_symmetry"].lower() == "se":
                for i in range(len(lattice_map)):
                    chromosome.extend(lattice_map[i][:(i+1)])
        elif problem_info["LWR_System"].lower() == "bwr": 
            if problem_info["lattice_symmetry"].lower() == "full":
                for i in range(len(lattice_map)):
                    chromosome.extend(lattice_map[i])
            if problem_info["lattice_symmetry"].lower() == "diagonal":
                for i in range(len(lattice_map)):
                    chromosome.extend(lattice_map[i][:(i+1)])
                

        chromosome = Lat_conversion_tools.process_lattice_chromosome(parameter_types, chromosome)

        return chromosome

class LP_conversion_tools():
    def is_assembly(assembly):
        # checks if an element is not refl or empty in the core map
        checker = True
        if assembly == 10  or assembly == 00:
            checker = False
        return checker 
    
    def process_single_cycle_chromosome(parameter_types, chromosome):
        genes, assemblies = zip(*parameter_types["assemblies"])
        for i in range(len(chromosome)):
            if chromosome[i] in assemblies: 
                chromosome[i] = genes[assemblies.index(chromosome[i])]
            else: 
                raise ValueError(f"'{chromosome[i]}' not found in {parameter_types["assemblies"]}")
        
        return chromosome

    def process_eq_cycle_chromosome(parameter_types, chromosome):
        chromosome_holder = []
        for gene in chromosome:
            new_gene = (gene[0], gene[1])
            chromosome_holder.append(str(new_gene))
        chromosome = chromosome_holder
        
        return chromosome

class Lat_conversion_tools():  
    def process_lattice_chromosome(parameter_types, chromosome):
        genes, assemblies = zip(*parameter_types["pins"])
        for i in range(len(chromosome)):
            if chromosome[i] in assemblies: 
                chromosome[i] = genes[assemblies.index(chromosome[i])]
            else: 
                raise ValueError(f"'{chromosome[i]}' not found in {parameter_types["pins"]}")
        
        return chromosome

Solution_to_Chromosome.distribute_conversion(problem_info, parameter_types)
print("Solution succesfully converted! ")