This is the Modular Optimization Framework repository. It utilizes inheritance, object-oriented, and functional programming to
create a simple, robust tool for solving optimization problems. It has been applied primarily to nuclear engineering design problems. 


# MOF: Modular Optimization Framework

MOF is designed to provide users with a variety of optimization methodologies to solve opimization problems with a focus on nuclear engineering design problems. Containing multiple optimization methodologies in a single package allows for the reuse of code in multiple ways leading to a shorter, simpler, and more versatile optimization package.

Current optimization methodologies supported by the MOF are:

* Genetic Algorithm
* Simulated Annealing
  

# Code Installation

It is highly advised to install Miniconda or Anaconda. This will allow you to create a controlled Python environment where you can
install the required packages, especially if you want to use it in a cluster with limited permissions. Go to the 
site: https://docs.conda.io/en/latest/miniconda.html and download the latest Python 3 installer. The installer is a bash file with an example name "miniconda_install.sh". Now install conda and the 
required dependencies entering the following commands:

    bash miniconda_install.sh

    pip install pyyaml 

    conda install numpy

    conda install matplotlib
	
	conda install h5py

    git clone https://github.com/gkdelipei/MOF.git

Congratulations. The code is now installed in your local machine.
  
# Running the Code 

Enter your local MOF directory and navigate to one of the sample problems. Type:

    python mofMain.py --input sample_problem_input.yaml --cpus 4 

Running the code is as sample as that. mofMain is the main function of the python code. There are two 
commands here. The first is __--input__. This command is used to designate the input yaml file that you want to run.
In the example it is just the input yaml file corresponding to the selected sample case. The second command
__--cpus__ is used to designate how many processors you want to apply to the optimization problem. In the example four are specified.

# YAML Input File Format

So now that you have run the code, let's look at what settings we used. The settings for the optimization are input 
through a yaml file. YAML is a high level markup file. YAML stands for Yeah, Aint Markup Language. It's a very simple,
easy to read and write file that can be directly loaded into python as a dictionary data tool, making it extremely easy
to use and perfect for our purposes here.

Open up the example.yaml file and take a look at it. Obviously use the file editor of your preference. If you don't 
know any file editors or are looking for one, we highly recommed VSCode. 

There is generally a wide range of flexibility in the YAML files from input to input, which can make it complicated to
use the optimization program. However, there are several things that will always be consistent. This section will 
detail the general outline of the yaml files to give you a brief exposure.

Yaml files may be written in any order. Like Python, indentation is used to seperate different layers in the file. 
Here, the markers in the example input file are explained:

    optimization: 
    "The top marker in the input file. Designates that the cards underneath are related to the optimization settings."

        methodology: genetic_algorithm
        "The methodology marker designates which type of optimization is being performed, in this case genetic algorithm."

        population_size: 8
        "Indicates population size used in the genetic algorithm."

        number_of_generations: 5
        "Indicates number of generations over which optimization will run." 

        mutation:
        "Marker indicating sub markers are related to how mutation is performed."
            
            method: [mutate_by_common,mutate_fixed]
            "Specifies how mutation is performed. In this case two mutation methods have been selected."

            initial_rate: 0.25
            "The starting percentage of population that undergo mutation."

            final_rate: 0.75
            "The final percentage of population that undergo mutation."

        fixed_problem: True
        "Indicates that the genome is held fixed in some way."

        number_fixed_groups: 4
        "Indicates that there are four genome groups used to fix the optimization problem."

        fixed_genes_per_group: [17,14,16,17]
        "How many genes must fit into each group in the solution. Note that the order of these groups is the order
        that the gene groups appear when the genes are specified."

        selection:
        "Indicates that the sub markers will be related to how selection is performed using the genetic algorithm." 
		
            fitness: ranked
            "Marker for selecting which fitness/scoring function to use to compare solutions. In this case the fitness is derived from how the solutions are ranked from best to worst."

            method: tournament
            "Indicates solutions are going to be compared using a tournament method." 

        data_type: loading_pattern
        "data_type designates the type of problem that is to be optimized. IN this case the fuel loading pattern of a nuclear reactor."

        objectives:
        "Markers under this marker indicate what objectives will be taken into account in the fitness functions."

            assembly_power:
            "Example of an optimization objective. In this case assembly radial peaking factors."
                goal: minimize
                "Says that radial peaking factors should be minimized in objective function."
            
	genome:
	"Markers under this marker are used to describe the actual problem that is to be solved."

	    chromosomes:
	    "Markers under this marker describe the genes used in the optimization."

	        Assembly_One: 
		    "The first gene in the example problem. Genes are directly under the chromosome marker, and may use any name."

		        gene_group: 2.0
		        "We specified earlier that our optimization problem is fixed. The gene_group is used to identify common genes that
		        are held fixed, i.e., labels that this falls into the first group that is only allowed to have 17 assemblies of 
		        this label. gene_groups can have any name, but the order in which they appear corresponds to the numbering in the 
		        fixed_genes_per_group marker."

		        type: 2
		        "For simulate problems, the type marker corresponds to an assembly designator used in simulate." 

		        name: 2.0_w/o
		        "The Name card is used for solution plotting and for designating casmo lattice types. If you are working on the NE
		        412/512 project, you don't need to worry about these cards and they can have any name you would like." 

		        map: &ID001
				  [1, 1, 1, 1, 1, 1, 1, 1, 0,
				   1, 1, 1, 1, 1, 1, 1, 1, 0,
				   1, 1, 1, 1, 1, 1, 1, 0, 0,     
				   1, 1, 1, 1, 1, 1, 1, 0,      
				   1, 1, 1, 1, 1, 1, 0, 0,            
				   1, 1, 1, 1, 1, 0, 0,                   
				   1, 1, 1, 1, 0, 0,                          
				   1, 1, 0, 0, 0,
				   0, 0, 0]
		       "The map marker specifies where genes may be expressed in the optimization, e.g. locations in the reactor core 
		       where assemblies may be placed. Binary markers are used to indicate whether a gene may or may not be expressed in
		       the location. A 1 indicates the gene may be expressed there, a 0 indicates it may not be expressed there."

		    Assembly_Two:
		        gene_group: 2.5
		        type: 3
		        serial: B300
		        name: 2.5_w/o_no_bp
		        map: *ID001
		    Assembly_Three:
		        gene_group: 3.2
		        type: 5
		        serial: C300
		        name: 3.2_w/o_no_bp
		        map: *ID001
		    Reflector:
		        type: 1
		        gene_group: reflector
		        serial: none
		        name: reflector
		        map: 
			     [0, 0, 0, 0, 0, 0, 0, 0, 1,
			      0, 0, 0, 0, 0, 0, 0, 0, 1,
			      0, 0, 0, 0, 0, 0, 0, 1, 1,     
			      0, 0, 0, 0, 0, 0, 0, 1,      
			      0, 0, 0, 0, 0, 0, 1, 1,            
			      0, 0, 0, 0, 0, 1, 1,                   
			      0, 0, 0, 0, 1, 1,                          
			      0, 0, 1, 1, 1,
			      1, 1, 1]

	    assembly_data:
	    "This marker is used to attach additional information required to run simulate." 

		    type: pwr
		    "The reactor type to be ran in simulate."

		    core_width: 15
		    "The size of the reactor core. If there are 152 assemblies the core width is 15."

		    load_point: 0.000
		    "The point that you want the restart file to load from. "

		    depletion: 20
		    "The max depletion step allowed in the simulate calculations. Note that just because you designate a maximum 
		    depletion doesn't mean you reach that depletion time step. "

		    batch_number: 0
		    "The cycle of the core. For initial loadings of the reactor core you can specify either 0 or 1."

		    pressure: 2250.
		    "The operating pressure of the reactor core."

		    boron: 900.
		    "A guess of the initial critical boron calculation. Need a guess. Doesn't matter what."

		    power: 100.
		    "The percent of rated power that the reactor is operating at."

		    flow: 100.
		    "The percent of rated flow that the reactor is operating at. "

		    inlet_temperature: 550.
		    "The inlet temperature of the coolant."

		    restart_file: s3.pwr.uo2.c02.depl.res
		    "The restart file being used in the simulate analysis."

		    cs_library: pwr.sim_one.lib
		    "The cms_link cross section library being used in the simulate analysis."

# Repository Structure

The repository is structured in the following way:

* samples: Directory including various sample cases for users to get familiarized with MOF framework. Users just need to navigate to the sample folder and run MOF from there with the corresponding input file or use the run.sh bash file on the RDFMG cluster. The samples consist of examples of genetic algorithm and simulated annealing optimization problems for a reduced number of code calculations. All the samples should finish within 10 minutes. Users are encourage to change the number of code evaluations (e.g. generations number) to experiment with optimization algorithms. Parallel execution is implemented only for genetic algorithm and the user could increase the number of allocated processors to evaluate the impact on the execution time.
  * sample_1: Fuel Lattice Optimization with CASMO4 code. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_2: First Cycle Core Loading Pattern Optimization with SIMULATE3. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_3: First Cycle Core Loading Pattern Optimization with SIMULATE3. Simulated annealing is used with 20 iterations. 
  * sample_4: Third Cycle Core Loading Pattern Optimization with SIMULATE3. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_5: Third Cycle Core Loading Pattern Optimization with SIMULATE3. Simulated annealing is used with 20 iterations.
  * dev_samples: Sample cases under development and not fully functional yet.

* rl_dev: Directory containing under development reinforcement learning optimization files. This is not functional yet. This folder will be deleted once reinforcement learning is integrated in the main MOF.

* documentation: Directory including all additional documentation.

* casmo.py: Python file that handles CASMO Lattice calculations evaluation and data extraction.

* crudworks.py: Python file used if CRUD machine learning predictions are required (TensorFlow should be installed).

* fitness.py: Python file for selecting and computing the objective function of the optimization. 

* geneticAlgorithm.py: Python file that stores all classes and functions for performing Genetic Algorithm optimization.

* metrics.py: Python file including tools for tracking solutions and storing the generated optimization data.

* mofMain.py: Python file that is the main body of MOF. In this file the interface between the input file and the optimization is performed by selecting the specified options and initializing all the necessary components.

* simulate.py: Python file that handles SIMULATE calculations evaluation and data extraction.

* simulateAnnealing.py: Python file that stores all classes and functions for performing Simulate Annealing optimization.

* solution_types.py: Python file for storing the solutions of the optimization together with some usefull functions.

* submission_script.sh: Bash file example for running MOF through SLURM on the RDFMG cluster.

