

# MIDAS: Modularly Integrated Design Assistance Suite

<p align="center">
  <img src="https://github.com/ardorlab/MIDAS/assets/61293473/8114a773-988b-493a-9831-8d3f7b80408b" width="300" height="300">
</p>


Welcome to the Modularly Integrated Design Assistance Suite (MIDAS) repository. MIDAS utilizes inheritance, object-oriented, and functional programming to
create a simple, robust tool for solving optimization problems. It has been applied primarily to nuclear engineering design problems, but it has been built to handle all kinds of problems.


MIDAS is designed to provide users with a variety of optimization methodologies to solve opimization problems with a focus on nuclear engineering design. Containing multiple optimization methodologies in a single package allows for the reuse of code in multiple ways leading to a shorter, simpler, and more versatile optimization package.


Current optimization methodologies supported in MIDAS are:

* Genetic Algorithm
* Simulated Annealing
* Parallel Simulated Annealing
* Bayesian Optimization  

  
# Code Installation

It is highly advised to install Miniconda or Anaconda. This will allow you to create a controlled Python environment where you can
install the required packages, especially if you want to use it in a cluster with limited permissions. Go to the 
site: https://docs.conda.io/en/latest/miniconda.html and download the latest Python 3 installer. The installer is a bash file with an example name "miniconda_install.sh". Now install conda and the 
required dependencies entering the following commands:

    bash miniconda_install.sh

	conda install matplotlib  
	conda install pyyaml  
	conda install scipy  
	conda install scikit-learn  
	conda install h5py  
	conda install pyarrow  
	conda install pandas  

    git clone https://github.com/ardorlab/MIDAS.git

<!-- 
If you want to use the newly added reinforcement learning algorithms, the python version in the environment should be 3.9 and some additional dependencies will need to be installed:

    pip3 install torch torchvision torchaudio

    pip install stable-baselines3[extra] 

An alternative way to configure the environment is to use the requirement files provided in the repository for pip and conda tools. This files are the "requirements_pip.txt" and "requirements_conda.txt".
-->

Congratulations. The code is now installed in your local machine.


# Running the Code 

Enter your local MIDAS directory and navigate to one of the sample problems. Type:

    midasmain.py --input midas_input.yaml --cpus 4 

midasmain is the main function of the python code and all secondary functions are disrtibuted from within this script. 
There are two commands here. The first is __--input__. This command is used to designate the input yaml file that you
want to run. In the example it is just the input yaml file corresponding to the selected sample case. The second command 
__--cpus__ is used to designate how many processors you want to apply to the optimization problem. In the example four 
are specified. For clean executions, it is recomended that you execute the code in a designated directory.

Depending on where MIDAS is installed, you will likely need to update the executable paths within MIDAS so that
the modules requiring theses executables can function. This can easily be done by updating the paths at the bottom of 
midas_data.py

# YAML Input File Format

All settings and options needed to execute an optimization can be accessed through the yaml input file. YAML is short for 
'Yeah, Aint Markup Language' and is a high level markup file format. It's a simple, easy to read and write file that can be 
directly loaded into python as a dictionary data tool, making it perfect for interfacing with MIDAS.

A number of sample input files are available in the samples directory to demonstrate how to execute optimizations for 
different types of problems and with different algorithms. Additionally, an extensive user manual is available in the Wiki
of this repository which goes over every input option availble through the yaml file.

There is generally a wide range of flexibility in the YAML files from input to input, which can make it complicated to
use the optimization program. However, there are several things that will always be consistent. This section will 
detail the general outline of the yaml files to give you a brief exposure.

Yaml files may be written in any order. Like Python, indentation is used to seperate different layers in the file. 
Here, the general block structure is explained:

	GENERAL Block
		The general block is used to define the optimization type and the algorithm to be used. Some general details 
		about the execution are also defined here like if statistical plots should be automatically generated when 
		the optimization is complete or the name of the file containing the user defined initial population.

	OPTIMIZATION Block
		The optimization block is used to define optimization parameters which are universal to all algorithms. The
		number of generations and population size are defined in this block as well as termination criteria. The 
		fitness function and it's objectives are also defined in this block.

	ALGORITHM Block
		This block is used to provide MIDAS with the hyperparameter values used with the respective algorithm during 
		optimizations. Each available algorithm has a unique set of hyperparameters that can be defined here. 
		For example, Genetic Algorithm's mutation rate and Simulated Annealing's intial temperature can be defined 
		in this block.

	OPTIONS Block 
		The options block is used to define each of the genes which are avaialble to the optimizers during the 
		optimization. For example, during loading pattern optimizations each unique assembly type is defined in this 
		block. Different problem types have different names given to the options block as its structure may be different
		depending on the needs of the optimization. For Loading Pattern optimizations the options block is called 
		'assembly options', lattice physics problems change the name to 'rod options' and all generic numerical and 
		combinatorial problem types change the name to 'gene options'.

	DECISION VARIABLES Block
		This block is used to further inform MIDAS on design constraints relating to the individual genes. Users use this 
		block to define where genes exist on the chromosome or impose limits on how many times a single gene can appear 
		within the chromosome. This block enables users to improse problem specific requirments to the design such as 
		a minimum number of batch 2 assemblies in an equilibrium core optimization or constraining a certain assembly type 
		by only allowing it to exist in specific locations in the core.

	DATA Block
		The data block is used to provide MIDAS with all problem specific information required to conduct the optimization
		but not necessarily relating to the optimization itself. For optimizations where PARCS is executed to retrieve 
		chromosome performances, this block will contain all the information needed to construct the parcs input files.
		Different problem types have different names given to the data block as the required information is different for 
		each problem type and code interface. For Loading Pattern optimizations the options block is called 'parcs data', 
		lattice physics problems change the name to 'polaris data' and all generic numerical and combinatorial problem types 
		change the name to 'optimization data'.
	

# Repository Structure

The repository is structured in the following way:

* samples: Directory including various sample cases to familiarize users with the MIDAS framework and assist in constructing input files. The available samples are complete and can be immediately run after installation. The samples consist of examples of each of the optimization algorithms for a variety of problem types. The samples contain reduced number of code calculations. Users are encourage to change the number of code evaluations (e.g. generations number or population size) and experiment to with the algorithm hyperparameters. All algorithms except simulated annealing support parallel execution and users are recommended to allocate multiple processors to speed up the optimization especially if the code evaluations are time intensive. For samples using PARCS for code executions, the required crosssections are held in the xslib folder located in the samples directory. 
  
  * BO_ipwr_database.yaml - Bayesian Optimization applied to the ipwr database 
  * GA_single_cycle.yaml - Genetic Algorithm applied to single cycle fresh fuel PWR loading pattern optimization
  * PSA_lattice_physics.yaml - Paralllel Simulated Annealing applied to PWR lattice physics optimization
  * SA_eq_cycle.yaml - Simulated ANnealing applied to equilibrium cycle PWR loading pattern optimization

* MIDAS: Directory containg human facing modules that users interact with while executing the code. This is the root directory of the MIDAS.
  
	* midasmain.py: Python file that is the main body of MIDAS. In this file, administrative tasks are performed such as reading in the input file, writing general information about the optimization to the output file, and initializing the optimization.
  	* midas_data.py: File used mainly for storing file paths to external code executables.

* midas: Directory containing core MIDAS modules which are required to execute an optimization regardless of optimization settings.
  
	* input_parser.py: Module for parsing the yaml input file and returning constructive error messages if the input file is incorrect.
	* optimizer.py: Main execution pathway for MIDAS. Contains the main 'loop' that all algorithms rely on. Distrubutes tasks and required variables to the optimization algorithms and distributes the generated populations to optimizer_tools.py to compute fitness values. All standard output file information is written here.

* midas/utils: Directory containing additional modules and functions to assist MIDAS in interfacing with aptimization algorithms and code interfaces.

	* LWR_fuelcyclecost.py: Module used to calculate fuel costs for LWR loading patterns if this is selected as a objective in the optimization.  
	* decorators.py: File containing additional decorators to increase the readability of midas output files.
	* optimizer_tools.py: Module containing general tools required to perform optimizations including: initial solution generators, solution validation checkers, statistics plot generators, and the fitness function. 
	* problem_preparation.py: Module containing problem specific functions for problems requiring special attention. This includes expanding high-symmetry chromosome representations of loading patterns into lower levels of symmetry (octant core chromosome -> full core map), and counting the number of gene instances in high-symmetry chromosomes representations.  
	* termination_criteria.py: Module for storing termination criteria methods.

* midas/codes: Directory containing code interface modules.

	* parcs343.py:  Module that handles PARCSv3.4.3 calculations evaluation and data extraction.
	* parcs342.py:  Module that handles PARCSv3.4.2 calculations evaluation and data extraction. This module is not actively updated.
	* ipwr_lut.py:  Module that handles reading the IPWR look up table and extracting data.
	* polaris624.py: Module that handles Polarisv6.2.4 calculations evaluation and data extraction.
	* trace50p5.py: Module that handles TRACE calculations evaluation and data extraction.

* midas/algorithms: Directory containing optimization algorithm modules.

	* genetic_algorithm.py: Python file that stores all classes and functions for performing Genetic Algorithm optimization.
	* simulated_annealing.py: Python file that stores all classes and functions for performing Simulated Annealing optimization.
   	* parallel_simulated_annealing.py: Python file that stores all classes and functions for performing Parallel Simulated Annealing optimization. Some fidentical functions are shared with simulated_annealing.py.
	* bayesian_optimization.py: Python file that stores all classes and functions for performing bayesian Optimization.

* midas_tools: Directory containing external tools to assist users in performing optimizations.

	* solution_to_chromosome_tool: Directory containing solution to chromosome tool and sample input files for the tool. The tool enables users to convert solutions into their chromosome representations so that the chromosome can be used as an initaila starting point for optimizations. For example, the tool can convert loading pattern maps as they exist in PARCS input files into their chromosome representations.



# Resources
MIDAS is an updated version of the MOF (Modular Optimization Framework) for which you can find more information about the framework structure, theory and applications in https://doi.org/10.48550/arXiv.2204.00141.

MIDAS previous version MOF has been used to control the Crud deposition in nuclear reactors (https://www.mdpi.com/2673-4117/3/4/36).
