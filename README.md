

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

* samples: Directory including various sample cases for users to get familiarized with MIDAS framework. Users just need to navigate to the sample folder and run MIDAS from there with the corresponding input file or use the run.sh bash file on the RDFMG cluster. The samples consist of examples of genetic algorithm and simulated annealing optimization problems for a reduced number of code calculations. All the samples should finish within 10 minutes. Users are encourage to change the number of code evaluations (e.g. generations number) to experiment with optimization algorithms. Parallel execution is implemented only for genetic algorithm and the user could increase the number of allocated processors to evaluate the impact on the execution time.
  * sample_0: First Cycle Core Loading Pattern Optimization with PARCS. Random selection is used. 
  * sample_1: Fuel Lattice Optimization with NCSU lattice simulator. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_2: First Cycle Core Loading Pattern Optimization with NCSU core simulator. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_3: First Cycle Core Loading Pattern Optimization with NCSU core simulator. Simulated annealing is used with 20 iterations. 
  * sample_4: Third Cycle Core Loading Pattern Optimization with NCSU core simulator. Genetic algorithm is used with 2 generations and 10 population per generation. 
  * sample_5: Third Cycle Core Loading Pattern Optimization with NCSU core simulator. Simulated annealing is used with 20 iterations.
  * sample_6: First Cycle Core Loading Pattern Optimization with PARCS. Genetic Algorithm is used. 
  * sample_7: First Cycle Core Loading Pattern Optimization with PARCS. Simulate Annealing is used. 
  * sample_8: First Cycle Core Loading Pattern Optimization with NCSU core simulator. Parallel Simulate Annealing is used. 
  * sample_9: First Cycle Core Loading Pattern Optimization with PARCS. Parallel Simulate Annealing is used.
  * sample_10: First Cycle Core Loading Pattern Optimization with PARCS. Reinforcement Learning is used. 
  * dev_samples: Sample cases under development and not fully functional yet.

* documentation: Directory including all additional documentation.

* ncsu_lattice.py: Python file that handles NCSU lattice calculations evaluation and data extraction.

* crudworks.py: Python file used if CRUD machine learning predictions are required (TensorFlow should be installed).

* fitness.py: Python file for selecting and computing the objective function of the optimization. 

* geneticAlgorithm.py: Python file that stores all classes and functions for performing Genetic Algorithm optimization.

* metrics.py: Python file including tools for tracking solutions and storing the generated optimization data.

* mofMain.py: Python file that is the main body of MIDAS. In this file the interface between the input file and the optimization is performed by selecting the specified options and initializing all the necessary components.

* ncsu_core.py: Python file that handles NCSU core simulator calculations evaluation and data extraction.

* parcs.py: Python file that handles PARCS calculations evaluation and data extraction.

* parcs_332.py: Python file that handles PARCS_332 calculations evaluation and data extraction.

* randomSolutions.py: Python file that stores all classes and functions for performing optimization with random solutions.

* lcoe.py: computation of levelized cost of electricity.

* simulateAnnealing.py: Python file that stores all classes and functions for performing Simulate Annealing optimization.

* reinforcement_learning.py: Python file that stores all classes and functions for performing Reinforcement Learning optimization.

* solution_types.py: Python file for storing the solutions of the optimization together with some usefull functions.

* submission_script.sh: Bash file example for running MIDAS through SLURM on the RDFMG cluster.

# Resources
MIDAS is an updated version of the MOF (Modular Optimization Framework) for which you can find more information about the framework structure, theory and applications in https://doi.org/10.48550/arXiv.2204.00141.

MIDAS previous version MOF has been used to control the Crud deposition in nuclear reactors (https://www.mdpi.com/2673-4117/3/4/36).
