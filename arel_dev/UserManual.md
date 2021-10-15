This is the Modular Optimization Framework. It utilizes inheritance, object-oriented, and functional programming to
create a simple, robust tool for solving nuclear engineering design problems. This file will explain each of the tools
available in the code, as well as how to run the code and set up the input yaml file.

# Install the code on the RDFMG cluster

The first thing to do is install Miniconda. This will allow you to create a controlled Python environment where you can
install the required packages, because you don't have permission to install them on the cluster by default. Go to the 
site: https://docs.conda.io/en/latest/miniconda.html and download the latest 64-bit Python 3 installer for Linux. Once 
the file has downloaded, upload it to the cluster.

Now enter the following commands:

    bash file_i_just_uploaded.sh

    pip install yaml

    pip install numpy

    pip install matplotlib

    git clone https://github.ncsu.edu/bdander3/modularOptimizationFramework.git

Congratulations. The code is now installed in your local section of the RDFMG cluster.
  
# Running the code. 

There are two ways to run the code on the RDFMG cluster. Interactively and through a batch script.

Enter the modularOptimizationFramework directory. 

To run interactively, type:

    module load workbench

    isqp

    python mofMain.py --input ne_512_project.yaml --cpus 4 or bash submission_script.sh

A batch job has already been written for the example input. To run it, enter the command:

    sbatch submission_script.sh
    
Running the code is as sample as that. mofMain is the main function of the python code obviously. There are two 
commands here though. The first is --input. This command is used to designate the input yaml file that you want to run.
In the case we just ran, we tried to minimize the enrichment of a PWR lattice in octant symmetry. The second command is
--cpus. This command is used to designate how many processors you want to apply to the optimization problem. In the 
case we just ran, four.

# Input Files

So now that you have run the code, let's look at what settings we used. The settings for the optimization are input 
through a yaml file. YAML is a high level markup file. YAML stands for Yeah, Aint Markup Language. It's a very simple,
easy to read and write file that can be directly loaded into python as a dictionary data tool, making it extremely easy
to use and perfect for our purposes here.

Open up the example.yaml file and take a look at it. Obviously use the file editor of your preference. If you don't 
know any file editors or are looking for one, I highly recommed VSCode. 

There is generally a wide range of flexibility in the YAML files from input to input, which can make it complicated to
use the optimization program. However, there are several things that will always be consistent. This section will 
detail the general outline of the yaml files to give you a brief exposure, then all options available will be described
 in the last section.

Yaml files may be written in any order. Like Python, indentation is used to seperate different layers in the file. 
Here, the markers in the example input file are explained:

    optimization: 
    The top marker in the input file. Designates that the cards underneath are related to the optimization settings.

        methodology: genetic_algorithm
        The methodology marker designates which type of optimization is being performed, in this case genetic algorithm.

        population_size: 8
        Indicates population size used in the genetic algorithm.

        number_of_generations: 5
        Indicates number of generations over which optimization will run. 

        mutation:
        Marker indicating sub markers are related to how mutation is performed.
            
            method: [mutate_by_common,mutate_fixed]
            Specifies how mutation is performed. In this case two mutation methods have been selected.

            initial_rate: 0.25
            The starting percentage of population that undergo mutation.

            final_rate: 0.75
            The final percentage of population that undergo mutation.

        fixed_problem: True
        Indicates that the genome is held fixed in some way.

        number_fixed_groups: 4
        Indicates that there are four genome groups used to fix the optimization problem.

        fixed_genes_per_group: [17,14,16,17]
        How many genes must fit into each group in the solution. Note that the order of these groups is the order
        that the gene groups appear when the genes are specified.

    selection:
    Indicates that the sub markers will be related to how selection is performed using the genetic algorithm. 
    
        fitness: ranked
        Marker for selecting which fitness/scoring function to use to compare solutions. In this case the fitness is 
        derived from how the solutions are ranked from best to worst. 

        method: tournament
        Indicates solutions are going to be compared using a tournament method. 

    data_type: loading_pattern
    data_type designates the type of problem that is to be optimized. IN this case the fuel loading pattern of a 
    nuclear reactor.

    objectives:
        Markers under this marker indicate what objectives will be taken into account in the fitness functions. 

        assembly_power:
        Example of an optimization objective. In this case assembly radial peaking factors.
            goal: minimize
            Says that radial peaking factors should be minimized in objective function.
            
genome:
  Markers under this marker are used to describe the actual problem that is to be solved.

  chromosomes:
  Markers under this marker describe the genes used in the optimization. 

    Assembly_One: 
    The first gene in the example problem. Genes are directly under the chromosome marker, and may use any name.

      gene_group: 2.0
      We specified earlier that our optimization problem is fixed. The gene_group is used to identify common genes that
      are held fixed, i.e., labels that this falls into the first group that is only allowed to have 17 assemblies of 
      this label. gene_groups can have any name, but the order in which they appear corresponds to the numbering in the 
      fixed_genes_per_group marker.

      type: 2
      For simulate problems, the type marker corresponds to an assembly designator used in simulate. 

      name: 2.0_w/o
      The Name card is used for solution plotting and for designating casmo lattice types. If you are working on the NE
      412/512 project, you don't need to worry about these cards and they can have any name you would like. 

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
      The map marker specifies where genes may be expressed in the optimization, e.g. locations in the reactor core 
      where assemblies may be placed. Binary markers are used to indicate whether a gene may or may not be expressed in
      the location. A 1 indicates the gene may be expressed there, a 0 indicates it may not be expressed there.

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
  This marker is used to attach additional information required to run simulate. 

    type: pwr
    The reactor type to be ran in simulate.

    core_width: 15
    The size of the reactor core. If there are 152 assemblies the core width is 15.

    load_point: 0.000
    The point that you want the restart file to load from. 

    depletion: 20
    The max depletion step allowed in the simulate calculations. Note that just because you designate a maximum 
    depletion doesn't mean you reach that depletion time step. 

    batch_number: 0
    The cycle of the core. For initial loadings of the reactor core you can specify either 0 or 1.

    pressure: 2250.
    The operating pressure of the reactor core.

    boron: 900.
    A guess of the initial critical boron calculation. Need a guess. Doesn't matter what. 

    power: 100.
    The percent of rated power that the reactor is operating at.

    flow: 100.
    The percent of rated flow that the reactor is operating at. 

    inlet_temperature: 550.
    The inlet temperature of the coolant.

    restart_file: s3.pwr.uo2.c02.depl.res
    The restart file being used in the simulate analysis.

    cs_library: pwr.sim_one.lib
    The cms_link cross section library being used in the simulate analysis.


For more information on every single option available to you in the Modular Optimization Framework, please consult the
wiki. 
