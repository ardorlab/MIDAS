The  Modular Optimization Framework (MOF) is written in Python 3. It utilizes both object oriented and functional programming. Class inheritance is also used extensively to reduce the amount of physical code.
If you are developing for this program, there are several concepts that you should be familiar with before working on the code. While a broad knowledge of Python will be helpful, specific topics that you should be familiar with before reading and working on the code include: classes, dictionaries, static methods, and parallel programming. 

General Python formatting should always be followed when editing the code. This document is to explain a few simple rules that should be followed when adding to MOF. This way, the code will remain clear and concise throughout. 

1. The code follows the underscore format. I don't know if it has a defined name. Variables should be formatted using underscores to seperate words. Proper variable names include Simulate_Class, some_unknown_value, specific_action_function. Do not format variables like GeneticAlgorithmClass, populationcount, theGOgetterFUNCTION.

2. This should be included in best coding practices, but a special note of it will be made here too. Class names should be nouns, and are always capitalized. Function/method names should general be verbs and or actions. They are not capitalized. Example class names are: Genetic_Algorithm, Optimization_Factory, Casmo_Solution. Example function names are return_peak_powers, mate_solutions, crossover.

3. Try to use as verbose of variable names as possible. The previously given examples of variable names are all good. they tell you exactly what the variable is for. This isn't Fortran. Don't use variable names like cmnt, ct, etc. Feel free to write things out.

4. In seemingly direct contrast to the previous rule, do not use lines longer than 120 characters. The 80 character variable limit is best to adhere to, but sometimes that isn't always possible. Lines should never be longer than 120 characters though.

5. Classes, functions, and methods should always have a docstring explaining what the purpose of it is. Below is an example of a function and the accompanying docstring. 

    def illustrate_docstring():
      
      """
      
      This function illustrates how to use a docstring to explain the purpose of the function.
    
      Parameters: None
    
      Written by Brian Andersen. 1/7/2019
      
      """
      
6. There should always be an error message associated with any user specified input that will cause the code to crash when the user input is incorrect or missing. For example      

       if 'gene_group' in chrom_settings[chrom]:
           stuff_happens...
       else:
           error = "Fixed Genome Group Problem has been specified, however gene "
           error += "{} does not have a specified gene_group.".format(chrom)
           raise KeyError(error)
           
7. Be cognizant of requesting user specified input. Do you want the user to have to specify it exactly as you have? Do they need to get    capitalization perfect as well? Using commands such as .upper() or .lower() are good ways to make it easier for users to use your        code and not have to fix as many typos. 
