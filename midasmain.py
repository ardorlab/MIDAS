#!/usr/bin/env python3

# # # # # # # # #
# Import Block  #
# # # # # # # # #
import os
import sys
import argparse
import logging
import logging.config
from copy import deepcopy
import midas_data
from midas.input_parser import  Input_Parser
from midas.optimizer import Optimizer
from midas.utils.problem_preparation import Prepare_Problem_Values as prep_inp
from midas.utils.decorators import error_handler, timer, profiler
# # # # # # # # # # #
# Define Functions  #
# # # # # # # # # # #
def Parse_Args():
    """
    Function for parsing command line arguments.
    
    Updated by Nicholas Rollins. 10/03/2024
    """
    input_help_message = 'Input file containing all the information'
    input_help_message += 'needed to perform the optimization routine.'  
    input_help_message += '\n Input file extension needs to be a '
    input_help_message += '.yaml file.'

    cpu_help_message = 'Number of cpus wanted for solving the optimization '
    cpu_help_message += 'problem.'

    parser = argparse.ArgumentParser()
    parser.add_argument('--input',help=input_help_message,
                        required=True,type=str)
    parser.add_argument('--cpus',help=cpu_help_message,
                        required=True,type=int)
    parser.add_argument('--test',help="Used to test changes to the optimization algorithms.",
                    required=False,type=bool)
    parser.add_argument('--restart',help="Used to restart the optimization after pauses.",
                   required=False,type=bool)

    return parser.parse_args()

class Formatter(logging.Formatter): #!TODO: It might be cleaner to move this to a separate file.
    """
    Create dynamic formatting of the logger based on message level.
    
    Written by Nicholas Rollins. 10/03/2024
    """
    def format(self, record):
        if record.levelno == logging.INFO:
            self._style._fmt = "%(message)s"
        else:
            self._style._fmt = "%(levelname)s: %(message)s"
        return super().format(record)

@error_handler
@timer
@profiler
def main():
    """
    Primary execution pathway for MIDAS optimization.
    
    Written by Nicholas Rollins. 10/08/2024
    """
    exitcode = 0

## Set global variables and settings
    #!TODO: read global variables from external file.
    #!global midas_rng_seed = 

## Read command line arguments
    args = Parse_Args()
    if not '.yaml' in args.input:
        raise ValueError("Input File needs to be a valid .yaml file")

## Print MIDAS header
    logger.info(midas_data.__logo__)
    logger.info("MIDAS Version %s\n\n", midas_data.__version__)
    logger.info("====================================================================================================\n")

## Parse input file
    inp_lines = Input_Parser(args.cpus, args.input)
## Prepare input values for writing
    inp_lines = prep_inp.prepare_cycle(inp_lines)
    logger.info("Parsed input file: %s", str(args.input))

## Generate optimizer
    optimizer = Optimizer(inp_lines)
    optimizer.build_optimizer()
    logger.info("Completed Optimizer assembly.")

## Execute
    logger.info("Begin Optimization...\n")
    optimizer.main() #execute optimization through the algorithm class.
    logger.info("\nOptimization completed.")

    return exitcode


# # # # # # # # # # # # # # # #
#  Primary Execution Pathway  #
# # # # # # # # # # # # # # # #
if __name__ == "__main__":
    #Clear output file
    if os.path.exists(midas_data.__ofile__):
        os.remove(midas_data.__ofile__)
    
    #Initialize logging
    logger = logging.getLogger("MIDAS_logger")
    filehandler = logging.FileHandler(filename=midas_data.__ofile__)
    consolehandler = logging.StreamHandler()
    filehandler.setFormatter(Formatter())
    consolehandler.setFormatter(Formatter())
    logger.setLevel(logging.INFO)
    logger.addHandler(filehandler)
    logger.addHandler(consolehandler)
    

    #Execute MIDAS

    exitcode = main()

    
    #Clean up
    logger.info("MIDAS execution completed.")
    logging.shutdown()
    sys.exit(exitcode)