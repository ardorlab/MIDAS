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
import random
import pickle
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
    input_help_message = 'Input file containing all the information needed to perform the optimization routine.\n'
    input_help_message += ' Input file extension needs to be a .yaml file.'

    cpu_help_message = 'Number of cpus wanted for solving the optimization problem.'

    parser = argparse.ArgumentParser()
    parser.add_argument('--input',help=input_help_message,
                        required=False,type=str,default=None)
    parser.add_argument('--cpus',help=cpu_help_message,
                        required=False,type=int,default=1)
    #!parser.add_argument('--test',help="Used to test changes to the optimization algorithms.",
    #!                    required=False,type=bool,default=None) #!TODO: this option is not in use.
    parser.add_argument('--restart',help="A valid '.rst' MIDAS restart file name, used to continue a previous optimization routine.",
                        required=False,type=str,default=None)

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
def restart(args):
    """
    Secondary execution pathway for MIDAS optimization allowing the 
    calculation to continue from a previously completed iteration.
    
    Written by Nicholas Rollins. 11/26/2024
    """
    exitcode = 0
    
## Print MIDAS header
    logger.info(midas_data.__logo__)
    logger.info("MIDAS Version %s\n\n", midas_data.__version__)
    logger.info("====================================================================================================\n")
    
## Open the restart file in binary read mode
    logger.info("Reading restart file %s...\n", args.restart)
    with open(args.restart, "rb") as f:
        # Load the object from the file
        optimizer = pickle.load(f)
    
## Seed the global RNG
    optimizer.input.set_seed = random.randrange(sys.maxsize) # generate an artibtrary RNG seed
    random.seed(optimizer.input.set_seed)
    logger.info("Set global RNG Seed: %s", optimizer.input.set_seed)    
    
## Execute
    logger.info("Restart Optimization...\n")
    optimizer.main() #execute optimization through the algorithm class.
    logger.info("\nOptimization completed.")
    
    return exitcode

@error_handler
@timer
@profiler
def main(args):
    """
    Primary execution pathway for MIDAS optimization.
    
    Written by Nicholas Rollins. 10/08/2024
    """
    exitcode = 0

## Print MIDAS header
    logger.info(midas_data.__logo__)
    logger.info("MIDAS Version %s\n\n", midas_data.__version__)
    logger.info("====================================================================================================\n")

## Parse input file
    inp_lines = Input_Parser(args.cpus, args.input)
## Prepare input values for writing
    inp_lines = prep_inp.prepare_cycle(inp_lines)
    logger.info("Parsed input file: %s", str(args.input))
    
## Seed the global RNG
    if not inp_lines.set_seed:
        inp_lines.set_seed = random.randrange(sys.maxsize) # generate an artibtrary RNG seed
    random.seed(inp_lines.set_seed)
    logger.info("Set global RNG Seed: %s", inp_lines.set_seed)

## Update logger level
    if inp_lines.debug_mode:
        logger.setLevel(logging.DEBUG)

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
    
    #Read command line arguments
    args = Parse_Args()
    if args.input and not '.yaml' in args.input:
        raise ValueError("Input File needs to be a valid '.yaml' file.")
    if args.restart and not '.rst' in args.restart:
        raise ValueError("Restart File needs to be a valid '.rst' file.")
    if not args.input and not args.restart: #no execution type specified
        raise NameError("One of the following arguments is required: '--input', '--restart'.")
    

    if not args.restart: # Execute MIDAS optimization
        exitcode = main(args)
    
    else: # perform restart
        exitcode = restart(args)

    
    #Clean up
    logger.info("MIDAS execution completed.")
    logging.shutdown()
    sys.exit(exitcode)