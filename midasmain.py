# # # # # # # # #
# Import Block  #
# # # # # # # # #
import os
import sys
import argparse
import logging
import logging.config
from copy import deepcopy
from midas import midas_version as version
from midas.input_parser import  Input_Parser
from midas.optimizer import Optimizer


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

def main():
    exitcode = 0
    try:
    #!TODO: set global variables from external file.
    
## Read command line arguments
        args = Parse_Args()
        if not '.yaml' in args.input:
            raise ValueError("Input File needs to be a valid .yaml file")

## Print MIDAS header
        logger.info(version.__logo__)
        logger.info("MIDAS Version %s\n", version.__version__)
    
## Parse input file
        inp_lines = Input_Parser(args.cpus, args.input)
        logger.info("Parsed input file: %s", str(args.input))
    
## Generate optimizer
        optimizer = Optimizer(inp_lines)
        optimizer.build_optimizer()
        logger.info("Completed Optimizer assembly.")
    
## Execute
        logger.info("Begin Optimization...\n")
        optimizer.main() #execute optimization through the algorithm class.
        logger.info("\nOptimization completed.")
    
    except Exception as e:
        exitcode = e.value

    return exitcode


# # # # # # # # # # # # # # # #
#  Primary Execution Pathway  #
# # # # # # # # # # # # # # # #
if __name__ == "__main__":
    #Clear output file
    if os.path.exists(version.__ofile__):
        os.remove(version.__ofile__)
    
    #Initialize logging
    logger = logging.getLogger()
    handler = logging.FileHandler(filename=version.__ofile__)
    handler.setFormatter(Formatter())
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    
    #Execute MIDAS
    exitcode = main()
    
    #Clean up
    logger.info("MIDAS execution completed.")
    logging.shutdown()
    sys.exit(exitcode)