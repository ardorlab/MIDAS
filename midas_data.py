#!/usr/bin/env python3

# Current MIDAS version
__version__ = "0.2.2"


# MIDAS ASCII art
__logo__ = '''
 .----------------.  .----------------.  .----------------.  .----------------.  .----------------.   
| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |  
| | ____    ____ | || |     _____    | || |  ________    | || |      __      | || |    _______   | |  
| ||_   \  /   _|| || |    |_   _|   | || | |_   ___ `.  | || |     /  \     | || |   /  ___  |  | |  
| |  |   \/   |  | || |      | |     | || |   | |   `. \ | || |    / /\ \    | || |  |  (__ \_|  | |  
| |  | |\  /| |  | || |      | |     | || |   | |    | | | || |   / ____ \   | || |   '.___`-.   | |  
| | _| |_\/_| |_ | || |     _| |_    | || |  _| |___.' / | || | _/ /    \ \_ | || |  |`\____) |  | |  
| ||_____||_____|| || |    |_____|   | || | |________.'  | || ||____|  |____|| || |  |_______.'  | |  
| |              | || |              | || |              | || |              | || |              | |  
| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |  
 '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  
'''


# # # # # # # # # # #
# Global Variables  #
# # # # # # # # # # #
"""
The below values may be edited by the User as needed.

	ofile	:	Name of the output file that will contain progress messages, results, and error messages.

Written by Nicholas Rollins. 10/03/2024
"""

__ofile__ = "midas.out"


# Interface code executable paths
__parcs342exe__ = "/cm1/codes/parcs_342/Executables/Linux/parcs-v342-linux2-intel-x64-release.x"
__parcs343exe__ = "/cm/shared/nuclearCodes/parcs-3.4.3/PARCS-v343_Exe/Executables/Linux/parcs-v343-linux2-intel-x64-debug.x"