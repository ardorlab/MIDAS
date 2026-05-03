#!/usr/bin/env python3

# Current MIDAS version
__version__ = "0.2.3"


# MIDAS ASCII art
__logo__ = r'''
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

    ofile    :    Name of the output file that will contain progress messages, results, and error messages.

Written by Nicholas Rollins. 10/03/2024
"""

__ofile__ = "midas.out"


# Interface code executable paths
__parcs342exe__ = "/cm1/codes/parcs_342/Executables/Linux/parcs-v342-linux2-intel-x64-release.x"
# __parcs343exe__ = "/cm/shared/nuclearCodes/parcs-3.4.3/PARCS-v343_Exe/Executables/Linux/parcs-v343-linux2-intel-x64-release.x"
__parcs343exe__ = "/data/cm/shared/nuclearCodes/parcs-3.4.3/PARCS-v343_Exe/Executables/Linux/parcs-v343-linux2-intel-x64-release.x"
__trace50p5exe__ = "/cm1/apps/ncsu/TRACE_PARCS/TRACE-V50P5-Exe/Executables/trace-V50p5-linux2-lahey-x64-release.exe"
__polaris624exe__ = "/cm/shared/codes/scale/SCALE-6.2.4/bin/scalerte"
#### Path for data/ surrogate model loading
__path_base_model__ = [r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/core193rpf-case05/MIONet_PWR3D_07-50000.ckpt',
                       r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/core193coreparam-case05/MIONet_PWR3D_core-200000.ckpt', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/base_model/pinmodel/pin_power_unet_deep_noscale_mod01.h5']
__path_base_data__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_core193/',
                      '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_coredata_core193/',
                      ]
__path_new_data__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall/',
                      '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_coredata/',
                      ]
__path_xs_pickle__ = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/updateXS.joblib'
__path_to_store_retrain_data__ = './retraindata/'
__path_to_restore_model__ = [r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/core193rpf-case05/MIONet_PWR3D_07-20000.ckpt', 
                       r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/core193coreparam-case05/MIONet_PWR3D_core-200000.ckpt', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale_mod01.h5']
__path_to_save_model__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/PWR-model07/MIONet_PWR3D_07', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/PWR-modelcoredata/MIONet_PWR3D_core', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale_mod01.h5']
__training_epochs_model_1__ = 20000
__training_epochs_model_2__ = 20000