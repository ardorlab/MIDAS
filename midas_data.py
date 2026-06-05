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
'''path for load the base models for prediction without retrainning [RPFmodel, coreparametermodel, CNNpinmodel]'''
__path_base_model__ = [r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/NE512-RPF-model/MIONet_PWR3D_07-50000.ckpt',
                       r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/NE512-coreparam-model/MIONet_PWR3D_core-50000.ckpt', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale_updated_ne512.h5']
''' path for scaler and trunk data [RPFmodel, coreparametermodel]'''
__path_base_data__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_NE512/',
                      '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_coredata_NE512/',
                      ]
'''path for saving new data when perform adaptive retraining [RPFmodel, coreparametermodel]'''
'''Should be the same as __path_base_data__ if one wants to use surrogate models without retraining '''
__path_new_data__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_NE512/',
                      '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/traindataall_coredata_NE512/',
                      ]
''' path for the cross-section library in dictionary format'''
__path_xs_pickle__ = '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/xsdata/updateXSMay26.joblib'
''' path to store raw data in numpy array format for retraining'''
__path_to_store_retrain_data__ = './retraindata/'
''' path to restore which models to retrain from '''
__path_to_restore_model__ = [r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/base_model/PWR-model07/MIONet_PWR3D_07-200000.ckpt',
                       r'/home/khnguy22/Deeponet-midas/MIDAS/surmodel/base_model/PWR-modelcoredata/MIONet_PWR3D_core-200000.ckpt', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale_updated_157.h5']
''' path to save retrained models '''
__path_to_save_model__ = ['/home/khnguy22/Deeponet-midas/MIDAS/surmodel/PWR-model07/MIONet_PWR3D_07', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/PWR-modelcoredata/MIONet_PWR3D_core', 
                       '/home/khnguy22/Deeponet-midas/MIDAS/surmodel/pinmodel/pin_power_unet_deep_noscale_updated_157.h5']
''' deteremine how many epochs you want to perform retrainning, typically 20,000 for a known models,
    50,000 for initial calibration phase in case of unknown problem
'''
__training_epochs_model_1__ = 50000
__training_epochs_model_2__ = 50000
