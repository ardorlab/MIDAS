import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path
import yaml 
sys.path.append('/home/gkdelipe/codes/mof/MOF/')
from parcs_332 import Loading_Pattern_Solution
from lcoe import LCOE, LCOE2

# LCOE calculation for specific core loading pattern
# Both Once-Through and Twice-Through cycles are supported.

parameters = {'cycle_length': {'value': 0.0, 'weight':1.0},
                'PinPowerPeaking': {'value': 0.0, 'weight':400.0},
                'FDeltaH': {'value': 0.0, 'weight':400.0},
                'max_boron': {'value': 0.0, 'weight':-1.0}}

core_param = {   'xs_library':'/home/gkdelipe/codes/mof/xslib',
                'power': 3800.0,
                'flow': 18231.89,
                'inlet_temperature': 565.0,
                'number_assemblies':193,
                'symmetry': 'octant',
                'map_size': "quarter",
                'symmetry_axes': ((8,8),(16,16),(16,8))}

parcs_lp = Loading_Pattern_Solution()
settings={'optimization':{'reproducer':'test'},
          'genome': {"parcs_data":core_param }}


core_inventory={'FE200': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'200', 'Cross_Section':'xs_g200_gd_0_wt_0'},
                'FE220': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'220', 'Cross_Section':'xs_g220_gd_0_wt_0'},
                'FE250': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'250', 'Cross_Section':'xs_g250_gd_0_wt_0'},
                'FE280': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'280', 'Cross_Section':'xs_g280_gd_0_wt_0'},
                'FE300': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'300', 'Cross_Section':'xs_g300_gd_0_wt_0'},
                'FE320': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'320', 'Cross_Section':'xs_g320_gd_0_wt_0'},
                'FE350': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'350', 'Cross_Section':'xs_g350_gd_0_wt_0'},
                'FE400': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'400', 'Cross_Section':'xs_g400_gd_0_wt_0'},
                'FE400GD4WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'401', 'Cross_Section':'xs_g400_gd_4_wt_2'},
                'FE400GD8WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'402', 'Cross_Section':'xs_g400_gd_8_wt_2'},
                'FE400GD12WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'403', 'Cross_Section':'xs_g400_gd_12_wt_2'},
                'FE400GD16WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'404', 'Cross_Section':'xs_g400_gd_16_wt_2'},
                'FE400GD20WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'405', 'Cross_Section':'xs_g400_gd_20_wt_2'},
                'FE400GD24WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'406', 'Cross_Section':'xs_g400_gd_24_wt_2'},
                'FE400GD4WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'407', 'Cross_Section':'xs_g400_gd_4_wt_4'},
                'FE400GD8WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'408', 'Cross_Section':'xs_g400_gd_8_wt_4'},
                'FE400GD12WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'409', 'Cross_Section':'xs_g400_gd_12_wt_4'},
                'FE400GD16WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'410', 'Cross_Section':'xs_g400_gd_16_wt_4'},
                'FE400GD20WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'411', 'Cross_Section':'xs_g400_gd_20_wt_4'},
                'FE400GD24WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'412', 'Cross_Section':'xs_g400_gd_24_wt_4'},
                'FE400GD4WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'413', 'Cross_Section':'xs_g400_gd_4_wt_6'},
                'FE400GD8WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'414', 'Cross_Section':'xs_g400_gd_8_wt_6'},
                'FE400GD12WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'415', 'Cross_Section':'xs_g400_gd_12_wt_6'},
                'FE400GD16WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'416', 'Cross_Section':'xs_g400_gd_16_wt_6'},
                'FE400GD20WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'417', 'Cross_Section':'xs_g400_gd_20_wt_6'},
                'FE400GD24WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'418', 'Cross_Section':'xs_g400_gd_24_wt_6'},
                'FE400GD4WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'419', 'Cross_Section':'xs_g400_gd_4_wt_8'},
                'FE400GD8WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'420', 'Cross_Section':'xs_g400_gd_8_wt_8'},
                'FE400GD12WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'421', 'Cross_Section':'xs_g400_gd_12_wt_8'},
                'FE400GD16WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'422', 'Cross_Section':'xs_g400_gd_16_wt_8'},
                'FE400GD20WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'423', 'Cross_Section':'xs_g400_gd_20_wt_8'},
                'FE400GD24WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'424', 'Cross_Section':'xs_g400_gd_24_wt_8'},
                'FE450': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'450', 'Cross_Section':'xs_g450_gd_0_wt_0'},
                'FE450GD4WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'451', 'Cross_Section':'xs_g450_gd_4_wt_2'},
                'FE450GD8WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'452', 'Cross_Section':'xs_g450_gd_8_wt_2'},
                'FE450GD12WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'453', 'Cross_Section':'xs_g450_gd_12_wt_2'},
                'FE450GD16WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'454', 'Cross_Section':'xs_g450_gd_16_wt_2'},
                'FE450GD20WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'455', 'Cross_Section':'xs_g450_gd_20_wt_2'},
                'FE450GD24WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'456', 'Cross_Section':'xs_g450_gd_24_wt_2'},
                'FE450GD4WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'457', 'Cross_Section':'xs_g450_gd_4_wt_4'},
                'FE450GD8WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'458', 'Cross_Section':'xs_g450_gd_8_wt_4'},
                'FE450GD12WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'459', 'Cross_Section':'xs_g450_gd_12_wt_4'},
                'FE450GD16WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'460', 'Cross_Section':'xs_g450_gd_16_wt_4'},
                'FE450GD20WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'461', 'Cross_Section':'xs_g450_gd_20_wt_4'},
                'FE450GD24WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'462', 'Cross_Section':'xs_g450_gd_24_wt_4'},
                'FE450GD4WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'463', 'Cross_Section':'xs_g450_gd_4_wt_6'},
                'FE450GD8WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'464', 'Cross_Section':'xs_g450_gd_8_wt_6'},
                'FE450GD12WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'465', 'Cross_Section':'xs_g450_gd_12_wt_6'},
                'FE450GD16WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'466', 'Cross_Section':'xs_g450_gd_16_wt_6'},
                'FE450GD20WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'467', 'Cross_Section':'xs_g450_gd_20_wt_6'},
                'FE450GD24WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'468', 'Cross_Section':'xs_g450_gd_24_wt_6'},
                'FE450GD4WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'469', 'Cross_Section':'xs_g450_gd_4_wt_8'},
                'FE450GD8WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'470', 'Cross_Section':'xs_g450_gd_8_wt_8'},
                'FE450GD12WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'471', 'Cross_Section':'xs_g450_gd_12_wt_8'},
                'FE450GD16WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'472', 'Cross_Section':'xs_g450_gd_16_wt_8'},
                'FE450GD20WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'473', 'Cross_Section':'xs_g450_gd_20_wt_8'},
                'FE450GD24WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'474', 'Cross_Section':'xs_g450_gd_24_wt_8'},
                'FE495': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'500', 'Cross_Section':'xs_g495_gd_0_wt_0'},
                'FE495GD4WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'501', 'Cross_Section':'xs_g495_gd_4_wt_2'},
                'FE495GD8WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'502', 'Cross_Section':'xs_g495_gd_8_wt_2'},
                'FE495GD12WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'503', 'Cross_Section':'xs_g495_gd_12_wt_2'},
                'FE495GD16WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'504', 'Cross_Section':'xs_g495_gd_16_wt_2'},
                'FE495GD20WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'505', 'Cross_Section':'xs_g495_gd_20_wt_2'},
                'FE495GD24WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'506', 'Cross_Section':'xs_g495_gd_24_wt_2'},
                'FE495GD4WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'507', 'Cross_Section':'xs_g495_gd_4_wt_4'},
                'FE495GD8WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'508', 'Cross_Section':'xs_g495_gd_8_wt_4'},
                'FE495GD12WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'509', 'Cross_Section':'xs_g495_gd_12_wt_4'},
                'FE495GD16WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'510', 'Cross_Section':'xs_g495_gd_16_wt_4'},
                'FE495GD20WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'511', 'Cross_Section':'xs_g495_gd_20_wt_4'},
                'FE495GD24WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'512', 'Cross_Section':'xs_g495_gd_24_wt_4'},
                'FE495GD4WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'513', 'Cross_Section':'xs_g495_gd_4_wt_6'},
                'FE495GD8WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'514', 'Cross_Section':'xs_g495_gd_8_wt_6'},
                'FE495GD12WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'515', 'Cross_Section':'xs_g495_gd_12_wt_6'},
                'FE495GD16WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'516', 'Cross_Section':'xs_g495_gd_16_wt_6'},
                'FE495GD20WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'517', 'Cross_Section':'xs_g495_gd_20_wt_6'},
                'FE495GD24WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'518', 'Cross_Section':'xs_g495_gd_24_wt_6'},
                'FE495GD4WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'519', 'Cross_Section':'xs_g495_gd_4_wt_8'},
                'FE495GD8WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'520', 'Cross_Section':'xs_g495_gd_8_wt_8'},
                'FE495GD12WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'521', 'Cross_Section':'xs_g495_gd_12_wt_8'},
                'FE495GD16WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'522', 'Cross_Section':'xs_g495_gd_16_wt_8'},
                'FE495GD20WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'523', 'Cross_Section':'xs_g495_gd_20_wt_8'},
                'FE495GD24WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'524', 'Cross_Section':'xs_g495_gd_24_wt_8'},
                'FE550': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'550', 'Cross_Section':'xs_g550_gd_0_wt_0'},
                'FE550GD4WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'551', 'Cross_Section':'xs_g550_gd_4_wt_2'},
                'FE550GD8WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'552', 'Cross_Section':'xs_g550_gd_8_wt_2'},
                'FE550GD12WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'553', 'Cross_Section':'xs_g550_gd_12_wt_2'},
                'FE550GD16WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'554', 'Cross_Section':'xs_g550_gd_16_wt_2'},
                'FE550GD20WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'555', 'Cross_Section':'xs_g550_gd_20_wt_2'},
                'FE550GD24WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'556', 'Cross_Section':'xs_g550_gd_24_wt_2'},
                'FE550GD4WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'557', 'Cross_Section':'xs_g550_gd_4_wt_4'},
                'FE550GD8WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'558', 'Cross_Section':'xs_g550_gd_8_wt_4'},
                'FE550GD12WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'559', 'Cross_Section':'xs_g550_gd_12_wt_4'},
                'FE550GD16WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'560', 'Cross_Section':'xs_g550_gd_16_wt_4'},
                'FE550GD20WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'561', 'Cross_Section':'xs_g550_gd_20_wt_4'},
                'FE550GD24WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'562', 'Cross_Section':'xs_g550_gd_24_wt_4'},
                'FE550GD4WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'563', 'Cross_Section':'xs_g550_gd_4_wt_6'},
                'FE550GD8WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'564', 'Cross_Section':'xs_g550_gd_8_wt_6'},
                'FE550GD12WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'565', 'Cross_Section':'xs_g550_gd_12_wt_6'},
                'FE550GD16WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'566', 'Cross_Section':'xs_g550_gd_16_wt_6'},
                'FE550GD20WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'567', 'Cross_Section':'xs_g550_gd_20_wt_6'},
                'FE550GD24WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'568', 'Cross_Section':'xs_g550_gd_24_wt_6'},
                'FE550GD4WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'569', 'Cross_Section':'xs_g550_gd_4_wt_8'},
                'FE550GD8WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'570', 'Cross_Section':'xs_g550_gd_8_wt_8'},
                'FE550GD12WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'571', 'Cross_Section':'xs_g550_gd_12_wt_8'},
                'FE550GD16WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'572', 'Cross_Section':'xs_g550_gd_16_wt_8'},
                'FE550GD20WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'573', 'Cross_Section':'xs_g550_gd_20_wt_8'},
                'FE550GD24WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'574', 'Cross_Section':'xs_g550_gd_24_wt_8'},
                'FE600': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'600', 'Cross_Section':'xs_g600_gd_0_wt_0'},
                'FE600GD4WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'601', 'Cross_Section':'xs_g600_gd_4_wt_2'},
                'FE600GD8WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'602', 'Cross_Section':'xs_g600_gd_8_wt_2'},
                'FE600GD12WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'603', 'Cross_Section':'xs_g600_gd_12_wt_2'},
                'FE600GD16WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'604', 'Cross_Section':'xs_g600_gd_16_wt_2'},
                'FE600GD20WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'605', 'Cross_Section':'xs_g600_gd_20_wt_2'},
                'FE600GD24WT2': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'606', 'Cross_Section':'xs_g600_gd_24_wt_2'},
                'FE600GD4WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'607', 'Cross_Section':'xs_g600_gd_4_wt_4'},
                'FE600GD8WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'608', 'Cross_Section':'xs_g600_gd_8_wt_4'},
                'FE600GD12WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'609', 'Cross_Section':'xs_g600_gd_12_wt_4'},
                'FE600GD16WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'610', 'Cross_Section':'xs_g600_gd_16_wt_4'},
                'FE600GD20WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'611', 'Cross_Section':'xs_g600_gd_20_wt_4'},
                'FE600GD24WT4': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'612', 'Cross_Section':'xs_g600_gd_24_wt_4'},
                'FE600GD4WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'613', 'Cross_Section':'xs_g600_gd_4_wt_6'},
                'FE600GD8WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'614', 'Cross_Section':'xs_g600_gd_8_wt_6'},
                'FE600GD12WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'615', 'Cross_Section':'xs_g600_gd_12_wt_6'},
                'FE600GD16WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'616', 'Cross_Section':'xs_g600_gd_16_wt_6'},
                'FE600GD20WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'617', 'Cross_Section':'xs_g600_gd_20_wt_6'},
                'FE600GD24WT6': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'618', 'Cross_Section':'xs_g600_gd_24_wt_6'},
                'FE600GD4WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'619', 'Cross_Section':'xs_g600_gd_4_wt_8'},
                'FE600GD8WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'620', 'Cross_Section':'xs_g600_gd_8_wt_8'},
                'FE600GD12WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'621', 'Cross_Section':'xs_g600_gd_12_wt_8'},
                'FE600GD16WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'622', 'Cross_Section':'xs_g600_gd_16_wt_8'},
                'FE600GD20WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'623', 'Cross_Section':'xs_g600_gd_20_wt_8'},
                'FE600GD24WT8': {'Max_Limit':np.inf, 'In_Design':0, 'Cost':0, 'Tag':'624', 'Cross_Section':'xs_g600_gd_24_wt_8'}
        }

core_inventory_groups={'E200': {'Values':['FE200'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E220': {'Values':['FE220'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E250': {'Values':['FE250'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E280': {'Values':['FE280'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E300': {'Values':['FE300'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E320': {'Values':['FE320'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E350': {'Values':['FE350'], 'Limit': 'Max', 'Limit_Value':np.inf},
                        'E400': {'Values':['FE400',
                                            'FE400GD4WT2','FE400GD8WT2','FE400GD12WT2','FE400GD16WT2','FE400GD20WT2','FE400GD24WT2',
                                            'FE400GD4WT4','FE400GD8WT4','FE400GD12WT4','FE400GD16WT4','FE400GD20WT4','FE400GD24WT4',
                                            'FE400GD4WT6','FE400GD8WT6','FE400GD12WT6','FE400GD16WT6','FE400GD20WT6','FE400GD24WT6',
                                            'FE400GD4WT8','FE400GD8WT8','FE400GD12WT8','FE400GD16WT8','FE400GD20WT8','FE400GD24WT8'],
                                'Limit': 'Max', 'Limit_Value':np.inf},
                        'E450': {'Values':['FE450',
                                            'FE450GD4WT2','FE450GD8WT2','FE450GD12WT2','FE450GD16WT2','FE450GD20WT2','FE450GD24WT2',
                                            'FE450GD4WT4','FE450GD8WT4','FE450GD12WT4','FE450GD16WT4','FE450GD20WT4','FE450GD24WT4',
                                            'FE450GD4WT6','FE450GD8WT6','FE450GD12WT6','FE450GD16WT6','FE450GD20WT6','FE450GD24WT6',
                                            'FE450GD4WT8','FE450GD8WT8','FE450GD12WT8','FE450GD16WT8','FE450GD20WT8','FE450GD24WT8'],
                                'Limit': 'Max', 'Limit_Value':np.inf},
                        'E495': {'Values':['FE495',
                                            'FE495GD4WT2','FE495GD8WT2','FE495GD12WT2','FE495GD16WT2','FE495GD20WT2','FE495GD24WT2',
                                            'FE495GD4WT4','FE495GD8WT4','FE495GD12WT4','FE495GD16WT4','FE495GD20WT4','FE495GD24WT4',
                                            'FE495GD4WT6','FE495GD8WT6','FE495GD12WT6','FE495GD16WT6','FE495GD20WT6','FE495GD24WT6',
                                            'FE495GD4WT8','FE495GD8WT8','FE495GD12WT8','FE495GD16WT8','FE495GD20WT8','FE495GD24WT8'],
                                'Limit': 'Max', 'Limit_Value':np.inf},
                        'E550': {'Values':['FE550',
                                            'FE550GD4WT2','FE550GD8WT2','FE550GD12WT2','FE550GD16WT2','FE550GD20WT2','FE550GD24WT2',
                                            'FE550GD4WT4','FE550GD8WT4','FE550GD12WT4','FE550GD16WT4','FE550GD20WT4','FE550GD24WT4',
                                            'FE550GD4WT6','FE550GD8WT6','FE550GD12WT6','FE550GD16WT6','FE550GD20WT6','FE550GD24WT6',
                                            'FE550GD4WT8','FE550GD8WT8','FE550GD12WT8','FE550GD16WT8','FE550GD20WT8','FE550GD24WT8'],
                                'Limit': 'Max', 'Limit_Value':np.inf},
                        'E600': {'Values':['FE600',
                                            'FE600GD4WT2','FE600GD8WT2','FE600GD12WT2','FE600GD16WT2','FE600GD20WT2','FE600GD24WT2',
                                            'FE600GD4WT4','FE600GD8WT4','FE600GD12WT4','FE600GD16WT4','FE600GD20WT4','FE600GD24WT4',
                                            'FE600GD4WT6','FE600GD8WT6','FE600GD12WT6','FE600GD16WT6','FE600GD20WT6','FE600GD24WT6',
                                            'FE600GD4WT8','FE600GD8WT8','FE600GD12WT8','FE600GD16WT8','FE600GD20WT8','FE600GD24WT8'],
                                'Limit': 'Max', 'Limit_Value':np.inf}}

settings['genome']['inventory']=core_inventory
settings['genome']['inventory_groups']=core_inventory_groups

parcs_lp.add_additional_information(settings)
parcs_lp.name = 'run'
parcs_lp.parameters=parameters
state0_file = 'state0.yml'
with open(state0_file) as f:
    state0 = yaml.safe_load(f)

parcs_lp.genome=list(state0.values())
start = time.time()
parcs_lp.evaluate()
end = time.time()
print('Total Parcs Evaluation and Processing Time: {} s'.format(end-start))
lcoe = parcs_lp.additional_parameters["LCOE"]
bu = parcs_lp.additional_parameters['Discharge_Burnup']
print(f"The LCOE once-through fuel cycle cost is: {lcoe} $/MWh")
print(f"The discharge burnup is: {bu} MWd/kgU")



