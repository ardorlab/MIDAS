import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path
sys.path.append('/home/gkdelipe/codes/mof/MOF/')
from lcoe import LCOE, LCOE2

# Manual test case for demonstartion of how LCOE can be computed

cycle_param={'EFPD': 388,
                'Batches': 3,
                'Thermal_Power': 3800,
                'Efficiency': 0.33,
                'Fuel_Assemblies': 193}

lcoe_param={'Cycle': 'Once_Through',
            'Discount_Rate': 0.07,
            'Uranium_Ore_Price': 80,
            'Conversion_Price': 10,
            'Enrichment_Price': 160,
            'Fabrication_Price': 250,
            'Uranium_Ore_Loss': 0.002,
            'Conversion_Loss': 0.002,
            'Enrichment_Loss': 0.002,
            'Fabrication_Loss': 0.002,
            'Enrichment_Feed': 0.00711,
            'Enrichment_Tail': 0.003,
            'Storage_Price': 200,
            'Disposal_Price': 463,
            'Uranium_Ore_Time': -2.0,
            'Conversion_Time': -1.5,
            'Enrichment_Time': -1.0,
            'Fabrication_Time': -0.5,
            'Storage_Time': 5.0+cycle_param['EFPD']*cycle_param['Batches']/365.25,
            'Disposal_Time': cycle_param['EFPD']*cycle_param['Batches']/365.25}

asb_fe200 = {'Number': 61,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.02,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 250,
            }

asb_fe250bp = {'Number': 32,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.025,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 260,
            }

asb_fe250 = {'Number': 24,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.025,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 250,
            }


asb_fe320bp = {'Number': 32,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.032,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 260,
            }

asb_fe320 = {'Number': 44,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.032,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 250,
            }

asb_fe450 = {'Number': 193,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.045,
            'Fuel_Density': 10.23,
            'Fabrication_Price': 250,
            }



asb_param=  {'FE200': asb_fe200,
            'FE250Bp': asb_fe250bp,
            'FE250': asb_fe250,
            'FE320Bp': asb_fe320bp,
            'FE320': asb_fe320}

#asb_param=  {'FE450': asb_fe450}

lcoe,bu, asb_cost = LCOE(cycle_param,lcoe_param, asb_param)
print(f"The LCOE once-through fuel cycle cost is: {lcoe} $/MWh")
print(f"The discharge burnup is: {bu} MWd/kgU")

lcoe2_param=copy.deepcopy(lcoe_param)
lcoe2_param['Reprocess_Price'] = 1600
lcoe2_param['Rep_Waste_Price'] = 185
lcoe2_param['RepU_fraction'] = 0.93
lcoe2_param['RepU_sell_price'] = 108.4
lcoe2_param['Pu_fraction'] = 0.011
lcoe2_param['RepU_buy_price'] = 10
lcoe2_param['MOX_fabrication'] = 2400
lcoe2_param['MOX_RepU_fraction'] = 0.913
lcoe2_param['Reprocess_Time'] = 5.0+cycle_param['EFPD']*cycle_param['Batches']/365.25
lcoe2_param['Rep_Waste_Time'] = 5.0+cycle_param['EFPD']*cycle_param['Batches']/365.25
lcoe2_param['RepU_Time'] = 1.0 + 5.0 + cycle_param['EFPD']*cycle_param['Batches']/365.25

lcoe,lcoe1, lcoe2, pu_price, bu = LCOE2(cycle_param,lcoe2_param, asb_param)
print(f"The LCOE twice-through fuel cycle cost is: {lcoe} $/MWh")
print(f"The discharge burnup is: {bu} MWd/kgU")