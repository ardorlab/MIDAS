import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path

# Early LCOE implementation without links to a core loading pattern

def LCOE(core_param, lcoe_param):

    efpd = core_param['EFPD']
    nbatches = core_param['Batches']
    tpower = core_param['Thermal_Power']  # MWth
    efficiency = core_param['Efficiency'] 
    nass = core_param['Fuel_Assemblies'] # int
    fuel_rad = core_param['Fuel_Radius'] # cm
    fuel_height = core_param['Fuel_Height'] # cm
    nfuel = core_param['Fuel_Rods'] # int
    enr = core_param['Enrichment'] # fraction
    fuel_dens=core_param['Fuel_Density'] # g/cm^3
    epower = tpower*efficiency

    discount_rate = lcoe_param['Discount_Rate']# fraction
    comp_discount_rate = np.log(1+discount_rate)
    natu_price = lcoe_param['Uranium_Ore_Price']  # $/kgU
    conv_price = lcoe_param['Conversion_Price']  # $/kgHM
    enr_price = lcoe_param['Enrichment_Price']  # $/SWU
    fab_price = lcoe_param['Fabrication_Price'] # $/kgHM
    natu_loss = lcoe_param['Uranium_Ore_Loss']  # fraction
    conv_loss = lcoe_param['Conversion_Loss']  # fraction
    enr_loss = lcoe_param['Enrichment_Loss']  # fraction
    fab_loss = lcoe_param['Fabrication_Loss'] # fraction
    natu_enr  = lcoe_param['Enrichment_Feed']
    tail_enr = lcoe_param['Enrichment_Tail']
    storage_price = lcoe_param['Storage_Price'] # $/kgiHM
    disposal_price = lcoe_param['Disposal_Price'] # $/kgiHM

    fuel_mass = (nass*nfuel*fuel_dens*fuel_height*np.pi*fuel_rad**2)/1000

    discharge_bu = tpower*efpd*nbatches/fuel_mass

    # Fabrication cost

    u2uo2 = 270/238
    fb_u = fuel_mass/u2uo2
    fb_cost = fab_price*fb_u

    # Enrichment cost


    enr_u = fb_u*(1+fab_loss)*(1+enr_loss)
    enr_nat = enr_u*(enr-tail_enr)/(natu_enr-tail_enr)
    sv_nat = (2*natu_enr-1)*np.log(natu_enr/(1-natu_enr))
    sv_enr = (2*enr-1)*np.log(enr/(1-enr))
    sv_tail = (2*tail_enr-1)*np.log(tail_enr/(1-tail_enr))
    uenr_diff = enr_nat-enr_u
    enr_swu = (enr_u*sv_enr+uenr_diff*sv_tail-enr_nat*sv_nat)
    enr_cost = enr_price*enr_swu

    # Converion

    u2uf6 = 352/238
    nat_uf6 = enr_nat
    conv_u = nat_uf6/u2uf6
    conv_cost = conv_u*conv_price

    # Raw ore

    u2u3o8 = 842/714
    raw_u = conv_u*(1+conv_loss)
    raw_u3o8=raw_u*u2u3o8
    raw_cost = natu_price*raw_u

    # Disposal cost 
    store_cost = storage_price*fb_u
    disp_cost = disposal_price*fb_u

    # Net present value

    cost_npv = raw_cost/(1+discount_rate)**(-2) + conv_cost/(1+discount_rate)**(-1.5) + enr_cost/(1+discount_rate)**(-1) + fb_cost/(1+discount_rate)**(-0.5) + store_cost/(1+discount_rate)**(5+efpd*nbatches/365.25) + disp_cost/(1+discount_rate)**(efpd*nbatches/365.25)
    qener_npv = epower*8766*(1-np.exp(-comp_discount_rate*efpd*nbatches/365.25))/comp_discount_rate
    lcoe = cost_npv/qener_npv
    return(lcoe,discharge_bu)

lcoe_param={'Discount_Rate': 0.07,
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
            'Disposal_Price': 0.00711}

core_param={'EFPD': 400,
            'Batches': 3,
            'Thermal_Power': 3800,
            'Efficiency': 0.33,
            'Fuel_Assemblies': 193,
            'Fuel_Rods': 264,
            'Fuel_Radius': 0.45,
            'Fuel_Height': 350,
            'Enrichment': 0.045,
            'Fuel_Density': 10.23}

lcoe,bu = LCOE(core_param,lcoe_param)
print(f"The LCOE fuel cycle cost is: {lcoe}")
print(f"The discharge burnup is: {bu}")