import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import shutil
import copy
from pathlib import Path
import yaml 

# LCOE calculation for specific core loading pattern
# Both Once-Through and Twice-Through cycles are supported.

def get_asb_mass(asb_dict):
    """
    Calculation of the assembly uo2 mass based on 
    geometric quantitities

    Parameters:
        asb_dict: Dictionary containing the 
        necessary assembly information for the computation.
    """
    nfuel = asb_dict["Fuel_Rods"] # fraction
    fuel_dens = asb_dict["Fuel_Density"] # g/cm^3
    fuel_height = asb_dict["Fuel_Height"] # cm
    fuel_rad = asb_dict["Fuel_Radius"] # cm
    res = (nfuel*fuel_dens*fuel_height*np.pi*fuel_rad**2)/1000
    return(res)


def get_fab_cost(uo2_mass,price):
    """
    Calculation of the UO2 fabrication cost.

    Parameters:
        uo2_mass: Float - mass of UO2 fuel.
        price: Float - fabrication price in $/kgHM
    """
    u2uo2 = 270/238
    u_mass = uo2_mass/u2uo2
    cost = price*u_mass
    return(cost,u_mass)

def get_enr_cost(u_mass,price,nat_enr,enr,tail_enr,fab_loss,enr_loss):
    """
    Calculation of the UF6 enrichment cost.

    Parameters:
        u_mass: Float - mass of enriched U fuel to be produced.
        price: Float - enrichment price in $/SWU.
        nat_enr: Float - feed enrichement.
        enr: Float - produced enrichment.
        tail_enr: Float - tail assay enrichment.
        fab_loss: Float - fabrication loss.
        enr_loss: Float - enrichment loss.
    """
    u2uf6 = 352/238
    enr_u_mass = u_mass*(1+fab_loss)*(1+enr_loss)
    enr_uf6_mass = enr_u_mass*u2uf6
    enr_nat = enr_uf6_mass*(enr-tail_enr)/(nat_enr-tail_enr)
    sv_nat = (2*nat_enr-1)*np.log(nat_enr/(1-nat_enr))
    sv_enr = (2*enr-1)*np.log(enr/(1-enr))
    sv_tail = (2*tail_enr-1)*np.log(tail_enr/(1-tail_enr))
    uenr_diff = enr_nat-enr_uf6_mass
    swu = (enr_uf6_mass*sv_enr+uenr_diff*sv_tail-enr_nat*sv_nat)
    cost = price*swu
    return(cost,enr_nat)


def get_conv_cost(u_mass,price):
    """
    Calculation of the UF6 conversion cost
    from U3O8 to UF6.

    Parameters:
        u_mass: Float - mass of U feed needed for the enrichment process.
        price: Float - conversion price in $/kgHM
    """
    u2uf6 = 352/238
    nat_uf6 = u_mass
    conv_u_mass = nat_uf6/u2uf6
    cost = conv_u_mass*price
    return(cost,conv_u_mass)

def get_raw_cost(u_mass,price,conv_loss):
    """
    Calculation of the raw U purchase cost cost.

    Parameters:
        u_mass: Float - mass of U in the UF6 after conversion.
        price: Float - purchase price in $/kgU
    """
    u2u3o8 = 842/714
    raw_u_mass = u_mass*(1+conv_loss)
    raw_u3o8_mass=raw_u_mass*u2u3o8
    cost = price*raw_u_mass
    return(cost,raw_u3o8_mass)

def get_cost(u_mass,price):
    """
    Generic cost calculator usually for back-end costs.

    Parameters:
        u_mass: Float - mass of U in UO2.
        price: Float -  price in $/kgHM
    """
    cost = u_mass*price
    return(cost)


def get_asb_cost(asb_mass,asb_dict,lcoe_dict,efpd,nbatches,no_back=False):
    """
    Calculate the cost of a fuel assembly.

    Parameters:
        asb_mass: Float - mass of UO2 in the assembly.
        asb_dict: Dictionary - Dictionary with assembly related information.
        lcoe_dict: Dictionary - Dictionary with costs related information.
        efpd: Float - cycle length in equivalent full power days.
        nbatches: Float - number of batches in the core. 
    """

    natu_price = lcoe_dict['Uranium_Ore_Price']  # $/kgU
    conv_price = lcoe_dict['Conversion_Price']  # $/kgHM
    enr_price = lcoe_dict['Enrichment_Price']  # $/SWU
    fab_price = asb_dict['Fabrication_Price'] # $/kgHM
    natu_loss = lcoe_dict['Uranium_Ore_Loss']  # fraction
    conv_loss = lcoe_dict['Conversion_Loss']  # fraction
    enr_loss = lcoe_dict['Enrichment_Loss']  # fraction
    fab_loss = lcoe_dict['Fabrication_Loss'] # fraction
    natu_enr  = lcoe_dict['Enrichment_Feed']
    tail_enr = lcoe_dict['Enrichment_Tail']
    enr = asb_dict['Enrichment'] # fraction
    storage_price = lcoe_dict['Storage_Price'] # $/kgiHM
    disposal_price = lcoe_dict['Disposal_Price'] # $/kgiHM
    nat_time = lcoe_dict['Uranium_Ore_Time'] 
    conv_time = lcoe_dict['Conversion_Time'] 
    enr_time = lcoe_dict['Enrichment_Time'] 
    fab_time = lcoe_dict['Fabrication_Time'] 
    store_time = lcoe_dict['Storage_Time'] 
    disp_time = lcoe_dict['Disposal_Time'] 
    discount_rate = lcoe_dict['Discount_Rate']# fraction

    # Fabrication cost
    fb_cost, u_mass = get_fab_cost(asb_mass,fab_price)
    # Enrichment cost
    enr_cost, enr_nat_mass= get_enr_cost(u_mass,enr_price,natu_enr,enr,tail_enr,fab_loss,enr_loss)
    # Converion
    conv_cost, conv_u_mass = get_conv_cost(enr_nat_mass,conv_price)
    # Raw ore
    raw_cost, raw_u3o8_mass = get_raw_cost(conv_u_mass,natu_price,conv_loss)
    # Disposal cost 
    if no_back:
        cost_npv = raw_cost/(1+discount_rate)**(nat_time) + conv_cost/(1+discount_rate)**(conv_time) + enr_cost/(1+discount_rate)**(enr_time) \
         + fb_cost/(1+discount_rate)**(fab_time)
    else:
        store_cost = get_cost(u_mass,storage_price)
        disp_cost = get_cost(u_mass,disposal_price)
        # Net present value
        cost_npv = raw_cost/(1+discount_rate)**(nat_time) + conv_cost/(1+discount_rate)**(conv_time) + enr_cost/(1+discount_rate)**(enr_time) \
            + fb_cost/(1+discount_rate)**(fab_time) + store_cost/(1+discount_rate)**(store_time) \
            + disp_cost/(1+discount_rate)**(disp_time)
    return(cost_npv)

def LCOE(core_param, lcoe_param, asb_param):
    """
    Calculation of the levelized cost of electricity (LCOE) for the core
    in a Once-Through cycle. This LCOE does not include capital and O&M costs.

    Parameters:
        core_param: Dictionary - Dictionary with core related information.
        lcoe_param: Dictionary - Dictionary with costs related information.
        asb_param: Dictionary - Dictionary with assembly related information.
    """

    efpd = core_param['EFPD']
    nbatches = core_param['Batches']
    tpower = core_param['Thermal_Power']  # MWth
    efficiency = core_param['Efficiency'] 
    nass = core_param['Fuel_Assemblies'] # int
    epower = tpower*efficiency
    
    discount_rate = lcoe_param['Discount_Rate']# fraction
    comp_discount_rate = np.log(1+discount_rate)

    nasb = 0
    asb_cost = []
    asb_mass = []
    asb_mult =  []
    for key, value in asb_param.items():
        nasb+=value['Number']
        asb_mult.append(value['Number'])
        mass = get_asb_mass(value)
        asb_mass.append(mass)
        cost = get_asb_cost(mass,value,lcoe_param,efpd,nbatches)
        asb_cost.append(cost)

    if nasb ==nass:
        pass
    else:
        raise ValueError("Wrong number of assemblies")

    asb_mult = np.array(asb_mult)
    asb_mass = np.array(asb_mass)
    asb_cost = np.array(asb_cost)
   
    fuel_mass = np.dot(asb_mass,asb_mult)
    fuel_cost_npv = np.dot(asb_cost,asb_mult)

    discharge_bu = tpower*efpd*nbatches/fuel_mass

    qener_npv = epower*8766*(1-np.exp(-comp_discount_rate*efpd*nbatches/365.25))/comp_discount_rate
    lcoe = fuel_cost_npv/qener_npv
    return(lcoe,discharge_bu,asb_cost)

def LCOE_MCYC(core_param, lcoe_param, asb_param):
    """
    Calculation of the levelized cost of electricity (LCOE) for the core
    in a Once-Through cycle. Adaptation for multi-cycle calculation.
    This LCOE does not include capital and O&M costs.

    Parameters:
        core_param: Dictionary - Dictionary with core related information.
        lcoe_param: Dictionary - Dictionary with costs related information.
        asb_param: Dictionary - Dictionary with assembly related information.
    """

    efpd = core_param['EFPD']
    nbatches = core_param['Batches']
    tpower = core_param['Thermal_Power']  # MWth
    efficiency = core_param['Efficiency'] 
    nass = core_param['Fuel_Assemblies'] # int
    epower = tpower*efficiency
    
    discount_rate = lcoe_param['Discount_Rate']# fraction
    comp_discount_rate = np.log(1+discount_rate)

    nasb = 0
    asb_cost = []
    asb_mass = []
    asb_mult =  []
    for key, value in asb_param.items():
        nasb+=value['Number']
        asb_mult.append(value['Number'])
        mass = get_asb_mass(value)
        asb_mass.append(mass)
        cost = get_asb_cost(mass,value,lcoe_param,efpd,nbatches)
        asb_cost.append(cost)

    if nasb ==nass:
        pass
    else:
        raise ValueError("Wrong number of assemblies")

    asb_mult = np.array(asb_mult)
    asb_mass = np.array(asb_mass)
    asb_cost = np.array(asb_cost)
   
    fuel_mass = np.dot(asb_mass,asb_mult)
    fuel_cost_npv = np.dot(asb_cost,asb_mult)

    discharge_bu = tpower*efpd*nbatches/fuel_mass

    qener_npv = epower*8766*(1-np.exp(-comp_discount_rate*efpd*nbatches/365.25))/comp_discount_rate
    lcoe = fuel_cost_npv/qener_npv
    return(lcoe,discharge_bu,asb_cost)


def LCOE2(core_param, lcoe_param, asb_param):
    """
    Calculation of the levelized cost of electricity (LCOE) for the core
    in a Twice-Through cycle. This LCOE does not include capital and O&M costs.

    Parameters:
        core_param: Dictionary - Dictionary with core related information.
        lcoe_param: Dictionary - Dictionary with costs related information.
        asb_param: Dictionary - Dictionary with assembly related information.
    """

    efpd = core_param['EFPD']
    nbatches = core_param['Batches']
    tpower = core_param['Thermal_Power']  # MWth
    efficiency = core_param['Efficiency'] 
    nass = core_param['Fuel_Assemblies'] # int
    epower = tpower*efficiency
    
    discount_rate = lcoe_param['Discount_Rate']# fraction
    comp_discount_rate = np.log(1+discount_rate)

    nasb = 0
    asb_cost = []
    asb_mass = []
    asb_mult =  []
    for key, value in asb_param.items():
        nasb+=value['Number']
        asb_mult.append(value['Number'])
        mass = get_asb_mass(value)
        asb_mass.append(mass)
        cost = get_asb_cost(mass,value,lcoe_param,efpd,nbatches,no_back=True)
        asb_cost.append(cost)

    if nasb ==nass:
        pass
    else:
        raise ValueError("Wrong number of assemblies")

    asb_mult = np.array(asb_mult)
    asb_mass = np.array(asb_mass)
    asb_cost = np.array(asb_cost)
   
    fuel_mass = np.dot(asb_mass,asb_mult)
    fuel_u_mass = fuel_mass*238/270
    fuel_cost_npv = np.dot(asb_cost,asb_mult)
    discharge_bu = tpower*efpd*nbatches/fuel_mass
    qener_npv = epower*8766*(1-np.exp(-comp_discount_rate*efpd*nbatches/365.25))/comp_discount_rate
    # fuel front-end cost for reactor1
    f21 = fuel_cost_npv/qener_npv
    # reprocessing cost for reactor1
    proc_cost1 = lcoe_param["Reprocess_Price"]*fuel_u_mass/(1+lcoe_param["Discount_Rate"])**lcoe_param["Reprocess_Time"]
    s21 = proc_cost1/qener_npv
    # wast disposal cost for reactor1
    waste_cost1 = lcoe_param["Rep_Waste_Price"]*fuel_u_mass*((1+discount_rate)**5)/(1+lcoe_param["Discount_Rate"])**lcoe_param["Rep_Waste_Time"]
    w21 = waste_cost1/qener_npv
    # reprocessed uranium credit
    repu_mass = fuel_u_mass*lcoe_param["RepU_fraction"]
    repu_credit = lcoe_param["RepU_sell_price"]*repu_mass/(1+lcoe_param["Discount_Rate"])**lcoe_param["RepU_Time"]
    u21 = repu_credit/qener_npv
    # plutonium reactor1 factor
    pu_mass = fuel_u_mass*lcoe_param["Pu_fraction"]
    pu_factor= pu_mass/(1+lcoe_param["Discount_Rate"])**lcoe_param["RepU_Time"]
    pz21 = pu_factor/qener_npv
    # energy produced per kgHM
    nener_npv = qener_npv/fuel_u_mass
    # purchase of RepU cost for reactor2
    mox_repU_mass_pkg = 1.0*lcoe_param["Fabrication_Loss"]*lcoe_param["MOX_RepU_fraction"]
    pU_cost2 = lcoe_param["Reprocess_Price"]*mox_repU_mass_pkg/(1+lcoe_param["Discount_Rate"])**(-1)
    u22 = pU_cost2/nener_npv
    # plutonium reactor2 factor
    mox_Pu_mass_pkg = 1.0 - mox_repU_mass_pkg
    mox_pu_factor= mox_Pu_mass_pkg/(1+lcoe_param["Discount_Rate"])**(-1)
    pz22 = mox_pu_factor/nener_npv
    # mox fabrication cost for reactor2
    fb_cost2 = lcoe_param["MOX_fabrication"]/(1+lcoe_param["Discount_Rate"])**(-0.5)
    b22 = fb_cost2/nener_npv
    # back-end cost for reactor2
    store_cost2 = lcoe_param["Storage_Price"]/(1+lcoe_param["Discount_Rate"])**lcoe_param["Storage_Time"]
    d22 = store_cost2/nener_npv + 1/0.15
    # Pu price
    Pu_price = (f21 + s21 + w21 - u21 - u22 - b22 -d22)/(pz21+pz22)
    # LCOE calculation for both reactors and the whole cycle
    lcoe = lcoe1 = f21 + s21 + w21 -u21 -pz21*Pu_price
    lcoe2 = u22 + pz22*Pu_price + b22 + d22
    return(lcoe, lcoe1, lcoe2, Pu_price, discharge_bu)