import numpy as np
from midas.utils.problem_preparation import LWR_Core_Shapes
from midas.utils import optimizer_tools as optools
"""
    Method for calculating the fuel cycle cost of an LWR loading pattern.
    If the chromosome is formatted as a shuffling scheme, the full-core loading pattern is first extracted.
    
    Written by Nicholas Rollins. 12/11/2024
"""


def get_fuelcycle_cost(soln, input):
    """
    Prepares and calls the necessary calculations for evaluating the LWR fuel cycle cost.
    
    Written by Nicholas Rollins. 12/11/2024
    """
    ## fetch the duplication multiplicity of each location when expanded to the full core.
    multdict = LWR_Core_Shapes.get_symmetry_multiplicity(input.nrow, input.ncol, input.num_assemblies, input.symmetry)
    
    ## Interpret loading pattern from chromosome
    if input.calculation_type in ['eq_cycle']:
        loading_pattern = optools.Constrain_Input.SS_decoder(soln.chromosome)
    else:
        loading_pattern = soln.chromosome

    ## calculate fuel cost for each FA (accounting for solution symmetry)
    fuelcycle_cost_total = 0.0 # USD
    for i in range(len(loading_pattern)):
        FA_type = loading_pattern[i]
    
        x_p = input.fa_options['fuel'][FA_type]['enrichment'] # w.t.
        M_p = input.fa_options['fuel'][FA_type]['hm_loading'] # kg
        fuelcycle_cost_total += calc_fuelcost(x_p,M_p)*multdict[i] # USD
        
        try:
            if input.fa_options['fuel'][FA_type]['blanket']:
                blanket_type = input.fa_options['fuel'][FA_type]['blanket']
                x_p = input.fa_options['blanket'][blanket_type]['enrichment'] # w.t.
                M_p = input.fa_options['fuel'][blanket_type]['hm_loading'] # kg
                fuelcycle_cost_total += calc_fuelcost(x_p,M_p)*multdict[i] # USD
        except:
            pass # FA doesn't include blanket, continue.

    soln.parameters['cost_fuelcycle']['value'] = fuelcycle_cost_total # USD
    
    return soln.parameters


def calc_fuelcost(x_p, mass_product):
    """
    x_p: weight fraction (not percent) of U235 in the enrichment product stream.
    mass_product: initial heavy metal loading of the enrichment product stream.

    Written by Nicholas Rollins. 12/11/2024
    """
    x_f = 0.00711 # w.t.
    
    x_t = 0.00300 # w.t.
    
    a = 2.6 # lb_U3O8/kg_U; conversion factor of U feed
    
    t = 1.0 # years
    t_b = 1.0 # years
    
    S_feed = 3.0 # %/year
    S_conv = 3.0 # %/year
    S_enr = 3.0 # %/year
    S_fab = 3.0 # %/year
    
    l_feed = 0.0 # w.t.
    l_conv = 0.005 # w.t.
    l_enr = 0.0005 # w.t.
    l_fab = 0.0 # w.t.
    
    f_feed = (1 + l_conv)*(1 + l_fab)
    f_conv = (1 + l_fab)
    f_enr = (1 + l_fab)
    f_fab = 1.0
    
    P_feed = 31.35 # USD/lb_U3O8
    P_conv = 5.00 # USD/kg_U3O8
    P_enr = 105 # USD/SWU
    P_fab = 210 # USD/kg_U3O8
    
    
    mass_feed = mass_product * ((x_p - x_t)/(x_f - x_t)) # kg
    mass_tails = mass_feed - mass_product # kg
    
    V_product = (2*x_p - 1) * np.log(x_p/(1 - x_p))
    V_tails = (2*x_t - 1) * np.log(x_t/(1 - x_t))
    V_feed = (2*x_f - 1) * np.log(x_f/(1 - x_f))
    
    SWU = mass_product * V_product + mass_tails * V_tails - mass_feed * V_feed
    
    cost_feed = mass_feed * a * f_feed * P_feed * (1 + S_feed)**(t - t_b) # USD
    cost_conversion = mass_feed * f_conv * P_conv * (1 + S_conv)**(t - t_b) # USD
    cost_enrichment = SWU * f_enr * P_enr * (1 + S_enr)**(t - t_b) # USD
    cost_fabrication = mass_product * P_fab * (1 + S_fab)**(t - t_b) # USD
    
    return cost_feed + cost_conversion + cost_enrichment + cost_fabrication # USD