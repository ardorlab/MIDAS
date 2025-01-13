from midas.utils.problem_preparation import LWR_Core_Shapes
from midas.utils import optimizer_tools as optools
from copy import deepcopy
"""
    Methods for calculating the average fuel enrichment of an LWR assembly-based reactor core.
    Supports both loading pattern and shuffling scheme chromosome types.
    
    Written by Nicholas Rollins. 1/13/2025
"""

def get_avfuelenrichment(soln, input):
    """
    Prepares and calls the necessary calculations for evaluating the average LWR fuel enrichment.
    
    Note: This method ignores the enrichment of axial blanket regions.
    
    Written by Nicholas Rollins. 1/13/2025
    """
    ## fetch the duplication multiplicity of each location when expanded to the full core.
    multdict = LWR_Core_Shapes.get_symmetry_multiplicity(input.nrow, input.ncol, input.num_assemblies, input.symmetry)
    
    ## Interpret loading pattern from chromosome
    if input.calculation_type in ['eq_cycle']:
        loading_pattern = optools.Constrain_Input.SS_decoder(soln.chromosome)
        if input.objectives['av_fuelenrichment']['settings'] == 'feed_batch_only':
            # extract the feed batch only
            LP_full = deepcopy(loading_pattern)
            loading_pattern = []
            for i in range(len(LP_full)):
                if type(soln.chromosome[i][1]) == str: #assume fresh fuel from feed batch.
                    loading_pattern.append(LP_full[i])
    else:
        loading_pattern = soln.chromosome

    ## calculate average fuel enrichment (accounting for solution symmetry)
    enrichments_list = []
    for i in range(len(loading_pattern)):
        FA_type = loading_pattern[i]
        
        x_p = input.fa_options['fuel'][FA_type]['enrichment'] # w.t.
        enrichments_list.append(x_p)
        
    soln.parameters['av_fuelenrichment']['value'] = sum(enrichments_list)/len(enrichments_list) # w.t.
    
    return soln.parameters