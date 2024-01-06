#!/usr/bin/env python
# coding: utf-8

# In[1]:
import alpha_finder
import setup
from alpha_finder import CocultureAlphaFinder
from setup import load_models
import pandas as pd
import multiprocessing
import json
import concurrent.futures
from functools import partial
import cometspy as c

potential_genes = alpha_finder.potential_genes
# potential_genes = ['folP']

alpha_table = pd.read_csv('./Data/alpha_table_m1.csv', index_col=0) 
initial_pop = 1e-8 

load_model_kwargs = {  
                'gal_scale': 3,
                'ac_scale': 10, 
                'w': True}

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 180)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

print(p.all_params['maxCycles'])

E0, S0, all_components = load_models(Egal_reuptake=False, **load_model_kwargs)

obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'

# obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'
file_suffix = 'm2'
# file_suffix = 'Jt2'

# file_suffix = 'acnb'

# start with gr_Normal 
def get_gr_Normal():     
    # return 0.11688402350837822 # Jan5 scale ac10 result
    
    AF = CocultureAlphaFinder(model=[E0, S0], search_alpha=None, current_gene = 'Normal', p=p,
                             exp_leap= 1.3, alpha_table=alpha_table, carbon_source_val=.1,
                            target_obj_val=.5, precision=1, initial_pop=initial_pop, obj_style=obj_style)
    gr_Normal = AF.calculate_gr_Normal()
    print(f'---------GRNROMAL {gr_Normal}')
    return gr_Normal
    
gr_Normal = get_gr_Normal()
# potential_genes = ['serC']
# potential_genes = ['dapF', 'acnB', 'pgk', 'pfkA', 'mrdA', 'serC', 'dadX', 'gltD', 'gapA']
# potential_genes = ['folP','folA','acnB', 'mrdA']
# potential_genes = ['folP', 'folA']
# search_alpha = 13
search_alpha = None

kwargs = {
    'exp_leap': 1.5, # for coculture gr inherently > 1, require large leap 
    'alpha_table': alpha_table,
    'p': p, 
    'carbon_source_val': .1,
    'target_obj_val': .5,
    'precision': 2, 
    'gr_Normal': gr_Normal,
    'initial_pop': initial_pop, 
    'obj_style': obj_style
}
 
opt_alpha_list, trace_dict = alpha_finder.run_coculture_search_mp(potential_genes, file_suffix, n_processor=10,
                                                                model=[E0, S0], search_alpha=search_alpha,
                                                                **kwargs)

