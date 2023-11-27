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

alpha_table = pd.read_csv('./Data/alpha_table_csearch.csv', index_col=0) 
initial_pop = 1e-8
initial_pop = 1e-6
E0, S0, all_components = load_models()

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 100)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

# obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'
obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'
file_suffix = 'MAXBM'

# start with gr_Normal 
def get_gr_Normal():    
    # return .000523
    AF = CocultureAlphaFinder(model=[E0, S0], search_alpha=None, current_gene = 'Normal', p=p,
                             exp_leap= 3, alpha_table=alpha_table, carbon_source_val=.1,
                            target_obj_val=.5, precision=1, initial_pop=initial_pop, obj_style=obj_style)
    gr_Normal = AF.calculate_gr_Normal()
    print(f'---------GRNROMAL {gr_Normal}')
    return gr_Normal
    
gr_Normal = get_gr_Normal()
# potential_genes = ['gltA']
# search_alpha = 3.2
search_alpha = None


opt_alpha_list, trace_dict = alpha_finder.run_coculture_search_mp(potential_genes, 'm2MAXBM', n_processor=10,
                                                                model=[E0, S0], search_alpha=search_alpha,
                                                                exp_leap= 3, alpha_table=alpha_table, p=p,
                                                                carbon_source_val=.1, target_obj_val=.5, 
                                                                precision=1, gr_Normal=gr_Normal, initial_pop=initial_pop, obj_style=obj_style)

