#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import cometspy as c
import numpy as np
import itertools
import ast
import logging
from time import time
from itertools import chain  
import concurrent.futures
from setup import (load_models, 
                    get_BM_df, 
                    rename_columns, 
                    gene_index_culture_col_df)

no_monoculture=2
initial_pop = 1e-8 # initial biomass in comets simulation
n_combos = 50 # first n gene combo 
n_processor = 17 # number of process run in parallel

E0, S0, all_components = load_models()

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 450)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

# alpha table keep E0, S0 column only
'''
alpha_table = pd.read_csv('./Data/alpha_table.csv', index_col='Gene_inhibition')
# alpha_table = pd.read_csv('./Data/alpha_table_wgal.csv', index_col='Gene_inhibition')

alpha_table.columns = ['E0', '_', 'S0']
alpha_table = pd.read_csv('./Data/checker_alpha_table.csv', index_col = 0)
'''

file_list = ['alpha_table_m1', 'alpha_table_m2', 'alpha_table_m3']
# file_list = ['alpha_table_m1', 'alpha_table_m3']

# def generate_csv(gene_combos: list, filename: str, alpha_table, mono=True):
def generate_csv(gene_combos: list, filename: str, **kwargs):
    
    # Multiprocessing double gene
    result_list = list()
     
    with concurrent.futures.ProcessPoolExecutor(n_processor) as executor:
        for i, current_gene in enumerate(gene_combos):
            future = executor.submit(get_BM_df, current_gene=current_gene, n_dir=i, **kwargs)
            result_list.append(future)
    
    result_df_list, result_dict_list = zip(*[r.result() for r in result_list])
    
    # with multiprocessing.Pool(n_processor) as pool:
    #     result_df_list, result_dict_list = zip(*pool.starmap(get_BM_df, zipped_arg)) 
    result_df = pd.concat(result_df_list,axis = 1)
    result_df.columns = rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename}.csv')
    
    # TODO: bypass DataFrame construction, write json directly
    analysis_df = gene_index_culture_col_df(pd.DataFrame(result_dict_list))
    analysis_df.to_json(f"./Data/flanalysis_{filename}.json") 

    return result_df, analysis_df 

# Run
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Trial")
    
    start = time()
    
    gene_combos = pd.read_csv('./Data/GeneCombos.csv',header=None)[0][:n_combos]
    gene_combos = list(ast.literal_eval(ele) for ele in gene_combos)
    
    kwargs = {'alpha_table': None, 
                'mono': True, 
                'p': p,
                'E0': E0, 
                'S0': S0,
                'return_sim': False, 
                'ko': False}
    
    for filename in file_list: 
        alpha_table = pd.read_csv(f'./Data/{filename}.csv', index_col=0) 
        kwargs['alpha_table'] = alpha_table
        gene50 = list(alpha_table.index)
        gene50.extend(['Normal'])
        
        if '_m2' in filename:
            alpha_table.iloc[:, :2] = alpha_table.iloc[:, :2].where(alpha_table.iloc[:, :2] <= 5e3, 1e5) # reset to 1e5 for nonessential genes
        
        method_n = filename.split('_')[-1]
        SGresult_df, SGanalysis_df  = generate_csv(gene50, f'BM_SG_{method_n}', **kwargs)
        DGresult_df, DGanalysis_df  =  generate_csv(gene_combos, f'BM_DG_{method_n}', **kwargs) 
    
    end = time() 
    print('Time Elapsed: ', end-start)

    # test    
    # _ = generate_csv(gene_combos[:1], 't1') 
    # _ = generate_csv(gene50[:1], 't2')


# load json
# analysis_df = pd.read_json('./Data/flanalysis.json')
# pd.read_json(analysis_df.coculture_media[0]).query("metabolite =='lcts_e'")
# # .plot(x='cycle', y=['conc_mmol'])
# analysis_df.plot(x='cycle', y = ['ACtex','EX_ac_e'])

