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
from itertools import product

no_monoculture=2
initial_pop = 1e-6 # initial biomass in comets simulation
n_combos = 50 # first n gene combo 
n_processor = 10 # number of process run in parallel
carbon_source_val = .05
mono = False
mono_S = True
co = True
test = False

# mono = True
mono_S = False
# co = False
test = True

E0, S0, all_components = load_models(Egal_reuptake=False)

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
# p.set_param("maxCycles", 450)
p.set_param("maxCycles", 80)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)
# p.set_param('FluxLogRate', 1)
# p.set_param('MediaLogRate', 1)

obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'
# obj_style = 'MAXIMIZE_OBJECTIVE_FLUX'


# alpha table keep E0, S0 column only
'''
alpha_table = pd.read_csv('./Data/alpha_table.csv', index_col='Gene_inhibition')
# alpha_table = pd.read_csv('./Data/alpha_table_wgal.csv', index_col='Gene_inhibition')

alpha_table.columns = ['E0', '_', 'S0']
alpha_table = pd.read_csv('./Data/checker_alpha_table.csv', index_col = 0)
'''

file_list = ['alpha_table_m1', 'alpha_table_m2', 'alpha_table_m3']
file_list = ['alpha_table_m2', 'alpha_table_m3']
file_list = ['alpha_table_m1']

# def generate_csv(gene_combos: list, filename: str, alpha_table, mono=True):
def generate_csv(gene_combos: list, filename: str, test=False, force_list = [0], **kwargs):
    
    test_suffix = '_test' if test else ''
    # Multiprocessing double gene
    result_list = list()
     
     
    with concurrent.futures.ProcessPoolExecutor(n_processor) as executor:
        for i, scale_force_gal_c in enumerate(force_list):
            future = executor.submit(get_BM_df, current_gene=gene_combos, n_dir=i,scale_force_gal_c = scale_force_gal_c, **kwargs)
            result_list.append(future)
    
    result_df_list, result_dict_list = zip(*[r.result() for r in result_list])
    
    # with multiprocessing.Pool(n_processor) as pool:
    #     result_df_list, result_dict_list = zip(*pool.starmap(get_BM_df, zipped_arg)) 
    result_df = pd.concat(result_df_list,axis = 1)
    result_df.columns = rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename+test_suffix}.csv')
    
    # TODO: bypass DataFrame construction, write json directly
    analysis_df = gene_index_culture_col_df(pd.DataFrame(result_dict_list))
    analysis_df.to_json(f"./Data/flanalysis_{filename+test_suffix}.json") 
    return result_df, analysis_df 

def SG_diff_alpha(alpha_table): # single genes with different alpha compared to alpha_table_m1(mostly nonessential)
    alpha_table_m1 = pd.read_csv(f'./Data/alpha_table_m1.csv', index_col=0)
    TF_df = alpha_table_m1 == alpha_table
    TF_df['product'] = TF_df.apply(lambda x: x['E0'] * x['S0'], axis=1)
    return list(TF_df.query('product == False').index), set(TF_df.query('product == True').index)

def fill_skipped(method_list=['m2','m3'], XG_list=['SG','DG'], test=True):
    def fill_skipped_cols(df_need_filled, df_m1):
        overlap_cols = df_m1.columns.difference(df_need_filled.columns)
        return df_need_filled.merge(df_m1[overlap_cols], left_index=True, right_index=True)
    test_suffix = '_test' if test else ''
    
    Biomass_m1 = {}
    for XG, method in product(XG_list, method_list):
        if method != 'm1':
            m1_str, mx_str = f'BM_{XG}_m1', f'BM_{XG}_{method}'
            if m1_str not in Biomass_m1:
                Biomass_m1[m1_str] = pd.read_csv(f'./Data/{m1_str}.csv', index_col=0)
            Biomass_mx = pd.read_csv(f'./Data/{mx_str+test_suffix}.csv', index_col=0)
            Biomass_mx = fill_skipped_cols(Biomass_mx, Biomass_m1[m1_str])
            Biomass_mx.to_csv(f'./Data/{mx_str}.csv')

def DG_diff_alpha(DG, repeated_SG):
    print('rp', repeated_SG)
    print('DDG', [gcomb for gcomb in DG if len(set(gcomb) & repeated_SG)!=2])
    return [gcomb for gcomb in DG if len(set(gcomb) & repeated_SG)!=2]
    
# Run
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Trial")
    
    start = time()
    
    gene_combos = pd.read_csv('./Data/GeneCombos.csv',header=None)[0][:n_combos]
    DG = list(ast.literal_eval(ele) for ele in gene_combos)
    kwargs = {'alpha_table': None, 
                'mono': mono, 
                'mono_S': mono_S, 
                'co' : co, 
                'p': p,
                'E0': E0, 
                'S0': S0,
                'return_sim': False, 
                'ko': False,
                'obj_style' : obj_style,
                'carbon_source_val': carbon_source_val
                
    }
    filename = 'alpha_table_m1'
    method_n = filename.split('_')[-1]
    # kwargs['alpha_table'] = alpha_table

    SG = ['Normal']
    # modify scaling at every E_model 

    # SGresult_df, SGanalysis_df  = generate_csv(SG, f'BM_SG_{method_n}', test=test, force_list = np.linspace(0, 1, 4), **kwargs)
    SGresult_df, SGanalysis_df  = generate_csv(SG, f'BM_SG_{method_n}', test=test, force_list = [0,.1,.5,1,3], **kwargs)
    
    end = time() 
    print('Time Elapsed: ', end-start)

    # TODO: merge repeated SGDG from m1 dfs 


# load json
# analysis_df = pd.read_json('./Data/flanalysis.json')
# pd.read_json(analysis_df.coculture_media[0]).query("metabolite =='lcts_e'")
# # .plot(x='cycle', y=['conc_mmol'])
# analysis_df.plot(x='cycle', y = ['ACtex','EX_ac_e'])

