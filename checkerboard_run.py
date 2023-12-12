from setup import (load_models, 
                    get_BM_df, 
                    rename_columns, 
                    gene_index_culture_col_df)
import multiprocessing
import pandas as pd
import json
from collections import defaultdict
import ast
import cometspy as c
import concurrent.futures

n_processor = 17
alpha_table = pd.read_csv('./Data/checkerboard_alpha_table.csv', index_col = 0)
# alpha_table = pd.read_csv('./Data/subc_alpha_table.csv', index_col = 0)
# first_n = 2
first_n = None
mono=True

E0, S0, all_components = load_models(Egal_reuptake=False)

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
# p.set_param("maxCycles", 30)
p.set_param("maxCycles", 1250)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

def checkerboard_simulation(gene_combo: tuple, alpha_table, filename='BM_checkerboard', mono=True, **kwargs):
    def convert_checker_alpha_table_to_dict(checker_alpha_table): # alpha table keys as append ''.join(['_'+ str(ele) for ele in _i_j as suffix 
        # d = defaultdict(dict)
        d= dict()
        for i, lv_pairs in enumerate(checker_alpha_table.lv_pairs.unique()):   
            d['_'+'.'.join([str(ele) for ele in ast.literal_eval(lv_pairs)])] = checker_alpha_table.query('lv_pairs == @lv_pairs')
            if first_n == i-1:
                return d
        return d
    
    def replace_i_to_suffix(colname, sub_alpha_dict):
        i = colname.split('_')[-1]
        print(i, sub_alpha_dict[i])
        return colname.replace(i, sub_alpha_dict[i]['suffix'])
    
    sub_alpha_dict = convert_checker_alpha_table_to_dict(alpha_table)

    # zipped_arg = [[('folP', 'folA'), i, sub_alpha_table, mono, p, k] # i as directory suffix
    #                 for i, (k, sub_alpha_table) in enumerate(sub_alpha_dict.items())] # also zip alpha_table, 'monoculture' to pass to get_BM_df for discrete concentration 
    
    # TODO: use function from Div_COMETS -difference in enumerate potential_genes VS subalpha_table
    result_list = list()
    with concurrent.futures.ProcessPoolExecutor(n_processor) as executor:
        for i, (checker_suffix, sub_alpha_table) in enumerate(sub_alpha_dict.items()):
            future = executor.submit(get_BM_df, current_gene=gene_combo, n_dir=i,E0=E0, S0=S0,
                                        alpha_table=sub_alpha_table, mono=mono, 
                                        p=p, checker_suffix=checker_suffix, **kwargs)
            result_list.append(future) 
    
    result_df_list, result_dict_list = zip(*[r.result() for r in result_list])
    
    result_df = pd.concat(result_df_list,axis = 1)
    result_df.columns = rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename}.csv')
    
    analysis_df = gene_index_culture_col_df(pd.DataFrame(result_dict_list))
    analysis_df.to_json(f"./Data/flanalysis_{filename}.json") 
    return None

checkerboard_gene_combo = ('folP', 'folA')
# checkerboard_simulation(checkerboard_gene_combo, alpha_table, 'test', mono=mono)
checkerboard_simulation(checkerboard_gene_combo, alpha_table, 'checkerboard_run3', mono=mono)
# checkerboard_simulation(checkerboard_gene_combo, DC.alpha_table, 'test_checker', mono=False)
