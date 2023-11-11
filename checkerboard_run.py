import Div_COMETS as DC
import multiprocessing
import pandas as pd
import json

n_processor = 10

alpha_table = pd.read_csv('./Data/checker_alpha_table.csv', index_col = 0)


def checkerboard_simulation(gene_combo: tuple, alpha_table, filename='BM_checkerboard', mono=True):
    def convert_checker_alpha_table_to_list(checker_alpha_table):
        l = list()
        for lv_pairs in checker_alpha_table.lv_pairs.unique():                                      
            l.append(checker_alpha_table.query('lv_pairs == @lv_pairs'))
        return l
    
    sub_alpha_list = convert_checker_alpha_table_to_list(alpha_table)[:2]
    
    zipped_arg = [[('folP', 'folA'), i%n_processor, sub_alpha_df, mono, '_'+str(i)] for i, sub_alpha_df in enumerate(sub_alpha_list)] # also zip alpha_table, 'monoculture' to pass to get_BM_df for discrete concentration 

    # Multiprocessing double gene
    with multiprocessing.Pool(n_processor) as pool:
        result_df_list, result_dict_list = zip(*pool.starmap(DC.get_BM_df, zipped_arg)) 
    result_df = pd.concat(result_df_list,axis = 1)
    result_df.columns = DC.rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename}.csv')
    
    analysis_df = DC.gene_index_culture_col_df(pd.DataFrame(result_dict_list))
    analysis_df.to_json(f"./Data/flanalysis_{filename}.json") 
    return None

checkerboard_gene_combo = ('folP', 'folA')
checkerboard_simulation(checkerboard_gene_combo, alpha_table, 'test_checker', mono=False)

# checkerboard_simulation(checkerboard_gene_combo, DC.alpha_table, 'test_checker', mono=False)
