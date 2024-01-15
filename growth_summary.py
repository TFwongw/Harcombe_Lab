import numpy as np
import pandas as pd

def get_XG_cycle_from(desired_cycle):
    if len(desired_cycle.index[-1].split('.')) <=2:
        SG_cycle = desired_cycle.loc[[len(ele.split('.')) ==1 for ele in desired_cycle.index]]
        DG_cycle = desired_cycle.loc[[len(ele.split('.')) ==2 for ele in desired_cycle.index]]
    else:
        SG_cycle = desired_cycle.loc[['0' in ele for ele in desired_cycle.index]]
        DG_cycle = desired_cycle.loc[['0' not in ele for ele in desired_cycle.index]]
    SG_cycle.index.name='SG'
    DG_cycle.index.name='DG'
    SG_cycle.columns.name=None
    DG_cycle.columns.name=None
    return SG_cycle, DG_cycle    

def search_gr_cycle_with_biomass(df_search, biomass_values):
    return [df_search[df_search >= biomass_value].idxmin() 
                for biomass_value in list(biomass_values)]

def get_maximum_growth_cycle(desired_BM):
    c_max_gr = desired_BM.iloc[1]+ (desired_BM.iloc[-1] - desired_BM.iloc[1])/2
    bool_growing = ((desired_BM.iloc[-1]-desired_BM.iloc[-5])/desired_BM.iloc[-1]).apply(lambda x: x > 1e-10)
    for k, bool_grow in bool_growing.items():
        if bool_grow:
            c_max_gr[k] = desired_BM[k].iloc[-6]
    biomass_diff = (desired_BM.iloc[-1]-desired_BM.iloc[0])
    start = desired_BM.iloc[0] + biomass_diff*0.1
    end = desired_BM.iloc[0] + biomass_diff*0.9
    return c_max_gr, start, end, bool_growing



def get_desired_cycle(Biomass_df, log_step=5):
    def correct_cycle(cycle): # 
        if cycle < log_step:
            return log_step
        return round(cycle / log_step) * log_step

    def get_growth_phase_length():
        return ((desired_cycle['end'] - desired_cycle['start'])*(1-desired_cycle.bool_growing) + # if not growing, not changing growth phase length
                1e4*(desired_cycle.bool_growing)) #if growing, set growth length to 999
    
    def split_index_to_cols(df):
        items = df.index.str.split('_')
        columns = ['Species', 'Gene_inhibition', 'culture']
        if len(items[0]) > 3:
            columns.extend(['alpha_lv_pairs'])
        return pd.DataFrame(items.tolist(), index=df.index, columns=columns)
    desired_biomass_df = pd.DataFrame(get_maximum_growth_cycle(Biomass_df), index=['c_max_gr', 'start', 'end', 'bool_growing'])
    
    desired_cycle = (desired_biomass_df.iloc[:-1]
                .apply(lambda x: 
                        search_gr_cycle_with_biomass(Biomass_df.loc[:,x.name],x))
                .T)
    desired_cycle['bool_growing'] = desired_biomass_df.T.bool_growing
    desired_cycle['cycle_max_gr'] = desired_cycle['c_max_gr'].apply(correct_cycle) # -> cycle_max_gr
    desired_cycle['growth_phase'] = desired_cycle[['start', 'end']].values.tolist()
    desired_cycle['growth_phase_length'] = get_growth_phase_length()
    desired_cycle['end_cycle'] = Biomass_df.index[-1]//5*5
    desired_cycle = desired_cycle.join(split_index_to_cols(desired_cycle))
#     .query('culture=="coculture"')
    if len(desired_cycle.Gene_inhibition.unique()) >1:
        desired_cycle = desired_cycle.set_index('Gene_inhibition')[['cycle_max_gr', 'bool_growing', 'growth_phase','growth_phase_length', 'Species','culture','end_cycle']]
    else:
        desired_cycle.index = ['_'.join([x[1],x[3]]) for x in desired_cycle.index.str.split('_')]
        desired_cycle.Gene_inhibition = desired_cycle.index

    return desired_cycle