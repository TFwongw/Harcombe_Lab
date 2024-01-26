import numpy as np
import pandas as pd
from setup import convert_arg_to_list
from flux_snapshot import set_GI_SP_as_MI

def get_Biomass_df(files):
    return pd.concat(
            [pd.read_csv(file, index_col='cycle')
             for file in convert_arg_to_list(files)]
        ,axis=1)

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

# BM to merge with p_o to form full flux dict
def get_BM_bin(BM): # separate S0 E0 bin, combine again
    sub_BM = BM.pivot(columns='Species', values='Biomass') #merge gene and species
    sub_BM['S0_BM_bin'] =  pd.cut(sub_BM['S0'], [-1,.001,.0025,1], labels = ['Low', 'Medium','High'])
    sub_BM['E0_BM_bin'] =  pd.cut(sub_BM['E0'], [-1,.001,.0065,1], labels = ['Low', 'Medium','High'])
    return sub_BM.melt(value_vars=['S0_BM_bin','S0_BM_bin'],value_name='BM_bin',ignore_index=False).BM_bin

def get_end_BM(Biomass_df):
    BM = pd.DataFrame(Biomass_df.loc[:, [col for col in Biomass_df.columns if 'coculture' in col]].iloc[-1]            )
    BM.columns = ['Biomass']
    BM['Species'] = list(pd.Series(BM.index).apply(lambda x: str(x).split('_')[0]))
    BM['Gene_inhibition'] = list(pd.Series(BM.index).apply(lambda x: str(x).split('_')[1]))
    BM = BM.set_index(['Gene_inhibition'])
    BM['Total_BM'] = BM.groupby('Gene_inhibition').Biomass.sum()

    BM['Total_BM_bin'] = pd.cut(BM['Total_BM'], [-1,.004,.008,1], labels = ['Low', 'Medium','High'])
    BM['BM_bin'] = get_BM_bin(BM)
    
    end_BM = set_GI_SP_as_MI(BM).join(set_GI_SP_as_MI(add_to_end_BM(BM)))
    # lcts_used = .1-metab_end.lcts_e.fillna(0)
    # end_BM = end_BM.join(lcts_used).rename(columns={'lcts_e': 'lcts_consumed'})
    # end_BM['Standzrdized_Carbon_source_conversion_efficiency'] = end_BM.apply(lambda x: x.Total_BM/x.lcts_consumed,axis=1)
    # end_BM['Standzrdized_Carbon_source_conversion_efficiency'] = end_BM.Standzrdized_Carbon_source_conversion_efficiency/end_BM.loc['Normal', 'Standzrdized_Carbon_source_conversion_efficiency'][0] # series of 2 element same value for E0, S0
    return end_BM

def separate_Species_df(df, model_id, inc_Species = False):
    if not isinstance(model_id, str):
        model_id = model_id.id
    def get_species_loc(model_id):
        return list(temp_df.Species == model_id)
    df = df.loc[:,~df.columns.str.contains('Ndiff|ESdiff')]
    temp_df = pd.concat([df['Species'],df.select_dtypes(include = ['float'])], axis=1)
    result_df = temp_df.loc[get_species_loc(model_id)]
    return result_df if inc_Species else result_df.drop('Species',axis=1)

def add_to_end_BM(end_BM):
    def get_ratioNstd_col(model_id):
        temp_df = separate_Species_df(end_BM, model_id, inc_Species=True)
#         temp_df[f'{model.id}_BM_ratio'] = temp_df.Biomass/temp_df.Total_BM
#         temp_df[f'{model.id}_standardized_BM'] = temp_df.Biomass/temp_df.loc['Normal', 'Biomass']
        temp_df[f'BM_consortia_frac'] = temp_df.Biomass/temp_df.Total_BM
        temp_df[f'standardized_BM'] = temp_df.Biomass/temp_df.loc['Normal', 'Biomass']
        return temp_df
    return pd.concat([get_ratioNstd_col('E0'), get_ratioNstd_col('S0')]).set_index('Species', append=True).drop(['Total_BM', 'Biomass'], axis=1)

