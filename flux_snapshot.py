from setup import (load_models,
                   convert_arg_to_list,
                   get_gene_id,
                   remove_Zero_col)
from dataclasses import dataclass, field
from typing import List, Dict, Any
import pandas as pd
import cobra
import itertools
import re
import json
import ast
import numpy as np
# from growth_summary import (get_maximum_growth_cycle)

def get_rcts_list(model, gcomb_list): 
    rcts_list = list()
    rcts_set = set()
    for i, gene in enumerate(gcomb_list):
        gene_rcts = [rct.id for rct in model.genes.get_by_id(get_gene_id(model, gene)).reactions]
        if i > 0:
            gene_rcts = list(set(gene_rcts) - rcts_set)
        rcts_list.append(gene_rcts)
        rcts_set = rcts_set | set(gene_rcts)
    return rcts_list

def adjust_flux_df(model, df, gene_combo: list, alpha_table:pd.DataFrame): 
    # Don't use v1 cols to indicate gene_inhibition, will skip unidirectional reaction
    def query_alpha(gene_combo, alpha_table):
        splitted = gene_combo.split('.')
        if len(splitted) == 2:
            gcomb_alpha = {gene: alpha_table.loc[gene, f'{model.id}'] for gene in gene_combo.split('.')}
            # return gcomb_alpha

        # alphas
        splitted = gene_combo.split('_')
        # Gene_inhibition, lv_pair = [ele.split('.') for ele in 'folP.folA_1.2'.split('_')]
        Gene_inhibition, lv_pair = [ele.split('.') for ele in gene_combo.split('_')]
        lv_pair = tuple([int(ele) for ele in lv_pair])
        gcomb_alpha = dict()
        alpha_table = alpha_table.query('lv_pairs == @lv_pair')
        for current_gene in Gene_inhibition:
            # print(alpha_table.head())
            gcomb_alpha.update({current_gene: 
                                alpha_table.loc[current_gene, f'{model.id}']})
        return gcomb_alpha, alpha_table

    print(gene_combo)


    if 'Normal' not in gene_combo:
#         gene_combo_dict = get_gcomb_alpha_dict(gene_combo) 
        v1_cols = df.filter(regex='v1').columns
        
        # gcomb_alpha = {gene: alpha_table.loc[gene, f'{model.id}'] for gene in gene_combo.split('.')}
        gcomb_alpha, alpha_table = query_alpha(gene_combo, alpha_table)
        # print(gcomb_alpha)
        gcomb_alpha = dict(sorted(gcomb_alpha.items(), key=lambda item: item[1], reverse=True))
        
        rcts_list = get_rcts_list(model, gcomb_alpha.keys()) # exclude repeated rct for two gene
        scaled_rcts = list()
    # return df
        for gene, rcts in zip(gcomb_alpha.keys(), rcts_list):  
#             rcts = [rct for rct in rcts if rct in orig_cols]
            for orig_col in rcts:
                alpha = alpha_table.loc[f'{gene}', f'{model.id}']    
                v1_col = orig_col + "_v1"
                reversible = v1_col  in v1_cols
                if reversible: 
                    df[f'{orig_col}'] = (df[f'{orig_col}'] + df[f'{v1_col}'])/alpha # only forward or backward != 0  
                    df = df.drop(f'{v1_col}', axis=1)  
                else:
                    df[f'{orig_col}'] = (df[f'{orig_col}'])/alpha # only forward or backward != 0              
    return df

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


def get_growth_phase_summary(model, culture, df_from_js, desired_cycle, current_gene_combo, is_flux=True):
    
    growth_phase = desired_cycle.query('Species == @model.id & culture == @culture').loc[current_gene_combo, 'growth_phase']
    growth_df = df_from_js.query('cycle>@growth_phase[0] & cycle<@growth_phase[1]').drop('cycle', axis=1)
#     growth_df = growth_df.drop(['cycle','x','y'], axis=1) if is_flux else growth_df.drop(['cycle','x','y'])
    summary_df = (growth_df.apply(['mean','std'])
         .unstack().to_frame().T)
    summary_df.columns = summary_df.columns.map('_'.join).str.strip('_')
    summary_df.index= [current_gene_combo]
    return summary_df


def set_GI_SP_as_MI(end_BM, manual_SI = False, multi_index = None):
    if type(end_BM.index[0]) is not int: # Gene_inhibition is not index # Series.dtype -> dtype('int64')
        if end_BM.index.name in [None, 'DG']:
            end_BM.index.name = 'Gene_inhibition'
        if end_BM.index.name == 'Gene_inhibition':
            if 'Gene_inhibition' in end_BM.columns:
                end_BM = end_BM.drop('Gene_inhibition', axis=1)
            end_BM = end_BM.reset_index()
    if multi_index:
        return end_BM.set_index(multi_index)
    if 'Species' in end_BM.columns and manual_SI == False:
        return end_BM.set_index(['Gene_inhibition', 'Species'])
    else:
        return end_BM.set_index('Gene_inhibition')


def join_dfs_using_MI(df_list, how='left', multi_index = None):
    df_list = [set_GI_SP_as_MI(df, multi_index=multi_index) for df in df_list]
    result_df = df_list[0]
    for df in df_list[1:]:
        result_df = result_df.join(df, how=how)
    return result_df


all_culture_items = ['coculture_media', 'E0_coculture_flux', 'S0_coculture_flux', 'E0_monoculture_media', 'E0_monoculture_flux', 'S0_monoculture_media', 'S0_monoculture_flux']
def get_desired_keys(culture_options, item_options, all_culture_items=all_culture_items):
    culture_options, item_options = '|'.join(culture_options), '|'.join(item_options)
    search_str = re.compile(f'{culture_options}.*[{item_options}]')
    desired_keys = [ele for ele in all_culture_items if re.search(search_str, ele)]
    return desired_keys


def retrive_specific_culture(file_list, culture_options, item_options, XG_df_list):
    out_dict = dict()
    # # mono_dict = dict()
    # file_list = j
    # file_list, XG_df_list = convert_arg_to_list(file_list), [SG_cycle, DG_cycle]
    first_file = file_list[0]
    for file, XG_df in itertools.zip_longest(file_list, XG_df_list, fillvalue=first_file):
        if isinstance(file, str):
            file = json.load(open(file))
        XG_list = XG_df.index

        for outer_key in get_desired_keys(culture_options=culture_options, item_options=item_options):
            out_dict['_'.join([outer_key, XG_df.index.name])] = {current_gene: file[outer_key][current_gene] for current_gene in XG_list}                                            
    return out_dict

def get_alpha_to_merge(E0, S0, alpha_table): # TODO: check necessity- use in normal drug comb
    def get_alpha(l: list[str, str], alpha_table):
        l=l[0]
    #     print(l)
        Species = [E0.id, S0.id]
        df = pd.DataFrame(alpha_table.loc[l, Species]).T
        return df.values.tolist()
    # a.apply(get_alpha, axis=1, result_type='broadcast')

    S = [ele for ele in alpha_table.columns if 'S0' in ele]
    ES_cols = ['E0', S[0]]
    alpha_wide = get_alpha_wide(alpha_table)
    df = (get_SG_DG(explode=False)
          .apply(lambda x: pd.Series(get_alpha(x, alpha_table), index=ES_cols), axis=1) # df with E0 S0 alpha columns
          .melt(value_vars=ES_cols, var_name='Species', value_name='alpha', ignore_index=False)
          .combine_first(alpha_wide))
    return df

def modify_checkerboard_flux_compare_df(alpha_table, df):
    def get_ICX_from_lv_pair_list(alpha_table, lv_pair_list=[4.2, 5.3]):
        ICX = dict()
        corrected_lv_pair_list = list()
        print(lv_pair_list)
        for lv_pair in convert_arg_to_list(lv_pair_list):
            if isinstance(lv_pair, str):
                lv_pair = str(lv_pair).split('.')
                lv_pair = tuple([int(x) for x in lv_pair])
            corrected_lv_pair_list.append(lv_pair)
            # print(alpha_table.query('lv_pairs==@lv_pair').ICX.unique())
            ICX.update({lv_pair: alpha_table.query('lv_pairs==@lv_pair').ICX.to_dict()})
        return ICX 

    def get_lv_pairs_col(flux_compare_df):
        lv_pairs = [lv_pair[1].split('.') for lv_pair in flux_compare_df.index.str.split('_')]
        lv_pairs = [tuple([int(ele) for ele in lv_pair]) for lv_pair in lv_pairs]
        return lv_pairs
    df['lv_pairs'] = get_lv_pairs_col(df)
    cols = df.columns.tolist()
    df['ICX'] = df.apply(lambda x : get_ICX_from_lv_pair_list(alpha_table, [x.lv_pairs]).values(), axis=1)
    cols = ['Species', 'XG', 'lv_pairs', 'ICX'] + [col for col in cols if col not in ['Species', 'XG', 'lv_pairs']]
    df = df[cols]
    return df

@dataclass(kw_only=True)
class FluxCompare:
    desired_cycle : pd.DataFrame
    file_list : None
    E0 : cobra.Model = None
    S0 : cobra.Model = None
    # species_options : List = field(default_factory=lambda: [E0, S0])
    culture_options : str = field(default_factory=lambda: ['coculture'])
    item_options : str = field(default_factory=lambda: ['media', 'flux'])
    flux_dict : Dict = field(default_factory=dict)
    flux_compare_dict : Dict = field(default_factory=dict)
    flux_compare_df : pd.DataFrame = field(default_factory=pd.DataFrame)
    desired_cycle_col : str = 'cycle_max_gr'
    is_checkerboard :bool = False
    alpha_table :pd.DataFrame = None
    count = 0

    def __post_init__(self):
        self.SG_cycle, self.DG_cycle = get_XG_cycle_from(self.desired_cycle)
        self.SG_list, self.DG_list = self.SG_cycle.index.to_list(), self.DG_cycle.index.to_list()
        self.XG_df_list = [self.SG_cycle, self.DG_cycle]
        self.culture_options, self.item_options = convert_arg_to_list(self.culture_options), convert_arg_to_list(self.item_options)
        self.file_list = convert_arg_to_list(self.file_list)
        if 'lv_pairs' in self.alpha_table.columns:
            self.is_checkerboard = True

    def retrive_specific_culture(self):
        out_dict = dict()
        first_file = self.file_list[0]
        for file, XG_df in itertools.zip_longest(self.file_list, self.XG_df_list, fillvalue=None):
            if isinstance(file, str):
                print('load from json file', str)
                file = json.load(open(file))
                previous_js = file
            elif file is None:
                file = previous_js

            XG_list = XG_df.index

            for outer_key in get_desired_keys(culture_options=self.culture_options, item_options=self.item_options):
                out_dict['_'.join([outer_key, XG_df.index.name])] = {current_gene: file[outer_key][current_gene] for current_gene in XG_list}                                            
        self.flux_dict = out_dict
        return None
    
    def get_flux_dict_per_culture_item(self, model, sub_flux_dict, culture_item, desired_cycle, log_step=globals().get('log_step', 5)): # desired_cycle as SG_cycle, DG_cycle
        def get_flux_snapshot(s : pd.Series, no_grow, temp_flux_dict, culture):
            cycle=s[self.desired_cycle_col]
            current_gene_combo = s.name
            print(model.id, culture_item, self.count)
            if cycle < 10:
                no_grow.append(current_gene_combo)
                temp_flux_dict[current_gene_combo] = pd.DataFrame(index = [f'{current_gene_combo}'])
            else:
                desired_js = sub_flux_dict.get(f'{current_gene_combo}')
                # self.s = desired_js
                # coculture_dict['E0_monoculture_flux_SG']['folP.folA_1.0'] 

                df_from_js = pd.DataFrame.from_dict(desired_js)
                is_flux = '_flux' in culture_item

                if is_flux:
                    df_from_js = adjust_flux_df(model, df_from_js, current_gene_combo,alpha_table=self.alpha_table).drop(['x','y'], axis=1)
                        
                else:
                    df_from_js = df_from_js.pivot(index=['cycle'],columns=['metabolite'], values='conc_mmol').reset_index()
                growth_phase_summary = get_growth_phase_summary(model, culture, df_from_js, desired_cycle, current_gene_combo, is_flux=is_flux)
                
                invalid_snapshot, max_correction = True, 3
                while invalid_snapshot:
                    print('cyc', cycle, 'mxc', max_correction, max_correction <= -1)
                    temp_df = df_from_js.query("cycle == @cycle")
                    temp_df.index = [current_gene_combo]
                    temp_df = pd.concat([temp_df, growth_phase_summary], axis=1) # mean, sd for each reaction -> ncols*3
                    print('snap&m',float(abs(temp_df['EX_bulk_ac_e'])),float(abs(temp_df['EX_bulk_ac_e_mean'])))
                    if np.isnan(temp_df['EX_bulk_ac_e_mean'][0]) or float(abs(temp_df['EX_bulk_ac_e'])) > float(abs(temp_df['EX_bulk_ac_e_mean']))*.7:
                        invalid_snapshot = False
                    else:
                        max_correction -= 1
                        cycle -= log_step                        
                        if max_correction <= -1 or cycle < 10:
                            invalid_snapshot = True
                            temp_df.loc[current_gene_combo, 'EX_bulk_ac_e'] = temp_df.loc[current_gene_combo, 'EX_bulk_ac_e_mean']
                        
                temp_flux_dict.update({current_gene_combo: temp_df}) 
                
                print(temp_flux_dict.keys())
                return temp_flux_dict, no_grow
        temp_flux_dict, no_grow = dict(), list()
        culture = culture_item.split('_')[0]
        for current_gene_combo, cycle_row in (desired_cycle
                                              .query('Species == @model.id & culture == @culture')
                                            #   .iloc[:1]
                                              .iterrows()):
            print(current_gene_combo,'currr')
            temp_flux_dict, no_grow = get_flux_snapshot(cycle_row, no_grow, temp_flux_dict, culture)
            
        # print(temp_flux_dict)
        if no_grow:
            print('Zero growth: ', ', '.join(no_grow))
        if not temp_flux_dict:
            return pd.DataFrame()

        return pd.concat(temp_flux_dict.values()).copy()

    def get_flux_compare_df(self): # generation of mono_flux, coculture_flux, exclude media
        temp_flux_dict = dict()
        if not self.flux_dict:
            self.flux_dict = self.retrive_specific_culture()

        for key, sub_flux_dict in self.flux_dict.items(): 
            print(key)
            if 'media' not in key:
                current_species, current_XG, culture = key.split('_')[0], key.split('_')[-1], key.split('_')[1]
                XG_cycle = self.SG_cycle if current_XG == 'SG' else self.DG_cycle
                model = self.E0 if current_species == 'E0' else self.S0
                current_species = model.id
                culture_item = culture + '_flux' 

                print(current_species, current_XG)
                flux_df = self.get_flux_dict_per_culture_item(model, sub_flux_dict, culture_item, XG_cycle) # flux snapshot for all GI as index 
                
                flux_df['Species'] = current_species
                flux_df['XG'] = current_XG
                flux_df['culture'] = culture
                temp_flux_dict[key] = flux_df # flux_df all gene_inhibition as index
                print(key, '--------', flux_df.index)
                
        flux_compare_df = pd.concat(temp_flux_dict.values()).copy()
        flux_compare_df.index.name = 'Gene_inhibition'
        flux_compare_df = remove_Zero_col(flux_compare_df)
    #     flux_compare_df = get_ESdiff(remove_Zero_col(flux_compare_df)) # ES diff 41 columns ?only common
    #     flux_compare_df = get_Ndiff(flux_compare_df)
        if not self.is_checkerboard: #  alpha value not joined for checkerboard, SG format different 
            flux_compare_df = join_dfs_using_MI(get_alpha_to_merge(self.alpha_table), flux_compare_df, how='right') # DG only
        else:
            flux_compare_df = modify_checkerboard_flux_compare_df(self.alpha_table, flux_compare_df)
        self.flux_compare_df = flux_compare_df
        return flux_compare_df