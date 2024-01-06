#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import cobra
from cobra.io import load_model
import cometspy as c
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools
from functools import partial
from collections.abc import Iterable
import re
from dataclasses import dataclass, field
from typing import Callable, List, Any, Dict
import sys
import itertools
import json
import logging
import multiprocessing
from functools import partial
import copy
import re

initial_pop = 1e-8
log_step = 5
current_gene = 'thrC'
precision = 2
n_processor = 10

if 'E0' in globals():
    # Code to execute if variable exists
    print("Variable exists, keep E0, S0")
else:
    # Code to execute if variable 'a' does not exist
    print('initialize E0, S0')
    E0, S0, all_components, media, alpha_table = None, None, None, None, None

# p = c.params()
# p.set_param("defaultKm", 0.00001) # M 
# p.set_param("defaultVmax", 10) #mmol/gDw/hr
# p.set_param("maxCycles", 100)
# # p.set_param("maxCycles", 20)
# # p.set_param("maxCycles", 160)
# p.set_param("timeStep", 1) 
# p.set_param('writeFluxLog', True)
# p.set_param('writeMediaLog', True)
# # function for iterate 3 species
p=None

# !module load java/openjdk-13.0.2
# !module load gurobi/9.0.2
# !unset PYTHONPATH
# !unset PYTHONHOME
# !unset PYTHONSTARTUP
# os.environ["GRB_LICENSE_FILE"]="/common/software/install/migrated/gurobi/license/gurobi.lic"

def load_models(Egal_reuptake=True, scale_force_gal_c=0,
                ac_scale=1, gal_scale=None, ACtex_ub=0, co2_scale=None, LCTStex_lb=0, cpath_scale = None, w=None):
    def get_all_components():
        all_metabolites = {
            model.id: model.metabolites for model in [E0, S0]
        }
        all_genes = {
            model.id: 
            [gene.name for gene in model.genes if type(gene) != str]
            for model in [E0, S0]
        }

        all_reactions = {
            model.id: model.reactions for model in [E0, S0]
        }

        all_components = {'metabolites': all_metabolites,
                        'genes': all_genes, 
                        'reactions': all_reactions}
        return all_components

    def rename_in_bulk(model):
        model.metabolites.ac_e.id = 'bulk_ac_e'
        model.metabolites.bulk_ac_e.name = 'bulk Acetate'
        model.metabolites.bulk_ac_e.formula = 'C20H30O20'
        model.reactions.EX_ac_e.id = 'EX_bulk_ac_e'
        model.reactions.EX_bulk_ac_e
    
    E0 = cobra.io.read_sbml_model("./models/iML1515_E0.xml")
    S0 = cobra.io.read_sbml_model("./models/STM_v1_0_S0.xml")
    E0.id, S0.id = 'E0', 'S0'
    
    # galactose reuptake in Diauxic environment result in disturbance in growth rate estimate after lcts being used up
    if Egal_reuptake == False:
        E0.reactions.GALtex.upper_bound = 0
        E0.reactions.ACtex.upper_bound = 0

    # alpha_table = pd.read_csv('./Data/alpha_table.csv',index_col=0)

    # # alpha_table = pd.read_csv('./Data/alpha_table_wgal.csv',index_col=0)
    # alpha_table.columns = ['E0', 'S0', 'S0']
    # alpha_table = alpha_table.iloc[:, [0,2]] # only E0, S0 columns

    nutrient_medium = {'EX_ca2_e': 10,
    'EX_cl_e': 10,
    'EX_cobalt2_e': 10,
    'EX_cu2_e': 10,
    'EX_fe2_e': 10,
    'EX_fe3_e': 10,
    'EX_k_e': 10,
    'EX_mg2_e': 10,
    'EX_mn2_e': 10,
    'EX_mobd_e': 10,
    'EX_ni2_e': 10,
    'EX_o2_e': 10,
    'EX_pi_e': 10,
    'EX_so4_e': 10,
    'EX_zn2_e': 10,
    'EX_nh4_e': 10}
    E0.medium = nutrient_medium 
    S0.medium = nutrient_medium 

    c_limiting_conc = 10 # corresponding to maximum uptake rate in carbon limited batch culture
    met_limiting_conc = 0.08 # DO not search momoculture alpha with met_limitng_conc

    change_medium(E0, ['EX_lcts_e', 'EX_met__L_e'], [c_limiting_conc , 10], True)
    change_medium(S0, 'EX_gal_e', c_limiting_conc, True)
    
    
    # add_10AC_export(E0, S0, all_components)
    E0.reactions.ACt2rpp.add_metabolites({E0.metabolites.ac_p: .9})
    S0.reactions.ACt2rpp.add_metabolites({S0.metabolites.ac_p: .9})
    
    all_components = get_all_components()
    if w:
        weight_carbon_byproduct(E0, S0,all_components, ac_scale=ac_scale, gal_scale=gal_scale, co2_scale=co2_scale, ACtex_ub=ACtex_ub, LCTStex_lb=LCTStex_lb, cpath_scale =cpath_scale) # reduce flux weight of carbon byproduct from parsimonious fba, weight of galactose fluxes as 2carbon  
        rename_in_bulk(E0), rename_in_bulk(S0)
        
    return E0, S0, all_components

def change_medium(model, metabolite, value, return_medium=False): # function for setting medium metabolite value
    medium = model.medium
    if not isinstance(value, Iterable):
        metabolite, value = [metabolite], [value]
#         print(value)
    for m, v in zip(metabolite, value):
        medium[m] = v
    model.medium = medium
    if return_medium:
        return model.medium

def substrate_only(rxn, met):
    return rxn.lower_bound>=0 and met in rxn.reactants        

def substrate_any(rxn, met):
    return (rxn.lower_bound>=0 and met in rxn.reactants) or rxn.lower_bound<0

def product_only(rxn, met):
    if re.search(r"_e$", rxn.id):
        return True
    return (rxn.lower_bound>=0 and met in rxn.products) or (rxn.upper_bound<=0 and met in rxn.reactants)   
    
def get_component(model, query, all_components):
    ID = model.id.split('.')[0]
    if query in all_components['genes'][ID]:
        return model.genes.get_by_id(get_gene_id(model, query))
    if query in all_components['metabolites'][ID]:
        return model.metabolites.get_by_id(query)
    return model.reactions.get_by_id(query)
    
def get_links_component(model, query, all_components, id_only=True, is_sub_only=None, is_any_sub=None, is_any_prod=None, is_prod_only=None):
    def get_result(model, query, all_components):        
        query = get_component(model, query, all_components)
        if type(query) in [cobra.Reaction]:
            result = set(query.metabolites)
        else: 
            result = query.reactions 
            if type(query) is cobra.Metabolite:
                if is_sub_only is not None:
                    result = set([rxn for rxn in result if substrate_only(rxn, query)==is_sub_only])
                elif is_any_sub is not None or is_any_prod is not None:
                    result = set([rxn for rxn in result if substrate_any(rxn, query)==is_any_sub])
                elif is_prod_only is not None:
                    result = set([rxn for rxn in result if product_only(rxn, query)==is_prod_only])
        return result
    
    result=set()
    for ele in convert_arg_to_list(query):
        result.update(get_result(model, ele, all_components))
    if id_only:
            result = set([ele.id for ele in result])
    return list(result)

def retrive_model_in_medium(key, media):
    model = E0 if 'E0' in key else S0
    model.medium = media[''.join([ele for ele in media.keys() if key in ele])]
    model.id = key.replace('_limited', '')
    return model 

def get_media():
    media = dict()
    c_limiting_conc = 1
    met_limiting_conc = 0.08

    # more met_limiting
    with E0, S0:
    #     media['E0_unlimited'] = lactose_met_medium
        with E0:
            media['E0.lcts_limited'] = change_medium(E0, ['EX_lcts_e', 'EX_met__L_e'], [c_limiting_conc, 10], True)
        media['E0.Met_limited'] = change_medium(E0, ['EX_lcts_e', 'EX_met__L_e'], [10, met_limiting_conc], True)
        with S0:
            media['S0.ac_limited'] = change_medium(S0, 'EX_ac_e', c_limiting_conc*6, True)
        media['S0.gal_limited'] = change_medium(S0, 'EX_gal_e', c_limiting_conc*2, True)
    model_ids = [ele.split('_')[0] for ele in media.keys()]
    return media

def get_list_target_obj_val():
    list_target_obj_val = {media_key: retrive_model_in_medium(media_key, media).slim_optimize()/2
                        for media_key in media.keys()}
    E0.id, S0.id = 'E0', 'S0'
    return list_target_obj_val
# list_target_obj_val = get_list_target_obj_val()

def search_gr_cycle_with_biomass(df_search, biomass_values):
    return [df_search[df_search >= biomass_value].idxmin() 
                for biomass_value in list(biomass_values)]

def get_maximum_growth_cycle(desired_BM, scale_diff=0.1):
    c_max_gr = desired_BM.iloc[1]+ (desired_BM.iloc[-1] - desired_BM.iloc[1])/2
    bool_growing = ((desired_BM.iloc[-1]-desired_BM.iloc[-5])/desired_BM.iloc[-1]).apply(lambda x: bool(x > 1e-10))
    for k, bool_grow in bool_growing.items():
        if bool_grow:
            c_max_gr[k] = desired_BM[k].iloc[-6]
    biomass_diff = (desired_BM.iloc[-1]-desired_BM.iloc[0])
    start = desired_BM.iloc[0] + biomass_diff*scale_diff
    end = desired_BM.iloc[0] + biomass_diff*(1-scale_diff)
    return c_max_gr, start, end, bool_growing

def get_desired_cycle(Biomass_df, log_step=5, scale_diff=0.1):
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
    desired_biomass_df = pd.DataFrame(get_maximum_growth_cycle(Biomass_df, scale_diff=scale_diff), index=['c_max_gr', 'start', 'end', 'bool_growing'])
    
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
    
    if len(desired_cycle.Gene_inhibition.unique()) >1 or len(desired_cycle)<2:
        desired_cycle = desired_cycle.set_index('Gene_inhibition')
    else:
        desired_cycle.index = ['_'.join([x[1],x[3]]) for x in desired_cycle.index.str.split('_')]
        desired_cycle.Gene_inhibition = desired_cycle.index
    # print(desired_cycle.columns)
    return desired_cycle   

def create_common_media(Species, carbon_source = "lcts_e", carbon_source_val = 5e-3, 
                        nutrients_val = 100, add_nutrient = '', add_nutrient_val = [100]):    # carbon_source = 'lcts_e', 'ac_e' or 'glc__D_e'
    # add_nutrient = 'met__L_e' for E.coli monoculture
    
    # convert into list for enumeration 
    if type(add_nutrient)  == str:
        add_nutrient = [add_nutrient]
    if type(add_nutrient_val)  == int:
        add_nutrient_val = [add_nutrient_val]
        
    l = c.layout(Species)
    
    base_nutrients = ["ca2_e", "cl_e", "cobalt2_e", "cu2_e","fe2_e", "fe3_e", "k_e","mg2_e",
                  "mn2_e", "mobd_e", "ni2_e", "o2_e", "pi_e", "so4_e", "zn2_e","nh4_e"]
    for nutrient in base_nutrients:
        l.set_specific_metabolite(nutrient, nutrients_val)
        
    if (add_nutrient != ['']):
        if (len(add_nutrient) == len(add_nutrient_val)):
            for _,i in enumerate(zip(add_nutrient, add_nutrient_val)): 
                l.set_specific_metabolite(i[0], i[1])
        else:
            print(f'Set all additional nutrients to {add_nutrient_val[0]}')
            for _,i in enumerate(add_nutrient): 
                l.set_specific_metabolite(i, add_nutrient_val[0])
    
    l.set_specific_metabolite(carbon_source, carbon_source_val)  
    return(l)
   
def get_gene_id(model, gene_name):
    for i in model.genes:
        if(i.name == gene_name):
            return(i.id)

# def sim_cultures(model,suffix='', p=p): # if len(model) = 1 -> mono, else coculture
    
def sim_cultures(model, layout, p=p, base = None, genes=''): # non_functional argument model
    # separate function into mono & coculture to prevent using wrong layer
    if type(layout) is list:
        layout = layout[0] # layout object stored inside list of one element unpacking from iter_species
    
    sim = c.comets(layout, p) 
    sim.working_dir = base
    print(sim.working_dir)
    
    try:
        sim.run()
        # if modify_bool==True:
        #     sim.total_biomass.rename(columns={'S0.ac':'S0.glc'},inplace=True)
    except:
        logging.exception(f"{sim.run_output}")
        print(f"{sim.run_output}")
    biomass_df = sim.total_biomass.set_index('cycle')
    biomass_df.columns =  rename_columns(biomass_df)
    return biomass_df, sim

def convert_arg_to_list(arg):
    if type(arg) in [pd.Series, pd.Index]:
        arg = list(arg) 
    elif type(arg) not in [list, tuple, set]:
        arg = [arg] # differ from list(arg) -> conversion of str
    return arg

def create_c(model, initial_pop=initial_pop, obj_style='MAX_OBJECTIVE_MIN_TOTAL'):
    model = c.model(model)
    model.open_exchanges()
    model.initial_pop = [0, 0, initial_pop]
    model.obj_style = obj_style
    return(model)

def iter_species(models,f,*args,**kwargs): 
    r_object = list()
    if type(models) is not zip:
        for model in models:
            r_object.append(f(model,*args,**kwargs))
    else:
        for model, *extra_objects in models: 
            r_object.append(f(model,extra_objects,*args,**kwargs))
    return(r_object) 

def get_alpha_steps(SG, steps=np.arange(0.2,2,0.2)):
    result_df = pd.DataFrame()
    for scaling in steps:
        temp_df = alpha_table.loc[[SG]]*scaling
        temp_df.index = temp_df.index + f'_{round(scaling,3)}'
        result_df = pd.concat([result_df, temp_df])
    return result_df

def get_alphas_from_tab(model, genes: list, alpha_table): # the returned df of one gcomb index only pass as series
    if isinstance(alpha_table, pd.Series):
        alpha_table = alpha_table.to_frame().T
    genes = convert_arg_to_list(genes)
    # print(model, genes, alpha_table.head())
    print([gene for gene in genes])
    alphas = [alpha_table.loc[gene ,model.id] for gene in genes]
    return alphas 

def scale_reaction(model, reaction_id, alpha, direction='forward'):
    # print(reaction_id, model.reactions.get_by_id(reaction_id).metabolites.values())
    if direction == 'forward':
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v >= 0)    # obetain 'end product'
    else:  # reverse reaction
        model.reactions.get_by_id(reaction_id).upper_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v <= 0)
    if isinstance(alpha, pd.Series):
        alpha = float(alpha)

    for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
        model.reactions.get_by_id(reaction_id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite
    # print(alpha, str(model.reactions.get_by_id(reaction_id)))
    return None

def separate_reaction(model, reaction_id, alpha):
    """
    decompose reversible reaction into forward and backward
    scale each end product by 1/alpha
    """
    (lb, ub) = model.reactions.get_by_id(reaction_id).bounds 
    rct_ids = [reaction_id] #? 
    if(lb < 0 and ub !=0): # only perform if reaction is bidirectional
        rev_reaction = model.reactions.get_by_id(reaction_id).copy() # copy of target reaction 
        rev_reaction.id = f'{reaction_id}_v1' # redefine id for copy of reaction
        model.add_reactions([rev_reaction]) # add to model
    
        scale_reaction(model, rev_reaction.id, alpha, direction='backward')      
        
        rct_ids.append(rev_reaction.id) # list of id for reaction and reaction_v1 
    scale_reaction(model, reaction_id, alpha, direction='forward')
    return(rct_ids)

def alter_Sij(model, alphas = 1, genes = 'folA', ko=False):  
    # get objective value for corresponding alpha
    
    # print(genes)
    # if 'folP' in genes: 
    #     alphas = 5
    
    alphas= convert_arg_to_list(alphas[0]) if type(alphas) is list and len(alphas)==1 else convert_arg_to_list(alphas)  # unlist one layer from zip comprehension 
    genes =  convert_arg_to_list(genes)
    # print(model, genes,alphas)
    
    genes_dict = {gene: alpha for gene, alpha in zip(genes, alphas)}
    genes_sorted = sorted(genes_dict.items(), key=lambda x:x[1], reverse=True) #sort by magnitude of alpha
    rct_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    for current_gene, alpha in genes_sorted:
        current_gene = current_gene.split('_')[0] # for step_alpha
        for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions: 
            if ko:
                # print('-----------ko-----------')
                model.reactions.get_by_id(rct.id).knock_out()
            elif (rct.id not in rct_ids):
                rct_ids.extend(separate_reaction(model, rct.id, alpha))# copy of reaction, forward_terminal_change = True
    return(rct_ids) 

def create_layout_object(E_model, S_model, carbon_source_val=5e-3, add_nutrient_val=[100], co_met=[0], co=True, mono=True, mono_S=True, Smono_carbon_source='bulk_ac_e'):
    add_nutrient_val = convert_arg_to_list(add_nutrient_val)
    partial_create_common_media = partial(create_common_media, carbon_source_val=carbon_source_val)
    
    # co_layout = partial_create_common_media([E_model, S_model], carbon_source='lcts_e') if co else None
    co_layout = partial_create_common_media([E_model, S_model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=co_met) if co else None
    E0_layout = partial_create_common_media([E_model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=add_nutrient_val) if mono else None
    
    Scarbon_source_val = carbon_source_val/10 if Smono_carbon_source=='bulk_ac_e' else carbon_source_val # ac_e as bulk 
    print('set Scarbon source as', Smono_carbon_source, Scarbon_source_val)
    
    S0_layout = create_common_media([S_model], carbon_source=Smono_carbon_source, carbon_source_val=Scarbon_source_val) if mono_S else None
    # return [co_layout, E0_layout, S0_ac_layout, S0_glc_layout]
    return [co_layout, E0_layout, S0_layout]

def unpack_output_object(output):
    df, *sim_object = zip(*output)
    return df, sim_object
# df, sim_object = unpack_output_object(mono_output)

def flux_correction(correction_factor, psol):
    print('flux correction', correction_factor)
        # correction_factor = {'ac': ac_scale, 'gal': gal_scale, 'co2': co2_scale}
    if correction_factor is None:
        return psol
    for key, scale in correction_factor.items():
        multiply_factor = scale if key == 'ac' else 1/scale
        
        print(f'EX_{key}_e' in psol, 'is in ')
        if f'EX_{key}_e' in psol:
            pass
            # print(key, 'not scaled',psol[f'EX_{key}_e'])
            # psol[f'EX_{key}_e'] = psol.loc[f'EX_{key}_e']*multiply_factor
        else: 
            print(f'EX_{key}_e', 'not in flux')
    return psol
        
def extract_dfs_from_sim_object(sim_object, flux_correction_factor):
    species_name = [ele for ele in sim_object.total_biomass.columns[1:]]
    if len(species_name)>1:
        culture = 'coculture' 
        out_dict = {f'{culture}_media' : sim_object.media}
    else: 
        culture = f'monoculture'
        out_dict = {f'{species_name[0]}_{culture}_media' : sim_object.media}
    for species in species_name:
        flux_df =  sim_object.fluxes_by_species[f'{species}']
        if species == 'E0':
            print('correct E')
            flux_df = flux_correction(flux_correction_factor, flux_df)
            print('maxac', flux_df['EX_bulk_ac_e'].abs().max())
        out_dict[f'{species}_{culture}_flux'] = flux_df
    return out_dict
        
# unpack_output_object
def extract_dfs(mono_sims, co_sim, flux_correction_factor=None):
    out_dict = extract_dfs_from_sim_object(co_sim, flux_correction_factor) if co_sim else dict() # initialize out_dict from coculture if have co_sim object 
    if mono_sims:
        for sim_object in mono_sims[0]:
            out_dict.update(extract_dfs_from_sim_object(sim_object, flux_correction_factor))
    out_dict = {k: v.to_dict() for k,v in out_dict.items()}
    return out_dict

# def force_gal_c(E0, scale_force_gal_c=0):
#     # tied GALt2pp to GALtex, force transport into cytosol
#     def get_metabolites_to_add(E0):
#         GALt2pp_metab = E0.reactions.GALt2pp.metabolites
#         reverse_GALt2pp_metab = {k:v*-1*scale_force_gal_c for k,v in GALt2pp_metab.items()}
#         return reverse_GALt2pp_metab
        
#     reverse_GALt2pp_metab = get_metabolites_to_add(E0)
#     E0.reactions.GALtex.add_metabolites(reverse_GALt2pp_metab)
#     # if scale_force_gal_c != 0:
#     if scale_force_gal_c != -1:
#         print(str(E0.reactions.GALtex))
#     return None
    
def additional_scales(E0, metab_scale=20, rxns=['LCTStex', 'LCTStpp', 'LACZ', 'HEX1', 
                                                'PGI', 'TALA', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK6', 'LDH_D', 'PYK3']):
    for rxn in rxns:
        if isinstance(rxn, str):
            rxn = E0.reactions.get_by_id(rxn)

        metab_to_scale = rxn.metabolites
        rxn.add_metabolites({k:v*metab_scale for k,v in metab_to_scale.items()})

def weight_carbon_byproduct(E0, S0, all_components, ac_scale=1, gal_scale=None, s_ac = True, co2_scale = None, ACtex_ub = 0, LCTStex_lb = 0,
    query_gal = ['gal_p','gal_c'], query_ac = ['ac_p', 'ac_e', 'ac_c'], query_co2 = ['co2_c', 'co2_p', 'co2_e'], 
    query_lac__D = ['lac__D_e', 'lac__D_c', 'lac__D_p'], query_for = ['for_e', 'for_c', 'for_p'], cpath_scale = None):
    
    print('LCTStex_lb', LCTStex_lb)
    
    # weight all exchange nigilible, req correction of flux
    # for rxn in [rxn for rxn in E0.reactions if 'tex' in rxn.id]:
    #     metab_to_scale = rxn.metabolites
    #     rxn.add_metabolites({k:v*5 for k,v in metab_to_scale.items()}) 
    
    
    print(ac_scale, gal_scale, co2_scale, ACtex_ub)
    print(gal_scale is not None)
    if gal_scale is not None:
        if gal_scale > 1:
            gal_scale = -1*(1-1/gal_scale) 
        for rxn in get_links_component(E0, query_gal, all_components, id_only=False, is_prod_only=True):
            metab_to_scale = rxn.metabolites
            rxn.add_metabolites({k:v*gal_scale for k,v in metab_to_scale.items()})
    # E0.reactions.GALtex.knock_out()
    # E0.reactions.GALtex.bounds = (1,1)
    
    if ac_scale is not None: # flux do not count in minflux
        add_scale = ac_scale-1 if ac_scale>=1 else ac_scale
        
        # for query_ac in [query_ac, query_lac__D, query_for]:
        for query_ac in [query_ac]:
            # print(get_links_component(E0, query_ac, all_components, is_prod_only=True))
            # for rxn in get_links_component(E0, query_ac, all_components, id_only=False, is_prod_only=True):
            # rxn_to_scale = get_links_component(E0, query_ac, all_components, id_only=False, is_prod_only=True)
            # rxn_to_scale.extend(['ACt2rpp'])
            # print(rxn_to_scale)
            # for rxn in rxn_to_scale:
            for rxn in [ 'ACt2rpp']:
                if isinstance(rxn, str):
                    rxn = get_component(E0, rxn, all_components)
                if not 'EX_' in rxn.id:
                    metab_to_scale = rxn.metabolites
                    rxn.add_metabolites({k:v*add_scale for k,v in metab_to_scale.items()}) 
                
        if ACtex_ub > 0:
            ACtex_ub = -1*ACtex_ub
        E0.reactions.ACtex.upper_bound = ACtex_ub/ac_scale if ac_scale>=0 else ACtex_ub
        # E0.reactions.ACtex.upper_bound = ACtex_ub
        # E0.reactions.ACtex.lower_bound = -10/ac_scale if ac_scale>=0 else -10
        
        # scale  cbyproduct bounds
        # E0.reactions.D_LACtex.upper_bound = ACtex_ub/ac_scale if ac_scale>=0 else ACtex_ub
        
        print(str(E0.reactions.ACtex), E0.reactions.ACtex.bounds)
        
    if co2_scale is not None:
        if co2_scale > 1:
            co2_scale = -1*(1-1/co2_scale) 
        for rxn in get_links_component(E0, query_co2, all_components, id_only=False, is_prod_only=True):
            metab_to_scale = rxn.metabolites
            rxn.add_metabolites({k:v*co2_scale for k,v in metab_to_scale.items()}) 
         
         
    cpath_rxns =  ['LACZ', 'HEX1', 'PGI', 'TALA', 'GAPD', 'PGK',
       'PGCD', 'PSERT', 'ACOTA', 'ACOTA', 'ACODA', 'ACt2rpp', 'ACtex', 'EX_ac_e']
    cpath_rxns =  ['ACODA', 'ACt2rpp']
    # cpath_rxns = ['LACZ', 'HEX1', 'PGI', 'TALA', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK6', 'LDH_D', 'PYK3']
    if cpath_scale is not None:
        cpath_scale = cpath_scale-1
        additional_scales(E0, metab_scale=cpath_scale, rxns=cpath_rxns)
        # additional_scales(E0, metab_scale=3, rxns=['ACtex', 'EX_ac_e']) # flux go through but no appraes
        # additional_scales(S0, metab_scale=3, rxns=['ACtex', 'EX_ac_e'])
        # E0.reactions.ACtex.upper_bound = ACtex_ub/cpath_scale if cpath_scale>=0 else ACtex_ub
   
    E0.reactions.ACtex.upper_bound = ACtex_ub
   
    # only actex by 10
    print(str(E0.reactions.ACtex), E0.reactions.ACtex.bounds)
    E0.reactions.LCTStex.lower_bound = LCTStex_lb
    print(E0.reactions.LCTStex.bounds)
    
    print(str(E0.reactions.EX_ac_e), str(E0.reactions.ACtex))
    
def get_BM_df(current_gene, n_dir, alpha_table, co=True, mono=True, mono_S=True, p=None, checker_suffix=None, E0=None, S0=None,return_sim=False, 
              ko=False, carbon_source_val=5e-3, add_nutrient_val=[100], initial_pop=initial_pop, obj_style='MAX_OBJECTIVE_MIN_TOTAL', no_monoculture=2,
              reduce_SG_cycle=False, Smono_carbon_source='bulk_ac_e'):
                  
    # correction_factor = {'ac': ac_scale, 'gal': gal_scale, 'co2': co2_scale}
    # correction_factor = {k:v for k,v in correction_factor.items() if v is not None} if w else {}
    if p is None:
        sys.exit('Need to provide COMETS params')
    
    def get_p_short(p):
        p_short = copy.deepcopy(p)
        if p.get_param("maxCycles") > 200:
            p_short.set_param("maxCycles", 200) 
        return p_short
    
    genes=convert_arg_to_list(current_gene) 
    base = f"/panfs/jay/groups/0/harcombe/wong0755/comets_RPS/rep_{n_dir}/"
    # base = f"../comets_RPS/rep_{n_dir}/"
    if not os.path.exists(base):
        os.mkdir(base)  
        
    with E0 as m_E0, S0 as m_S0:
        M0 = [m_E0, m_S0] 
        if not ('Normal' in genes):
            alphas = iter_species(M0, get_alphas_from_tab, genes=genes, alpha_table=alpha_table)
            zip_arg = zip(M0, alphas)
            iter_species(zip_arg, alter_Sij,genes=genes, ko=ko)
        
        E_model,S_model = iter_species(M0, create_c, initial_pop=initial_pop, obj_style=obj_style) # create comet object with altered stiochiometry
        
        co_layout, *mono_layout = create_layout_object(E_model, S_model, carbon_source_val=carbon_source_val, add_nutrient_val=add_nutrient_val,
                                                        co=co, mono=mono, mono_S=mono_S, Smono_carbon_source=Smono_carbon_source)

        if co:
            print(p.all_params['maxCycles'], 'co_p')
            co_df, co_sim = sim_cultures([E_model, S_model], co_layout, base=base, p=p)
            full_df = co_df.add_suffix(f'_{genes}_coculture')
            
        else: 
            co_df, co_sim = None, None
            full_df = pd.DataFrame()
        
        if mono: 
            if mono_S is False:
                no_monoculture = 1
            zip_arg = zip([E_model, S_model], mono_layout[:no_monoculture]) # only E_model layout if mono_S is False
            
            # change maxCycle of param for monoculture to 200
            p_short = get_p_short(p) if ((len(current_gene) == 1) and reduce_SG_cycle) else p # only single gene inhibition 
            mono_output = iter_species(zip_arg, sim_cultures, base=base, p=p_short)
            mono_df, mono_sim = unpack_output_object(mono_output) 
            for new_df in mono_df:
                full_df = pd.concat([full_df, new_df.add_suffix(f'_{genes}_monoculture')],axis=1)
        else:
            mono_sim = None
        # out_dict = extract_dfs(mono_sim, co_sim, flux_correction_factor=correction_factor)
        out_dict = extract_dfs(mono_sim, co_sim, flux_correction_factor={})
        genes_str = '.'.join(genes)
        
        # adjust for checker board
        if checker_suffix:
            full_df = full_df.add_suffix(checker_suffix) # ((G1_lv, G2_lv)
            genes_str = genes_str + checker_suffix
        
        out_dict.update( {'Gene_inhibition': genes_str}) # for DG
        print(genes_str)
        
        full_df.columns = rename_columns(full_df)
        # print(full_df.columns)

    return (full_df, out_dict) if not return_sim else (full_df, out_dict, co_sim)

def rename_columns(df):
    df.columns = [re.sub('S0_ac_','S0.ac_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub('S0_gal_','S0.gal_', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub(',','.',
           re.sub('\'|\(|\)| |\[|\]','',ele)) # ('gene1', 'gene2') -> gene1.gene2
           for ele in df.columns]
    return(df.columns)
    
def gene_index_culture_col_df(analysis_df): 
    analysis_df['Gene_inhibition'] =  ['.'.join(map(str, convert_arg_to_list(l))) for l in analysis_df.Gene_inhibition] # SG, DG, checkerboard g1.g2 form
    analysis_df = analysis_df.set_index('Gene_inhibition')
    return analysis_df

# %%
