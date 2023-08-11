#!/usr/bin/env python
# coding: utf-8

import cobra 
import cometspy as c
import multiprocessing
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools
import ast
import re
import logging
from time import time
from itertools import chain

n_processor = 4 # number of process run in parallel
initial_biomass = 1e-4 # initial biomass in comets simulation
n_combos = 3 # first n gene combo
# n_set = 4
genes = ('folA','folP')
alphas = [[1,2],[3,4],[5,6]]

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Trial")

os.environ["GUROBI_COMETS_HOME"] = os.environ["GUROBI_HOME"] # Gurobi path 
E0 = cobra.io.read_sbml_model("./models/iML1515_E0.xml")
S0 = cobra.io.read_sbml_model("./models/STM_v1_0_S0.xml")

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 200)
p.set_param("timeStep", 0.2) 

E0.id = 'E0'
S0.id = 'S0_ac'
FVA_bounds = pd.read_csv('./Data/FVA_bounds_full.csv', index_col= 0)

# # function for iterate 3 species
def create_common_media(Species, carbon_source = "lcts_e", carbon_source_val = 10, 
                        nutrients_val = 100, add_nutrient = '', add_nutrient_val = [100]):
    # carbon_source = 'lcts_e', 'ac_e' or 'glc__D_e'
    # add_nutrient = 'met__L_e' for E.coli monoculture
    
    # convert into list for enumeration 
    if type(add_nutrient)  == str:
        add_nutrient = [add_nutrient]
    if type(add_nutrient_val)  == int:
        add_nutrient_val = [add_nutrient_val]
        
    l = c.layout(Species)
    l.obj_style = 'MAX_OBJECTIVE_MIN_TOTAL' # parsimonious FBA
    
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
def sim_cultures(model,suffix='', p=p, base = None):
    modify_bool = False # naming of S0_glc
#     if model is list: #coculture
    if type(model) is list: #coculture
        model_list = model
        l = create_common_media(model_list, carbon_source='lcts_e')
    else:
        suffix = suffix[0]
        # print(suffix, model.id)
        if model.id == 'E0' or model.id == 'iML1515':
            l = create_common_media([model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=10)
        elif suffix in model.id:
            l = create_common_media([model], carbon_source='ac_e')
        else:
            modify_bool = True
            l = create_common_media([model], carbon_source='glc__D_e')
    sim = c.comets(l, p)
    sim.working_dir = base
    try:
        sim.run()
        sim.total_biomass.rename(columns={'S0_ac':'S0.ac'},inplace=True)
        
        if modify_bool==True:
            sim.total_biomass.rename(columns={'S0.ac':'S0.glc'},inplace=True)
    except:
        logging.exception(f"{sim.run_output}")
        print(f"{sim.run_output}")
    return(sim.total_biomass.set_index('cycle'))

def separate_reaction(model, reaction_id, forward_terminal_change = True):
    """
    decompose reversible reaction into forward and backward
    forward_terminal_change : classification of whether the forward reaction aligns with sign of FVA  
    if direction aligns, terminal metabolite is to be scaled -> original reaction
    if direction is opposite, quantity of both substrate and product in reaction is kept intact -> reaction_v1
    """
    (lb, ub) = model.reactions.get_by_id(reaction_id).bounds 
    rct_ids = [reaction_id] #? 
    if(lb < 0 and ub !=0): # only perform if reaction is bidirectional
        intact_reaction = model.reactions.get_by_id(reaction_id).copy() # copy of target reaction

        if (forward_terminal_change ==True): # forward reaction matches sign of FVA bound 
            model.reactions.get_by_id(reaction_id).bounds = (0,ub) # forward only for the orginal reaction
            intact_reaction.bounds = (lb,0) # backward only for reaction_v1
        else: # forward reaction is opposite of sign of FVA bound 
            model.reactions.get_by_id(reaction_id).bounds = (lb,0) 
            intact_reaction.bounds = (0,ub)

        intact_reaction.id = f'{reaction_id}_v1' # redefine id for copy of reaction
        model.add_reactions([intact_reaction]) # add to model
        rct_ids.append(intact_reaction.id) # list of id for reaction and reaction_v1 
#         print(rct_ids)
    return(rct_ids)

def convert_arg_to_list(arg):
    if type(arg) is not list and type(arg) is not tuple:
        arg = [arg]
    return(arg)

# def create_c(model=E0, initial_pop=[0, 0, 1initial_biomass]):
def create_c(model, initial_pop=[0, 0, initial_biomass]):
    model = c.model(model)
    model.open_exchanges()
    model.initial_pop = initial_pop
    model.obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'
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

alpha_table = pd.read_csv('./Data/alpha_table.csv', index_col='Gene_inhibition')

def get_alphas_from_tab(model, genes: list):
    genes = convert_arg_to_list(genes)
    alphas = [alpha_table[model.id][gene] for gene in genes]
    return alphas 

def alter_Sij(model, alphas = 1, genes = 'folA',FVA_bounds=FVA_bounds): 
    # get objective value for corresponding alpha
    alphas=alphas[0] # unlist one layer from zip comprehension 
    print(model, genes,alphas)
    rct_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    for current_gene, alpha in zip(convert_arg_to_list(genes), convert_arg_to_list(alphas)):
#         print(current_gene,alpha)
        for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions: 
            if (rct.id not in rct_ids):
                if(FVA_bounds[f'FVA_{model.id}'][rct.id]>=0):   # if forward reaction align with FVA, reduce product produced
                    rct_ids.extend(separate_reaction(model, rct.id, forward_terminal_change=True))# copy of reaction, forward_terminal_change = True
                    dict_product = dict((k, v) for k, v in model.reactions.get_by_id(rct.id).metabolites.items() if v >= 0)
                else:  # forward reaction oppose to FVA, reduce substrate consumed
                    rct_ids.extend(separate_reaction(model, rct.id, forward_terminal_change=False))# copy of reaction, forward_terminal_change = True           
                    dict_product = dict((k, v) for k, v in model.reactions.get_by_id(rct.id).metabolites.items() if v <= 0)
                for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
#                     print(rct.id)
                    model.reactions.get_by_id(rct.id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite
    return(rct_ids) 
    
no_monoculture=2
def get_BM_df(genes,n_dir=''):
    # m_E0, m_S0 = E0.copy(), S0.copy()
    print(n_dir)
    base = f"/panfs/jay/groups/0/harcombe/wong0755/comets_RPS/rep_{n_dir}/"
    try:
        os.mkdir(base)
    except:
        print("Base directory exists")
    with E0 as m_E0, S0 as m_S0:

        M0 = [m_E0, m_S0] 
        alphas = iter_species(M0, get_alphas_from_tab, genes=genes)
        zip_arg = zip(M0, alphas)
        iter_species(zip_arg, alter_Sij,genes=genes)

        E_model,S_model = iter_species(M0, create_c) # create comet object with altered stiochiometry
        zip_arg = zip([E_model, S_model, S_model],['','ac','glc'][:no_monoculture])
        mono_df = (pd.concat(iter_species(zip_arg, sim_cultures, base=base), axis=1) # df
                            .add_suffix(f'_{genes}_monoculture'))
        co_df = sim_cultures([E_model,S_model], base=base).add_suffix(f'_{genes}_coculture')
        full_df = pd.concat([co_df, mono_df],axis=1)
    return(full_df) 

def rename_columns(df):
    df.columns = [re.sub('S0_','S0.', ele) for ele in df]
    df.columns = [re.sub(',','.',
           re.sub('\'|\(|\)| |\[|\]','',ele)) 
     for ele in df.columns]
    
    return(df.columns)
    

start = time()

gene_combos = pd.read_csv('./Data/GeneCombos.csv',header=None)[0][:n_combos]
gene_combos = list(ast.literal_eval(ele) for ele in gene_combos)
gene_combos_set = list(set(chain(*gene_combos)))

def generate_csv(gene_combos: list, filename: str):
    zipped_arg = [[ele, i%n_processor] for i, ele in enumerate(gene_combos)]

    # Multiprocessing double gene
    with multiprocessing.Pool(n_processor) as pool:
        results = pool.starmap(get_BM_df, zipped_arg)

    result_df = pd.concat(results,axis = 1)
    result_df.columns = rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename}.csv')
    return(result_df)

# generate_csv(gene_combos, gene_combos_BM_DG)
# generate_csv(gene_combos_set, gene_combos_BM_SG)
generate_csv(gene_combos, 'g_DG')
generate_csv(gene_combos_set, 'g_SG')

# zipped_arg = [[ele, i%n_processor] for i, ele in enumerate(gene_combos)]
# start = time()

# # Multiprocessing double gene
# with multiprocessing.Pool(n_processor) as pool:
#     results = pool.starmap(get_BM_df, zipped_arg)

# result_df = pd.concat(results,axis = 1)
# result_df.columns = rename_columns(result_df)
# result_df.to_csv('./Data/gene_combos_BM_DG.csv')

    
# Single gene 
# gene_combos_set = list(set(chain(*gene_combos)))
# for gene_pair in gene_combos_set:
#     results.append(get_BM_df(gene_pair))
    
# result_df = pd.concat(results,axis = 1)
# result_df.to_csv('./Data/gene_combos_BM_SG.csv')



end = time() 
print('Time Elapsed: ', end-start)
