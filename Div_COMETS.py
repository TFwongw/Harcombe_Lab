#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

no_monoculture=2
initial_pop = 1e-8 # initial biomass in comets simulation
n_combos = 50 # first n gene combo 
genes = ('folA','folP')
alphas = [[1,2],[3,4],[5,6]]
n_processor = 26 # number of process run in parallel

os.environ["GUROBI_COMETS_HOME"] = os.environ["GUROBI_HOME"] # Gurobi path  

E0 = cobra.io.read_sbml_model("./models/iML1515_E0.xml")
S0 = cobra.io.read_sbml_model("./models/STM_v1_0_S0.xml")

p = c.params()
p.set_param("defaultKm", 0.00001) # M 
p.set_param("defaultVmax", 10) #mmol/gDw/hr
p.set_param("maxCycles", 450)
# p.set_param("maxCycles", 20)
p.set_param("timeStep", 1) 
p.set_param('writeFluxLog', True)
p.set_param('writeMediaLog', True)

E0.id = 'E0'
S0.id = 'S0.ac'
# FVA_bounds = pd.read_csv('./Data/FVA_bounds_full.csv', index_col= 0)
 
# # function for iterate 3 species
def create_common_media(Species, carbon_source = "lcts_e", carbon_source_val = 1e-1, 
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
    
def sim_cultures(model, layout, p=p, base = None, genes=''): # non_functional argument model
    # separate function into mono & coculture to prevent using wrong layer
    if type(layout) is list:
        layout = layout[0] # layout object stored inside list of one element unpacking from iter_species
    media_concentration = layout.media[['init_amount','metabolite']]
    # print(model, media_concentration[media_concentration["metabolite"]=="lcts_e"])
    
    sim = c.comets(layout, p)    
    sim.working_dir = base
    try:
        sim.run()
        # if modify_bool==True:
        #     sim.total_biomass.rename(columns={'S0.ac':'S0.glc'},inplace=True)
    except:
        logging.exception(f"{sim.run_output}")
        print(f"{sim.run_output}")
    return sim.total_biomass.set_index('cycle'), sim


def convert_arg_to_list(arg):
    if type(arg) in [pd.Series, pd.Index]:
        arg = list(arg) 
    elif type(arg) not in [list, tuple, set]:
        arg = [arg] # differ from list(arg) -> conversion of str
    return arg

def create_c(model, initial_pop=initial_pop):
    model = c.model(model)
    model.open_exchanges()
    model.initial_pop = [0, 0, initial_pop]
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
alpha_table.columns = ['E0', 'S0.ac', 'S0.glc']

def get_alpha_steps(SG, steps=np.arange(0.2,2,0.2)):
    result_df = pd.DataFrame([])
    for scaling in steps:
        temp_df = alpha_table.loc[[SG]]*scaling
        temp_df.index = temp_df.index + f'_{round(scaling,3)}'
        result_df = pd.concat([result_df, temp_df])
    return result_df

def get_alphas_from_tab(model, genes: list, alpha_table = alpha_table):
    genes = convert_arg_to_list(genes)
    alphas = [alpha_table[model.id][gene] for gene in genes]
    return alphas 

def scale_reaction(model, reaction_id, alpha, direction='forward'):
    if direction == 'forward':
        model.reactions.get_by_id(reaction_id).lower_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v >= 0)    # obetain 'end product'
    else:  # reverse reaction
        model.reactions.get_by_id(reaction_id).upper_bound = 0
        dict_product = dict((k, v) for k, v in model.reactions.get_by_id(reaction_id).metabolites.items() if v <= 0)
    
    for product, unit in dict_product.items():  # scale corresponding metabolite units involved in the reaction
        model.reactions.get_by_id(reaction_id).add_metabolites({product: -unit*(1-1/alpha)}) # only change unit of unidirection metabolite
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

def alter_Sij(model, alphas = 1, genes = 'folA'): 
    # get objective value for corresponding alpha
    alphas= convert_arg_to_list(alphas[0]) if type(alphas) is list and len(alphas)==1 else convert_arg_to_list(alphas)  # unlist one layer from zip comprehension 
    genes =  convert_arg_to_list(genes)
    print(model, genes,alphas)
    
    genes_dict = {gene: alpha for gene, alpha in zip(genes, alphas)}
    genes_sorted = sorted(genes_dict.items(), key=lambda x:x[1], reverse=True) #sort by magnitude of alpha
    rct_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    for current_gene, alpha in genes_sorted:
        current_gene = current_gene.split('_')[0] # for step_alpha
        for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions: 
            if (rct.id not in rct_ids):
#                 print(current_gene, alpha)
                rct_ids.extend(separate_reaction(model, rct.id, alpha))# copy of reaction, forward_terminal_change = True
    return(rct_ids) 

def create_layout_object(E_model, S_model):
    co_layout = create_common_media([E_model, S_model], carbon_source='lcts_e')
    E0_layout = create_common_media([E_model], carbon_source='lcts_e', add_nutrient='met__L_e', add_nutrient_val=[100])
    S0_ac_layout = create_common_media([S_model], carbon_source='ac_e')
    S0_glc_layout = create_common_media([S_model], carbon_source='glc__D_e')
    return [co_layout, E0_layout, S0_ac_layout, S0_glc_layout]

def unpack_output_object(output):
    df, *sim_object = zip(*output)
    return df, sim_object
# df, sim_object = unpack_output_object(mono_output)

def extract_dfs_from_sim_object(sim_object):
    species_name = [ele for ele in sim_object.total_biomass.columns[1:]]
    if len(species_name)>1:
        culture = 'coculture' 
        out_dict = {f'{culture}_media' : sim_object.media}
    else: 
        culture = f'monoculture'
        out_dict = {f'{species_name[0]}_{culture}_media' : sim_object.media}
    for species in species_name:
         out_dict[f'{species}_{culture}_flux'] = sim_object.fluxes_by_species[f'{species}']
    return out_dict
        
# unpack_output_object
def extract_dfs(mono_sims, co_sim):
    out_dict = extract_dfs_from_sim_object(co_sim)
    if mono_sims:
        for sim_object in mono_sims[0]:
            out_dict.update(extract_dfs_from_sim_object(sim_object))
    out_dict = {k: v.to_json() for k,v in out_dict.items()}
    return out_dict

def get_BM_df(genes,n_dir='',alpha_table=alpha_table,mono=True, E0=E0, S0=S0):
    genes=convert_arg_to_list(genes)
    print(genes)
    base = f"/panfs/jay/groups/0/harcombe/wong0755/comets_RPS/rep_{n_dir}/"
    if not os.path.exists(base):
        os.mkdir(base)
    with E0 as m_E0, S0 as m_S0:
        M0 = [m_E0, m_S0] 
        if not ('Normal' in genes):
            alphas = iter_species(M0, get_alphas_from_tab, genes=genes, alpha_table=alpha_table)
            zip_arg = zip(M0, alphas)
            iter_species(zip_arg, alter_Sij,genes=genes)
    
        E_model,S_model = iter_species(M0, create_c) # create comet object with altered stiochiometry
        
        co_layout, *mono_layout = create_layout_object(E_model, S_model)
        zip_arg = zip([E_model, S_model, S_model], mono_layout[:no_monoculture]) 
        
        co_df, co_sim = sim_cultures([E_model, S_model], co_layout, base=base)
        full_df = co_df.add_suffix(f'_{genes}_coculture')
        if mono:
            mono_output = iter_species(zip_arg, sim_cultures, base=base)
            mono_df, mono_sim = unpack_output_object(mono_output) 
            for new_df in mono_df:
                full_df = pd.concat([full_df, new_df.add_suffix(f'_{genes}_monoculture')],axis=1)
            out_dict = extract_dfs(mono_sim, co_sim)
        else:
            out_dict = extract_dfs(None, co_sim)
            
        out_dict.update( {'Gene_inhibition': '.'.join(genes)}) # for DG
    return full_df, out_dict   

def rename_columns(df):
    df.columns = [re.sub('S0_','S0.', ele) for ele in df] # S0_ac -> S0.ac
    df.columns = [re.sub(',','.',
           re.sub('\'|\(|\)| |\[|\]','',ele)) # ('gene1', 'gene2') -> gene1.gene2
           for ele in df.columns]
    return(df.columns)
    
def gene_index_culture_col_df(analysis_df): 
    analysis_df['Gene_inhibition'] =  ['.'.join(map(str, l)) for l in analysis_df.Gene_inhibition]
    analysis_df = analysis_df.set_index('Gene_inhibition')
    return analysis_df

def generate_csv(gene_combos: list, filename: str, alpha_table=alpha_table, mono=True):
    zipped_arg = [[ele, i%n_processor, alpha_table, mono] for i, ele in enumerate(convert_arg_to_list(gene_combos))] # also zip alpha_table, 'monoculture' to pass to get_BM_df for discrete concentration 

    # Multiprocessing double gene
    with multiprocessing.Pool(n_processor) as pool:
        result_df_list, result_dict_list = zip(*pool.starmap(get_BM_df, zipped_arg)) 
    result_df = pd.concat(result_df_list,axis = 1)
    result_df.columns = rename_columns(result_df)
    result_df.to_csv(f'./Data/{filename}.csv')
    
    analysis_df = gene_index_culture_col_df(pd.DataFrame(result_dict_list))
    analysis_df.to_json(f"./Data/flanalysis_{filename}.json") 

    return result_df, analysis_df 

# steps_thrB = get_alpha_steps('thrB', np.arange(0.2,1.2,0.2))

# result_df, analysis_df = generate_csv(teps_thrB.index, 'ttest', steps_thrB, mono=False)

    
# Run
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, filename="logfile", filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logging.info("Trial")
    
    start = time()
    
    gene_combos = pd.read_csv('./Data/GeneCombos.csv',header=None)[0][:n_combos]
    gene_combos = list(ast.literal_eval(ele) for ele in gene_combos)
    gene50 = list(alpha_table.index)
    gene50.extend(['Normal'])
    _ = generate_csv(gene50, 'BM_SG1')
    _ =  generate_csv(gene_combos, 'BM_DG1') 


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

