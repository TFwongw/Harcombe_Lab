#!/usr/bin/env python
# coding: utf-8

# # Initialization 

import cobra
from cobra.io import load_model
import cometspy as c
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import itertools

E0 = cobra.io.read_sbml_model("./models/iML1515_E0.xml")
S0_ac = cobra.io.read_sbml_model("./models/STM_v1_0_S0.xml")
S0_glc = S0_ac.copy()

lactose_medium = {'EX_ca2_e': 10,
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
 'EX_nh4_e': 10,
 'EX_lcts_e': 10}
lactose_met_medium = lactose_medium.copy()
lactose_met_medium["EX_met__L_e"] = 10 
E0.medium = lactose_met_medium 
E0.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" 

def change_medium(model, metabolite, value): # function for setting medium metabolite value
    medium = model.medium
    medium[metabolite] = value
    model.medium = medium

ac_medium = {'EX_ca2_e': 10,
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
 'EX_nh4_e': 10,
 'EX_ac_e': 10}
S0_ac.medium = ac_medium 

S0_glc.medium = ac_medium
change_medium(S0_glc, 'EX_ac_e', 0)
change_medium(S0_glc, 'EX_glc__D_e', 10)
S0_ac.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" 
S0_glc.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" 

E0.id = 'E0'
S0_ac.id = 'S0_ac'
S0_glc.id = 'S0_glc'

import pandas as pd
S2tab = pd.read_excel('./Data/S2Table.xlsx');
S2_Enzyme = list(S2tab['Enzyme'])
essential_genes = ['purU', 'pyrE', 'thrB', 'yrbG', 'folA', 'folP', 'pykF', 'rffG']
FVA_bounds = pd.read_csv('./Data/FVA_bounds_full.csv', index_col= 0)

# remove reactions involved in preexisting gene inhibition candidates
def get_gene_id(model, gene_name):
    for _,i in enumerate(model.genes):
        if(i.name == gene_name):
            return(i.id)

for i, Enzyme in enumerate(S2_Enzyme):
    mark_remove = False
    for _, current_gene in enumerate(essential_genes):
        for _, current_reaction in enumerate(E0.genes.get_by_id(get_gene_id(E0,current_gene)).reactions):           
            if (Enzyme==current_reaction.id):
                mark_remove = True
                print(current_reaction.id)
    if (mark_remove == True):
        S2_Enzyme.remove(Enzyme)            
            
# Only the first gene in case of degeneracy and choose only genes that also exist in S.enterica
# return first common gene name        
def extract_gene_name(gene_str, model):
    if (model.id =='iML1515' or model.id =='E0'):
        offset = 5
        pos = [pos for pos, char in enumerate(gene_str) if char == 'b']
    else:
        offset = 7
        pos = [pos for pos, char in enumerate(gene_str) if char == 'S']

    pos.reverse()
    gene_list = list()
    gene_str_temp = gene_str+' '
    for i,val in enumerate(pos):
        gene_id = gene_str[val:val+offset]
        gene_list.append(model.genes.get_by_id(gene_id).name)
#         print(gene_id)
        gene_str_temp = gene_str_temp[:val] + model.genes.get_by_id(gene_id).name + gene_str_temp[val+offset:]
    gene_list.reverse()
    return(gene_list,gene_str_temp)

def get_common_gene(Enzyme):
    Enzyme = Enzyme.replace('-','_')
    try: 
        gpr_S0 = str(S0_ac.reactions.get_by_id(Enzyme).gpr)
        gpr_E0 = str(E0.reactions.get_by_id(Enzyme).gpr)
        gene_list_S0,_ = extract_gene_name(gpr_S0, S0_ac)
        gene_list_E0,_ = extract_gene_name(gpr_E0, E0)

        desired_gene = next((ele for ele in gene_list_E0 if ele in set(gene_list_S0)), None) # first occurrence
        return(desired_gene, gene_list_E0, gene_list_S0)
    except:
        print(f'{Enzyme} does not co-exist in E.coli and S.enterica') 
        return(None,(),())

desired_genes = list()
none_common_gene_mat = list()
for _,Enzyme in enumerate(S2_Enzyme):
    Enzyme = Enzyme.replace('-','_')
#     print(Enzyme)
    gene_to_append,b,b = get_common_gene(Enzyme)
    if(gene_to_append == None):
        none_common_gene_mat.append(Enzyme)
    else:
        desired_genes.append(gene_to_append)
        
desired_genes = set(desired_genes)
potential_genes = desired_genes | set(essential_genes)   

# Perform gene inhibition
def set_bound(model, current_reaction, bound):
    old_bound = model.reactions.get_by_id(current_reaction).bounds
    if (bound>=0):
        ub = bound
        lb = 0 if old_bound[0]==0 else -bound
    else:
        lb = bound
        ub = 0 if old_bound[1]==0 else -bound 
    model.reactions.get_by_id(current_reaction).bounds = (lb,ub)
    
obj_flux_df = ['Normal']
obj_flux_df.extend([E0.slim_optimize()])
obj_flux_df.extend([S0_ac.slim_optimize()])
obj_flux_df.extend([S0_glc.slim_optimize()])
obj_flux_df = pd.DataFrame(obj_flux_df).T

# obj_flux_df = pd.concat([obj_flux_df, ], axis = 0)
# for _,current_gene in enumerate(potential_gene):
for _,current_gene in enumerate(list(potential_genes)[:]):
    new_row = list([current_gene])
    for _, (model, FVA_min) in enumerate(zip([E0,S0_ac,S0_glc],
                                          [FVA_bounds[['FVA_E0']],
                                          FVA_bounds[['FVA_S0_ac']],
                                          FVA_bounds[['FVA_S0_glc']]])):
        with model:
            for _,current_reaction in enumerate(model.genes.get_by_id(get_gene_id(model, current_gene)).reactions):
                set_bound(model, current_reaction.id, FVA_min.loc[current_reaction.id][0]/2)
            obj_val = model.slim_optimize()
            new_row.extend([obj_val])
    obj_flux_df = pd.concat([obj_flux_df, pd.DataFrame(new_row).T], axis=0)
                
obj_flux_df.columns = ['Gene_inhibition', 'E0', 'S0_ac', 'S0_glc']
obj_flux_df = obj_flux_df.set_index('Gene_inhibition')

# # Functions for getting objective values
alpha,obj_val = 1,0 

def convert_arg_to_list(arg):
    if type(arg) is not list and type(arg) is not tuple:
        arg = [arg]
    return(arg)

def get_summary_df(model = E0, alpha = alpha, obj_val = obj_val, rct_ids = 'DHFR', sol_sol = None): 
#     summarize solutions from optimization and alpha used  
# expect direction opposite with zero flux, to be consistent with FVA bound 
    if type(alpha) != int and type(alpha) != float:
        alpha = str(alpha)
    if sol_sol is not None:
        rct_dir, rct_rev, flux_dir, flux_rev = list(), list(), list(), list()
        for current_rct_id in rct_ids:
            append_str = str(model.reactions.get_by_id(current_rct_id))
            if ('_v1' not in current_rct_id): # direction align with FVA
                rct_dir.extend([append_str])
                try:
                    flux_dir.extend([round(model.reactions.get_by_id(current_rct_id).flux,5)])
                except: 
                    flux_dir.extend([sol_sol[current_rct_id]])
            else: # direction opposite
                rct_rev.extend([append_str])        
                try:
                    flux_rev.extend([round(model.reactions.get_by_id(current_rct_id).flux,5)])
                except:
                    flux_rev.extend([sol_sol[current_rct_id]])
        net_flux = sum(abs(element) for element in set(flux_dir) | set(flux_rev))
        net_flux_I = 'Zero Flux' if net_flux ==0 else 'Net Flux'
        summary_df = pd.DataFrame({f'div_opt_alpha': alpha,
                           f'div_opt_obj_val': obj_val,
                           f'Reactions_id': ', '.join(rct_ids),
                           f'Forward_Reactions_involved': ' '.join(rct_dir),                                    
                           f'Forward_Flux_values':[flux_dir],
                           f'Backward_Reactions_involved': ', '.join(rct_rev),     
                           f'Backward_Flux_values':[flux_rev],
                           f'Net_Flux_Boolean':[net_flux_I],
                           f'Net_Flux':[net_flux]}
                         )
    return(summary_df)

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

def Sij_obj(model = E0, alphas = 1, genes = 'folA', slim_opt = False): 
    # get objective value for corresponding alpha
    rct_ids = list() # store list of id of reaction and reaction_v1 that regulated by the same gene 
    alphas= convert_arg_to_list(alphas[0]) if type(alphas) is list and len(alphas)==1 else convert_arg_to_list(alphas)  # unlist one layer from zip comprehension 
    genes =  convert_arg_to_list(genes)
    
    with model:  # enumerate genes
        genes_dict = {gene: alpha for gene, alpha in zip(genes, alphas)}
        genes_sorted = sorted(genes_dict.items(), key=lambda x:x[1], reverse=True) #sort by magnitude of alpha
        
        for current_gene, alpha in genes_sorted:
            for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions: 
                rct_ids.extend(separate_reaction(model, rct.id, alpha))
                    # without abs, desired metabolites approach 0(Less -, less +)
        if (slim_opt == False):
            sol_sol = model.optimize()    
            obj_val = sol_sol.objective_value     
            summary_df = get_summary_df(model, alphas, obj_val, rct_ids, sol_sol.fluxes)  
            return(alpha, obj_val, summary_df)            
        else: 
#             return(pd.DataFrame([{f'div_opt_obj_val': model.slim_optimize()}]))
            return(model.slim_optimize())        

# functions for searching alpha
def evaluate_alpha(model, search_alpha, current_gene, flux_df_col, target_obj_val, opt_df, alpha_lb , alpha_ub, exp_leap, precision):
    _, obj_val, summary_df = Sij_obj(model, search_alpha, current_gene)

    if(opt_df is None and summary_df['Net_Flux'][0]<1e-8): # reinitialize if the first alpha to search gives zero flux
        search_alpha = 1+2e-3
        _, obj_val, summary_df = Sij_obj(model, search_alpha, current_gene)
        
    halt_condition = (summary_df['Net_Flux'][0]<1e-8 or summary_df['Net_Flux'][0]>100 or 
                      obj_val<target_obj_val or obj_val>1)
    termination_condition = (round(search_alpha,precision)==round(alpha_ub,precision) or
                             round(search_alpha,precision)==round(alpha_lb,precision) or
                            search_alpha>alpha_ub-1)
    # update optimal df
    if((summary_df['Net_Flux'][0]>0 and 
        halt_condition == False) or
       (opt_df is None or 
       summary_df['div_opt_obj_val'][0] > target_obj_val*0.9)):
        opt_df = summary_df

    return(search_alpha, termination_condition, halt_condition, opt_df)

def find_feasible_alpha(model, search_alpha, current_gene, flux_df_col = 'E0', target_obj_val = 0.5,
                        alpha_lb = 1+1e-7, alpha_ub = 1e5, exp_leap = 2, precision = 3, low_alpha_label = False):
    opt_df = None # store optimal summary for return function     
    search_alpha, termination_condition, halt_condition,opt_df = evaluate_alpha(model, search_alpha, current_gene, flux_df_col, target_obj_val, opt_df, alpha_lb , alpha_ub, exp_leap, precision)   
#     print(search_alpha,termination_condition, halt_condition)
    while (termination_condition == False): 
        if (halt_condition == True): # search alpha not accepted, adjust search towards alpha_lb
            # start binary search, in (alpha_lb, search_alpha)
            alpha_lb = alpha_lb
            alpha_ub = search_alpha
            
            search_alpha = max(min(1+(alpha_lb-1)*2,(search_alpha+alpha_lb)/2),1e-5)
            _, termination_condition, halt_condition, summary_df = evaluate_alpha(model, search_alpha, current_gene, flux_df_col, target_obj_val, opt_df, alpha_lb , alpha_ub, exp_leap, precision)
            
        else: # press alpha_ub
            alpha_lb = search_alpha
            
            if(search_alpha*exp_leap <alpha_ub): # exponential search(large step)
                while (search_alpha*exp_leap <alpha_ub and halt_condition == False): 
                    search_alpha *= exp_leap
                    _, termination_condition, halt_condition,summary_df = evaluate_alpha(model, search_alpha, current_gene, flux_df_col, target_obj_val, opt_df, alpha_lb , alpha_ub, exp_leap, precision)
            else: 
                # binary search (small step), adjust towards alpha_ub
                search_alpha = (search_alpha+alpha_ub)/2
                _, termination_condition, halt_condition, summary_df = evaluate_alpha(model, search_alpha, current_gene, flux_df_col, target_obj_val, opt_df, alpha_lb , alpha_ub, exp_leap, precision)
        if(summary_df['Net_Flux'][0]>0 and 
           summary_df['div_opt_obj_val'][0]>0):
            opt_df = summary_df
#     while (np.floor(opt_df['div_opt_alpha'][0]*10**(precision-3))/10**(precision-3)-1 == 0):
#         precision += 1
#         print(precisio)
    print(current_gene, opt_df['div_opt_alpha'][0], opt_df[f'div_opt_obj_val'][0])
#     opt_df['div_opt_alpha'] = np.floor(opt_df['div_opt_alpha'][0]*10**(precision))/10**(precision)
    
    opt_df.columns = [f'{flux_df_col}_'+element for element in opt_df.columns]
    return(opt_df[f'{flux_df_col}_div_opt_alpha'],opt_df[f'{flux_df_col}_div_opt_obj_val'],opt_df) # alpha_feasible, and upper bound of alpha

# ## workflow for generating objective value
list_target_obj_val = [E0.slim_optimize(),
                       S0_ac.slim_optimize(),
                       S0_glc.slim_optimize()]

list_target_obj_val = pd.DataFrame([E0.slim_optimize(),S0_ac.slim_optimize(),S0_glc.slim_optimize()]).T
list_target_obj_val.columns = ['E0','S0_ac','S0_glc']
def add_prefix_div_df(df, flux_df_col): # add prefix corresponding to given flux_df_col
    df.columns = [f'{flux_df_col}_{col}' for col in df.columns]
    
def get_div_obj_df(model, flux_df_col,target_obj_val,  
                   first_n_gene=len(obj_flux_df.index[1:]), alpha_df = None): 
    # generate full biomass dataframe for a single species
    # given alpha or automatic search for optimal with the given objective value
    with model:
        obj_div_df = pd.DataFrame()
        for i, current_gene in enumerate(obj_flux_df.index[1:][:first_n_gene]): # iter gene 
#         for i, current_gene in enumerate(obj_flux_df.index[1:][:1]): # iter gene 
            print(current_gene,'in cal')
            if alpha_df is None:
                if(current_gene=='yrbG'):
                    alpha,obj_value,temp_df = Sij_obj(model, alphas = 1.00001,
                                                      genes= 'yrbG')
                    temp_df.columns = [f'{flux_df_col}_'+element for element in temp_df.columns]
                else:
                    alpha,obj_value,temp_df = find_feasible_alpha(model, 1.02, current_gene, flux_df_col, 
                                                              target_obj_val, precision=8)
            else:
                alpha = alpha_df.loc[current_gene,f'{flux_df_col}_div_opt_alpha'] # get alpha for corresponding gene from alpha_df
                alpha,obj_value,temp_df = Sij_obj(model, alpha, #.pop?
                                              genes = current_gene)

            temp_df['Gene_inhibition'] = current_gene    
            obj_div_df = pd.concat([obj_div_df, temp_df.set_index('Gene_inhibition')],axis=0)
    return(obj_div_df)

def combine_div_obj_df(list_target_obj_val=list_target_obj_val, alpha_df = None):
    # join obj_dfs of three cultures
    # return a full dataframe and list containing the three separate cultures
    obj_div_df = pd.DataFrame()
    E0_df = get_div_obj_df(E0, 'E0', target_obj_val=list_target_obj_val['E0'][0], 
                           alpha_df=alpha_df)
    S0_ac_df = get_div_obj_df(S0_ac, 'S0_ac', target_obj_val=list_target_obj_val['S0_ac'][0],
                              alpha_df=alpha_df)
    S0_glc_df = get_div_obj_df(S0_glc, 'S0_glc', target_obj_val=list_target_obj_val['S0_glc'][0],
                               alpha_df=alpha_df)

    obj_div_df = E0_df.join(S0_ac_df)
    obj_div_df = obj_div_df.join(S0_glc_df)
    
    pathway_list = list() # add pathway column
    for _, current_gene in enumerate(obj_div_df.index):
        pathway_sub = list()
        for _, rct in enumerate(E0.genes.get_by_id(get_gene_id(E0, current_gene)).reactions):
            if rct.id in set(S2tab['Enzyme']):
                pathway_sub.extend([S2tab.query(f"Enzyme=='{rct.id}'")['Pathway'].to_string(index = False)])
        pathway_list.append(', '.join(ele for ele in pathway_sub))

    classify_gene = list() # add gene_source
    for current_gene in obj_div_df.index:
        if current_gene in essential_genes:
            classify_gene.append('Drug Similarity')
        else:
            classify_gene.append('div S2table')
    
    obj_div_df = pd.concat([pd.DataFrame(pathway_list, columns=['Pathway'], index=obj_div_df.index),
                            pd.DataFrame(classify_gene, columns=['Gene_Source'], index=obj_div_df.index),
                            obj_div_df], axis=1)
    return(obj_div_df, [E0_df, S0_ac_df, S0_glc_df])

def rearrange_rows(df,filename): # row order accoding to gene source and E0 biomass
    df = df.set_index('Gene_inhibition')
    df = df.sort_values(by = ['Gene_Source', 'E0_div_opt_obj_val'])
#     df = df.sort_values(by = ['Gene_Source', 'div_opt_obj_val_x'])
    df.to_csv(f'{filename}.csv')
    return(df)


# ## Implementations
obj_div_df_50percent, model_df_50percent = combine_div_obj_df(list_target_obj_val/2)
# obj_div_df_90percent, model_df_90percent = combine_div_obj_df(list_target_obj_val*0.9)

from ast import literal_eval
alpha_table = obj_div_df_50percent[['E0_div_opt_alpha','S0_ac_div_opt_alpha' ,'S0_glc_div_opt_alpha']]
alpha_table.applymap(lambda x: literal_eval(x)[0] ).to_csv('./Data/alpha_table.csv')

