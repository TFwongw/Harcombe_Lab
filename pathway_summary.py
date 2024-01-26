from setup import convert_arg_to_list, get_component, get_gene_id
from collections import defaultdict
import pandas as pd
import itertools
import numpy as np

if 'alpha_table' not in globals():
    print('set alpha table None')
    alpha_table = None

label_df = pd.read_excel(open('./Data/iML1515_GP.xlsx', 'rb'), # ? S0 rct pathway label necessary?
              sheet_name='iML1515_GP', index_col=0)

def record_rct(model, current_gene):
    return [rct.id for rct in model.genes.get_by_id(get_gene_id(model, current_gene)).reactions]

def pwy_dict_to_df(pwy_dict):
    df = pd.DataFrame.from_dict({k: ' + '.join(v) for k, v in pwy_dict.items()},
#     df = pd.DataFrame.from_dict({k: '+'.join(v) for k, v in pwy_dict.items()},
                       orient='index', columns=['Pathway'])
    df.index.name = 'Gene_inhibition'
    return df

def convert_dict_to_long_df(d: dict, cols: list):
    gcomb_pwy_list = list()
    for gcomb, pwys in d.items():
        for pwy in pwys:
            gcomb_pwy_list.append([gcomb, pwy])
#     index_col = [col for col in cols if col in 'Gene_inhibition|Reaction']
    df = pd.DataFrame.from_records(gcomb_pwy_list, columns=cols)
    return df

def get_SG_pwy_dict(model, alpha_table=alpha_table):
    gene_rct_dict = defaultdict(list)
    for current_gene in alpha_table.index:
        gene_rct_dict[f'{current_gene}'] = record_rct(model, current_gene)
    
    gene_pwy_dict = defaultdict(list)
    for current_gene,rcts in gene_rct_dict.items():
        for current_reaction in rcts:
            try:
                label = label_df.loc[f'{current_reaction}', 'm_subsystem']
                label = label if type(label) == str else label[0]
                gene_pwy_dict[f'{current_gene}'].append(label)
            except:
                print(current_gene, f'reaction {current_reaction} not exist in list')
        gene_pwy_dict[f'{current_gene}'] = set(convert_arg_to_list(gene_pwy_dict[f'{current_gene}'])) | set(gene_pwy_dict[f'{current_gene}'])
    return gene_pwy_dict

def get_gene_pathway_df(model, alpha_table=alpha_table):
    pwy_dict = get_SG_pwy_dict(model, alpha_table=alpha_table)
    DG_pwy_dict = get_DG_pwy_dict(pwy_dict)
    
    pwy_dict.update(DG_pwy_dict)
    
    gene_pathway_df = convert_dict_to_long_df(pwy_dict, ['Gene_inhibition', 'Pathway'])
    gene_pathway_df['XG'] = gene_pathway_df.Gene_inhibition.str.split('.').apply(
        lambda x: 'SG' if len(x)==1 else 'DG') 
    return gene_pathway_df

def antagonistic_count(genes:list):
    antagonistic_set = ['dadX.rffG', 'dadX.pyrD', 'acnB.thrB', 'acnB.thrB', 'acnB.thrB',
       'aroA.argD', 'aroA.argD', 'aroA.argD', 'dadX.guaB', 'aroA.dapB',
       'aroA.guaB']
    return list((set(genes).intersection(antagonistic_set)))

def get_ratio_col(gene_count_df):
    DG = gene_count_df.query('XG=="DG"').reset_index(drop=True).copy()
    DG['Antagonistic_gcomb'] = DG.Gene_list.apply(lambda x: antagonistic_count(x.split(', ')))
#     DG['Antagonistic_ratio'] = DG.Antagonistic_gcomb.str.split(', ').map(len)
    DG['Antagonistic_ratio'] = pd.Series(map(len, DG.Antagonistic_gcomb.tolist()))/DG.Gene_count
    DG['Antagonistic_gcomb'] = DG['Antagonistic_gcomb'].apply(lambda x: ', '.join(x))
    
    return DG[['Pathway', 'Antagonistic_ratio', 'Antagonistic_gcomb']]

def get_gene_count_df(gene_pathway_df):

    gene_count_df = (gene_pathway_df.groupby(['Pathway','XG'], as_index=False).agg(
                        Gene_list=('Gene_inhibition', lambda x: ', '.join(x.values)),
                        Gene_count=('Gene_inhibition', 'count')))
    gene_count_df = gene_count_df.merge(get_ratio_col(gene_count_df), on='Pathway')
    (gene_count_df.pivot(index=['Pathway','Antagonistic_ratio', 'Antagonistic_gcomb'], columns=['XG'], values=['Gene_list','Gene_count'])
         .sort_values(by=('Gene_count','DG'), ascending=False)).to_csv('./Data/gene_count.csv')
    return gene_count_df

def get_single_pathway_df(gene_pathway_df):
    # add p_o afterwards
    df = (gene_pathway_df.groupby(['Gene_inhibition','XG'], as_index=False).agg(
        Pathway_count=('Pathway', 'count'),
        Pathway_list = ('Pathway', lambda x: list(x.values))))
#         P = ('Pathway', lambda x: pd.DataFrame(x.values))))
    df['Single_pathway'] = pd.cut(df['Pathway_count'], [0,1,10], labels = ['Single_pathway', 'Multi_pathway'])
    return df.sort_values(['Pathway_count', 'XG'], ascending=True).set_index('Gene_inhibition')

def get_rct_pathway_df():
    pwy_rct_dict = defaultdict(list)
    for index, sub_df in label_df[['m_subsystem']].iterrows():
        pwy_rct_dict[sub_df[0]].append(index)
    rct_pathway_df = convert_dict_to_long_df(pwy_rct_dict,['Pathway','Reaction'])
    return rct_pathway_df.drop_duplicates()

def get_DG_pwy_dict(SG_pwy_dict):
    comb_list = list(pd.read_csv('./Data/gr_Div_DG_Blis_Aug31.csv').gene_inhibition[1:])
    DG_pwy_dict = dict()
    # for DG_pwy in comb_list:
    gene_pathway_row = list()
    for DG_pwy in comb_list:
        gene_pair = DG_pwy.split('.')
        DG_pwy_dict[DG_pwy] = SG_pwy_dict[gene_pair[0]] | SG_pwy_dict[gene_pair[1]] 
    return DG_pwy_dict

# reaction column
def get_SG_pwy_df(model, E0, S0, alpha_table):
    def remove_nan(x: pd.Series): #=function remove_nan_from
        return x.apply(lambda x: sorted(list(itertools.compress(x,[ele not in [None, np.nan] for ele in x]))))

    SG_pwy_df = pd.DataFrame.from_dict(get_SG_pwy_dict(model, alpha_table=alpha_table), orient='index')
    SG_pwy_df['Pathways'] = SG_pwy_df.values.tolist()
    SG_pwy_df = pd.DataFrame(remove_nan(SG_pwy_df.Pathways))
#     SG_pwy_df = pd.DataFrame(SG_pwy_df.Pathways
#                              .apply(lambda x: remove_nan(x)))
    SG_pwy_df['Reactions_E0'] = list(pd.Series(SG_pwy_df.index).apply(lambda x: record_rct(E0,x)))
    SG_pwy_df['Reactions_S0'] = list(pd.Series(SG_pwy_df.index).apply(lambda x: record_rct(S0,x)))    
    SG_pwy_df.index.name='Gene_inhibition'
    return SG_pwy_df

def find_reactions_not_exist_in(model, reactions: list):
    if len(reactions) ==0:
        return reactions
    reactions = convert_arg_to_list(reactions)
    reactions_all = [rct.id for rct in model.reactions]
    not_found_reactions = [rct for rct in reactions if rct not in reactions_all]
    return not_found_reactions

def get_pwy_rxn_df(model, E0, S0, alpha_table, DG_list):
    def get_pwys_from(SG_pwy_df, SG):
        return SG_pwy_df.loc[SG, 'Pathways']
    SG_pwy_df = get_SG_pwy_df(model, E0, S0, alpha_table)
    g_pwy_rct_df = pd.DataFrame(DG_list, columns=['Gene_inhibition'])
    g_pwy_rct_df['SG_list'] =  g_pwy_rct_df.apply(lambda x: x.Gene_inhibition.split('.'), axis=1)
    g_pwy_rct_df['Pathways'] = g_pwy_rct_df.apply(lambda x: list(set(SG_pwy_df.loc[x.SG_list[0],'Pathways']) | set(SG_pwy_df.loc[x.SG_list[1],'Pathways'])), axis=1)
    g_pwy_rct_df['Common_pathway'] = g_pwy_rct_df.apply(lambda x: set(get_pwys_from(SG_pwy_df, x.SG_list[0]))
                                                        .intersection(set(get_pwys_from(SG_pwy_df, x.SG_list[1])))
                                                        ,axis=1)
    g_pwy_rct_df['Reactions_E0'] = g_pwy_rct_df.apply(
        lambda x: list(set(SG_pwy_df.loc[x.SG_list[0],'Reactions_E0']) | set(SG_pwy_df.loc[x.SG_list[1],'Reactions_E0'])), axis=1)
    g_pwy_rct_df['Reactions_S0'] = g_pwy_rct_df.apply(
        lambda x: list(set(SG_pwy_df.loc[x.SG_list[0],'Reactions_S0']) | set(SG_pwy_df.loc[x.SG_list[1],'Reactions_S0'])), axis=1)
    g_pwy_rct_df['Common_Reactions'] = g_pwy_rct_df.apply(lambda x: set(x['Reactions_E0']).intersection(set(x['Reactions_S0'])), axis=1)
    g_pwy_rct_df['Reactions_E0_only'] = g_pwy_rct_df.apply(lambda x: set(x['Reactions_E0'])-set(x['Reactions_S0']), axis=1)
    g_pwy_rct_df['Reactions_S0_only'] = g_pwy_rct_df.apply(lambda x: set(x['Reactions_S0'])-set(x['Reactions_E0']), axis=1)
    g_pwy_rct_df['Reactions_not_in_S0'] = g_pwy_rct_df.apply(lambda x: find_reactions_not_exist_in(S0, x.Reactions_E0_only), axis=1)
    g_pwy_rct_df['Reactions_not_in_E0'] = g_pwy_rct_df.apply(lambda x: find_reactions_not_exist_in(E0, x.Reactions_S0_only), axis=1)
    g_pwy_rct_df = pd.concat([g_pwy_rct_df.set_index('Gene_inhibition'), SG_pwy_df])
    return g_pwy_rct_df

