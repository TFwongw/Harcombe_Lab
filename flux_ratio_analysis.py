import re
import pandas as pd
from setup import get_component
import itertools

if 'E0' not in globals():
    E0, S0, all_components = None, None,None

def get_total_carbon(rxn, model, all_components, is_substrate=True,  carbon_detail=False, atp_only=False, atp_exclude=True): # TODO: exclude EX_met__L_e
    def get_n_carbon(formula):
        if formula in ['CO2', 'CH1O2']:
            return 1
        value = re.findall(r'C(\d+)', formula)
        numeric_value = int(value[0]) if value else 0
        return(numeric_value)
    
    def eval_include_in_dict(metabolite, factor, is_substrate=True):
        if atp_exclude and metabolite.id == 'atp_c':
            return False
        if rxn.id == 'BIOMASS_Ec_iML1515_core_75p37M' and atp_only and metabolite.id != 'atp_c':
            return False
        return (factor<0) == is_substrate and get_n_carbon(metabolite.formula)!=0
    
    if isinstance(rxn, str):
        rxn = get_component(model, rxn, all_components)
    d = {m.id: {'formula': m.formula,
                'formula_carbon': get_n_carbon(m.formula),
                'factor': factor,
                'carbon': factor * get_n_carbon(m.formula),}
                for m, factor in rxn.metabolites.items() if eval_include_in_dict(m, factor, is_substrate)}
    carbon_only = {k:v['carbon'] for k,v in d.items()}
    if carbon_detail:
        return carbon_only
    carbon_total = sum(carbon_only.values())
    return carbon_total

def get_metabolite_summary(df, metabolite, index=0, top_n=10, model=E0, all_components=all_components, concat=False): # disregard factor, because of reversible
    # already separated into production_flux and consumption_flux 
    
    if isinstance(df, pd.DataFrame):
        if len(df)>1:
            solution = df.loc[index] if isinstance(index, str) else df.iloc[index]
        else: 
            solution = df.iloc[0]
    elif isinstance(df, pd.Series):
        solution = df
    if isinstance(metabolite, str):
        metabolite = get_component(model, metabolite, all_components)
    
    r_in_sol = [r for r in metabolite.reactions if r.id in solution.index]
    flux = pd.DataFrame( # code modified from cobra
        data=[
            (
                # r.id,
                solution[r.id],
                r.get_coefficient(metabolite.id),
            )
            for r in r_in_sol
        ],
        columns=["flux", "factor"],
        # columns=["reaction", "flux", "factor"],
        index=[r.id for r in r_in_sol],
    )
    # Scale fluxes by stoichiometric coefficient. (positive is production, negative is consumption)
    flux["flux"] *= flux["factor"]
    flux.index.name = "reaction"
    flux['Speices'] = model.id
    flux['Gene_inhibition'] = solution.name
    production_flux = flux.query("flux > 0").copy()
    consumption_flux = flux.query("flux < 0").copy()
    production, consumption = production_flux["flux"].dropna().abs(), consumption_flux["flux"].dropna().abs()
    production_flux["percent"] = (production / production.sum()).apply("{:.2%}".format)
    consumption_flux["percent"] = 1*(consumption / consumption.sum()).apply("-{:.2%}".format) # negative sign as consumption

    # return consumption_flux.sort_values('flux', ascending=True)
    out_list = (production_flux.sort_values('flux', ascending=False)[:top_n],
            consumption_flux.sort_values('flux', ascending=True)[:top_n])
    return pd.concat(out_list) if concat else out_list

def get_carbon_allocation_summary(s, model, all_components, detailed=False, carbon_dict= {}, format=True): # accept series as input and return allocation df   
    for rxn, flux_value in s.filter(regex='EX_.*e$|BIOMASS_Ec_iML1515_core_75p37M$').dropna().items():
        
        carbon_per_flux = get_total_carbon(rxn, is_substrate=True, model=model, all_components=all_components, atp_exclude=True)

        carbon_dict.update({rxn: {
                            'carbon_per_flux': carbon_per_flux,
                            'flux_quantity': flux_value,
                            'carbon_exchange': float(flux_value * carbon_per_flux),
                            }})
    if detailed:
        return pd.DataFrame(carbon_dict).T

    carbon_allocation = pd.DataFrame({rxn: {'total_carbon': v['carbon_exchange']}
                        for rxn, v in carbon_dict.items() if abs(v['carbon_exchange'])>1e-4}).T
    normalize_carbon = carbon_allocation.loc['EX_lcts_e','total_carbon']
    carbon_allocation['percent'] = (carbon_allocation.total_carbon/normalize_carbon)
    waste_portion = -1-carbon_allocation.query('percent<0').percent.sum()
    if abs(waste_portion)<0.01: # if waste portion is too small, ignore it
        waste_portion = 0

    waste_row = pd.DataFrame([-1*waste_portion*normalize_carbon, waste_portion] # inclusion of waste product
                             ,columns=['Waste'], index=carbon_allocation.columns).T
    
    carbon_allocation = pd.concat([carbon_allocation, waste_row])
    if format:
        carbon_allocation['percent'] = carbon_allocation['percent'].apply(lambda x: f'{x:.2%}')
    carbon_allocation['Gene_inhibition'] = s.name if isinstance(s.name, str) else s.name[0]
    carbon_allocation.index.name = 'reaction'
    return carbon_allocation

def get_syn_df(gr_XG, model, gr_DG):
    model_id = model.id
    additive_threshold = .01
    diff_bins = [-10,-1*additive_threshold,1*additive_threshold,10]
    gr_bins = [-1,.2,1,10]
    
    if gr_XG.index.name == 'gene_inhibition':
        gr_XG.index.name = 'Gene_inhibition'
    syn_df  = pd.DataFrame([], index=gr_DG.index)
    syn_df['Predicted_growth_rate'] = gr_XG[f'Predicted_additive_effect_{model.id}_coculture']
    syn_df['Observed_growth_rate'] = gr_XG[f'{model.id}_coculture']
    syn_df['P_O'] = syn_df.Predicted_growth_rate - syn_df.Observed_growth_rate
    syn_df['Drug_comb_effect'] = pd.cut(syn_df['P_O'], bins=diff_bins, labels=['Antagonistic', 'Additive', 'Synergistic'])
    
    syn_df['PGR_bin'] = pd.cut(syn_df['Predicted_growth_rate'], bins=gr_bins, labels=['Low', 'Normal', 'High'])
    syn_df['OGR_bin'] = pd.cut(syn_df['Observed_growth_rate'], bins=gr_bins, labels=['Low', 'Normal', 'High'])
    
    syn_df.loc[(syn_df.Predicted_growth_rate < 1.5e-8) & (syn_df.Observed_growth_rate < 1.5e-8),'P_O'] = 0
    syn_df['gene_sort'] = (syn_df['P_O'] < 0)*10 + abs(syn_df['P_O'])
    syn_df['Species'] = model.id
    return syn_df

def get_antagonistic_df(syn_df):
    antagonistic_list = syn_df.loc[syn_df['P_O']<0].query("P_O<-0.01").index 
    # ? func get pwy col
    potential_pwy = (gene_combo_pathway
                                      .loc[antagonistic_list,'Pathway']
                                      .str.split(' \+ ') # series string need \+
                                      .apply(lambda x: sorted(x))) 
    return pd.DataFrame(potential_pwy)

def remove_nan_from(x: pd.Series):
    return x.apply(lambda x: sorted(list(itertools.compress(x,[ele not in [None, np.nan] for ele in x]))))
  
def get_p_o_df(E0, S0, gr_DG):
    single_pathway_df = pd.read_csv('./Data/single_pathway_df.csv')
    gcomb_single_pathway_df = single_pathway_df.query('XG=="DG"') # select only DG
    p_o = pd.concat([get_syn_df(gr_DG, E0, gr_DG), get_syn_df(gr_DG, S0, gr_DG)])
    p_o = p_o.merge(gcomb_single_pathway_df.drop('XG', axis=1), left_index=True, right_on='Gene_inhibition', how='outer')
    return p_o.set_index('Gene_inhibition')

# functions for complete flux df
def add_SG_to_p_o(desired_cycle, p_o_full):
    SG_cycle, _ = get_XG_cycle_from(desired_cycle)
    SG_list = list(SG_cycle.index.unique())
    SG_list.remove('Normal')
    SG_empty = (pd.DataFrame(index=SG_list, columns=p_o_full.columns))
    SG_empty['Species'] = [['E0', 'S0']]*len(SG_empty)
    SG_empty = SG_empty.explode(['Species'])
    return pd.concat([p_o_full, SG_empty])

def get_full_df(desired_cycle, p_o_full, additional_df_list): # Missing metabolite
#     df_list = [add_SG_to_p_o(), desired_cycle, end_BM, pwy_rct_df, flux_compare_df]
    p_o_w_SG = add_SG_to_p_o(desired_cycle, p_o_full)
    # additional_df_list = end_BM, pwy_rxn_df, flux_compare_df
    df_list = [p_o_w_SG, desired_cycle.query('culture=="coculture"')]
    df_list.extend(additional_df_list)
    merged_df = df_list.pop(0)
    for i, next_df in enumerate(df_list):
        manual_SI = False

        if 'Species' not in next_df:
            manual_SI = True
        if 'culture' in next_df:
            print(set(next_df.culture))
            next_df = next_df.query('culture=="coculture"')
            print(set(next_df.culture))
            next_df = next_df.drop('culture', axis=1)
        merged_df = (set_GI_SP_as_MI(merged_df)
                     .join(set_GI_SP_as_MI(next_df), how='left')) # left join->DG only
        # inner_join SG info from desired cycle & alpha & flux
    return merged_df.reset_index(level='Species')