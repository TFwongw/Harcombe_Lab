# %%
# %run setup.ipynb
# !jupyter nbconvert --to script setup.ipynb

# %%
from setup import *
from dataclasses import dataclass, field
from typing import Callable, List, Any, Dict
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import concurrent.futures

@dataclass(kw_only=True)
class AlphaFinderConfig:
    alpha_table = None
    opt_df : pd.DataFrame = None # final output
    alpha_lb : float = 1+1e-7
    alpha_ub : float = 1e5 
    found_alpha : bool = None
    is_new_ub : bool = None
    exp_leap : float = 2
    target_obj_val : dict = None  
    ko_intercepted : bool = None
    ever_eval =  False
    iter_max = 25
    i_iter = 0    

    def __post_init__(self): # not post_inited
        self.is_monoculture = None
        self.out_fun = None
        self.eval_alpha_fun = None
    
    @staticmethod
    def get_next_alpha(search_alpha, alpha_lb, alpha_ub, is_new_ub, exp_leap=2, is_monoculture=True):
        # print('gna', search_alpha, alpha_lb, alpha_ub, is_new_ub)
        if (is_new_ub == False):  # raise lb, search higher alpha
            alpha_lb = search_alpha
            # if ((search_alpha*exp_leap <alpha_ub) and is_monoculture): # exponential search(large step)
            # if ((alpha_ub > 1e4) and (alpha_lb> 1e3)): # when alpha ub not get updated, exponential search(large step)
            if (search_alpha*exp_leap <alpha_ub and is_new_ub == False): 
                search_alpha *= exp_leap
            else: # search alpha not get updated, then binary search 
                search_alpha = (search_alpha+alpha_ub)/2        
        else: # search alpha not accepted, lower ub, lower search_alpha
            # start binary search, in (alpha_lb, search_alpha)
            alpha_ub = search_alpha
            search_alpha = max((search_alpha+alpha_lb)/2      
                               ,1+1e-5) # terminate at 1+1e-5
        # print(search_alpha, alpha_lb, alpha_ub)
        return search_alpha, alpha_lb, alpha_ub
    
    @staticmethod
    def is_found(search_alpha, alpha_lb, alpha_ub, precision):
        if (search_alpha<2) and (precision <3):
            precision = 3
        if alpha_lb > 15:
            precision = 1
        return (round(search_alpha,precision)==round(alpha_ub,precision) or
                                 round(search_alpha,precision)==round(alpha_lb,precision) or
                                (search_alpha>9e4))

    def classify_growth_switch(self, min_attain=.9, max_attain=.3):
        # eval all true
        alpha_range = (self.alpha_ub - self.alpha_lb)
        alpha_range_narrow = (((alpha_range < 0.15) and (self.alpha_ub < 3)) or # circumvent precision <2 in subsequent evaluation
                                ((alpha_range < .3) and (5<self.alpha_lb<10)) or
                                ((alpha_range < 1.3) and (self.alpha_lb > 10)))
        alpha_range_req =  (alpha_range_narrow # 
                            and (self.alpha_lb > 1.01)) # ever update lb and ub
        
        obj_req = (not any(0.3 < val < 0.8 for val in self.trace_obj_val) and # mrdA forever > 1 for low dose
                    (min(self.trace_obj_val) < min_attain) and 
                    (max(self.trace_obj_val) > max_attain))
        
        if alpha_range_narrow and all(val >1 for val in self.trace_obj_val): # small dose all > Normal
            obj_req = True
                
        # ? req more evalluation for alpha inside the bound?
        self.is_growth_switch = (alpha_range_req and obj_req) or (self.alpha_ub < 1.018)
        return self.is_growth_switch
    
    def find_feasible_alpha(self): 
        print(self.current_gene, self.ko_intercepted)
        def eval_continue():
            if self.is_monoculture: # as culture_flag
                return not self.found_alpha
            if self.i_iter>self.iter_max:
                self.early_stop = True
                print('Early stopped coculture')
            return ((not self.found_alpha) or (self.i_iter<2) 
                    or (self.alpha_lb < 1.01 and 1.018< self.alpha_ub)) # force iter -ensure optimal

        if self.ko_intercepted: # req run ko_gr first, otherwise cannot catch
            print(f'Intercept Non-essential: {self.current_gene}, calculating with current alpha')
            return self.out_fun() # sim_culture with current alpha_table already ran
                
        if not self.ever_eval:
            self.eval_alpha_fun() # initial ko is without alpha, only used for identify gr_ko

        stop = not (self.alpha_lb <self.search_alpha< self.alpha_ub)
        
        while eval_continue() and not stop: 
            self.search_alpha, self.alpha_lb, self.alpha_ub = self.get_next_alpha(
                self.search_alpha, self.alpha_lb, self.alpha_ub,  self.is_new_ub, 
                self.exp_leap, is_monoculture=self.is_monoculture)
            self.eval_alpha_fun()
            self.i_iter+=1
        if (not self.is_monoculture):
            print(f'Stopped at iter {self.i_iter}') if (self.i_iter<=self.iter_max or self.i_iter ==2) else print('Success search, end at: ', str(self.i_iter))
        
        return self.out_fun()

# %%
def get_summary_df(model, alpha, obj_val, rct_ids = 'DHFR', sol_sol = None): 
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
    #     rct_rev = fix_string(rct_rev)
    #     rct_dir = fix_string(rct_dir)
        net_flux = sum(abs(element) for element in set(flux_dir) | set(flux_rev))
        net_flux_I = 'Zero Flux' if net_flux ==0 else 'Net Flux'

    #     if type(alpha) is list or type(alpha) is tuple:

        summary_df = pd.DataFrame({f'div_opt_alpha': alpha,
                                   f'div_opt_obj_val': obj_val,
                                   f'FVAdir_Reactions_id': ', '.join(rct_ids),
                                   f'FVAdir_Reactions_involved': ' '.join(rct_dir),                                    
                                   f'FVAdir_Flux_values':[flux_dir],
                                   f'FVAopposite_Reactions_involved': ', '.join(rct_rev),     
                                   f'FVAopposite_Flux_values':[flux_rev],
                                   f'Net_Flux_Boolean':[net_flux_I],
                                   f'Net_Flux':[net_flux]}
                                 )
#     else:
#         summary_df = pd.DataFrame([{'div_opt_obj_val': obj_val}])        
    return(summary_df)

@dataclass(kw_only=True)
class MonocultureAlphaFinder(AlphaFinderConfig):
    model : str
    search_alpha : float
    current_gene : str
    target_obj_val : float = .5
    response_record = None # for iteratice uptate in checkerboard
    opt_df : pd.DataFrame = None 
    exp_leap : int = 2
    precision : int = 2 # precision of alpha
    acceptance_threshold_upper : float = None
    acceptance_threshold_lower : float = None
    trace_obj_val : List = field(default_factory=list)
        
    def __post_init__(self):
        super().__post_init__()
        self.is_monoculture = True
        self.eval_alpha_fun = self.evaluate_alpha_from_biomass 
        self.out_fun = self.modify_opt_df
        if not self.acceptance_threshold_lower:
            self.acceptance_threshold_lower = .85
        if not self.acceptance_threshold_upper:
            self.acceptance_threshold_upper = 1.1
        self.norm_obj_val = self.model.slim_optimize()
        if not self.response_record:       
            self.response_record = {self.current_gene: {self.model: {'ICX': convert_arg_to_list(self.target_obj_val), 'response': {}}}}

    def fill_response_record(self, obj_val, is_lowest_abs_diff=False):
        obj_val_interval = float(format(obj_val, '.2f'))
        record = self.response_record[self.current_gene][self.model]['response'].get(obj_val_interval) # closest alpha in the obj_interval
        interval_abs_diff = abs(obj_val - obj_val_interval) 

        if record and (record['interval_abs_diff'] >= interval_abs_diff): # erase record, update to current alpha
            record = None  
        if record and (obj_val_interval > .99) and (record['search_alpha']>self.search_alpha): # IC0 update most large alpha
            record = None

        if self.opt_df is not None:
            previous_min = min([ele['target_abs_diff'] for ele in 
                    self.response_record[self.current_gene][self.model]['response'].values()])
            # print(interval_abs_diff, previous_min)
            if interval_abs_diff <= previous_min:
                is_lowest_abs_diff = True

        if not record:
            self.response_record[self.current_gene][self.model]['response'][obj_val_interval] = {
                'search_alpha' : self.search_alpha,
                'precise_obj_val' : obj_val, 
                'interval_abs_diff' : interval_abs_diff,
                'target_abs_diff' : abs(obj_val - self.target_obj_val)}
        return is_lowest_abs_diff

    def evaluate_alpha_from_biomass(self):
        # print(self.alpha_ub)
        _, obj_val, summary_df = Sij_biomass(self.model, self.search_alpha, self.current_gene)
        
        if(self.opt_df is None and summary_df['Net_Flux'][0]<1e-8): # reinitialize if the first alpha to search gives zero flux
            self.search_alpha = 1+2e-3
            _, obj_val, summary_df = Sij_biomass(self.model, self.search_alpha, self.current_gene)

        # self.is_new_ub = (summary_df['Net_Flux'][0]<1e-8 or summary_df['Net_Flux'][0]>100 or 
        #                   obj_val<self.target_obj_val*.85 or obj_val>1) # overinhibition, search alpha too high
        self.trace_obj_val.append(obj_val)
        obj_val = obj_val/ self.norm_obj_val # to normalized 
        self.is_new_ub = (summary_df['Net_Flux'][0]<1e-8 or summary_df['Net_Flux'][0]>100 or 
                          obj_val<self.target_obj_val*self.acceptance_threshold_lower or obj_val>1.2) # overinhibition, search alpha too high 
        
        # print('obj', obj_val, self.target_obj_val)
        
        self.found_alpha = self.is_found(self.search_alpha, self.alpha_lb, self.alpha_ub, self.precision)
        
        # update optimal df
        net_flux_req = (summary_df['Net_Flux'][0]>0 or obj_val>self.norm_obj_val*.1)
        is_lowest_abs_diff = self.fill_response_record(obj_val)

        # obj_req = ((obj_val > self.target_obj_val*self.acceptance_threshold_lower) and 
        #            (obj_val < self.target_obj_val*self.acceptance_threshold_upper)) 
        
        obj_req = is_lowest_abs_diff
        if (net_flux_req and obj_req) or (self.opt_df is None): # store only qualified alpha
            self.opt_df = summary_df
        # print(summary_df['div_opt_obj_val'][0]/ self.norm_obj_val)
        
        # print(self.search_alpha, self.found_alpha, self.is_new_ub, self.alpha_lb, self.alpha_ub)

    def modify_opt_df(self):
        opt_df = self.opt_df
        opt_df['is_growth_switch'] = self.classify_growth_switch()
        print(self.current_gene, opt_df['div_opt_alpha'][0], opt_df[f'div_opt_obj_val'][0])
    #     opt_df['div_opt_alpha'] = np.floor(opt_df['div_opt_alpha'][0]*10**(precision))/10**(precision)
        opt_df.insert(loc=2, column='Percent_target_obj_val', value=opt_df['div_opt_obj_val']/self.norm_obj_val)
    #     opt_df.insert(loc=3, column='Essentiality', value=)
        opt_df.columns = [f'{self.model.id}_'+element for element in opt_df.columns]
        self.opt_df = opt_df
        return (opt_df[f'{self.model.id}_div_opt_alpha'],opt_df[f'{self.model.id}_div_opt_obj_val'],opt_df) # alpha_feasible, and upper bound of alpha

# %%
# with open("./Drug_gene_Similarity/Data/potential_genes.json", "r") as fp:
#     potential_genes = json.load(fp)
potential_genes = [
    'glyA', 'gltA', 'tktA', 'dapF', 'dapB', 'acnB', 'pgk', 'talB', 'mrcA', 'pyrE', 'dapD', 'pfkA', 'gdhA', 'folA', 'mrdA', 'thrB',
    'dapA', 'serC', 'argD', 'thrC', 'aceF', 'pykF', 'dadX', 'folC', 'pyrD', 'trpA', 'serB', 'fbp', 'eno', 'pgi', 'pheA',
    'gcvH', 'gnd', 'murA', 'aroA', 'guaB', 'glnA', 'yrbG', 'folP', 'purU', 'serA', 'gltD', 'purT', 'ackA', 'purN', 'rffG',
    'gapA'
]

##--------Perform Monoculture Secrch---------##
def Sij_biomass(model, alphas = 1, genes = 'folA', slim_opt = False, pfba = False): 
    with model:
        rct_ids = alter_Sij(model, alphas, genes) # same

        # need run sim and cal gr 
        if (slim_opt == False):
            sol_sol = model.optimize()   
            if pfba:
                return sol_sol, cobra.flux_analysis.pfba(model)

            obj_val = sol_sol.objective_value     
            summary_df = get_summary_df(model, alphas, obj_val, rct_ids, sol_sol.fluxes)  
            return(alphas, obj_val, summary_df)            
        else: 
#             return(pd.DataFrame([{f'div_opt_obj_val': model.slim_optimize()}]))
            return(model.slim_optimize()) 

def get_single_div_obj_df(model, target_obj_val, first_n_gene=None, alpha_df = None, search_alpha = None,precision=6,
                   potential_genes=potential_genes, acceptance_threshold_upper=1.1, acceptance_threshold_lower=.95): # media not required for GR

#     # generate full biomass dataframe for a single species
#     # given alpha or automatic search for optimal with the given objective value
#     if not first_n_gene: first_n_gene = len(potential_genes)
    obj_div_df = pd.DataFrame() # obj_val: biomass/ growth rate

    if first_n_gene is None: 
        query_gene_list = potential_genes
    elif isinstance(first_n_gene, (int, float)):
        query_gene_list = list(potential_genes)[:first_n_gene]
    else:
        query_gene_list = list(potential_genes)[first_n_gene[0]:first_n_gene[1]]

    for i, current_gene in enumerate(query_gene_list): # iter gene 
#         for i, current_gene in enumerate(obj_flux_df.index[1:][:1]): # iter gene 
        print(current_gene,'in cal')
        with model:
            
            if alpha_df is None: 
                if not search_alpha: search_alpha = 1.02 
                # AlphaFinder class for every model and SG
                AF = MonocultureAlphaFinder(model=model,
                            search_alpha = search_alpha,
                            current_gene = current_gene, 
                            target_obj_val = target_obj_val, 
                            exp_leap=3,
                            precision=precision,
                            acceptance_threshold_upper = acceptance_threshold_upper,
                            acceptance_threshold_lower = acceptance_threshold_lower)
                alpha, obj_value,temp_df = AF.find_feasible_alpha()
            else:
                alpha = alpha_df.loc[current_gene,f'{model.id}_div_opt_alpha'] # get alpha for corresponding gene from alpha_df
                alpha,obj_value,temp_df = Sij_biomass(model, alpha,  
                                                genes = current_gene, model_id = model.id)
            temp_df['Gene_inhibition'] = current_gene    
    #         print(alpha, current_gene, )
            obj_div_df = pd.concat([obj_div_df, temp_df.set_index('Gene_inhibition')],axis=0)
    return obj_div_df

def get_div_obj_df(model_list, target_obj_val, potential_genes=potential_genes, precision=4):
    alpha_obj_df_list = iter_species(model_list, get_single_div_obj_df,
                            target_obj_val=target_obj_val, potential_genes=potential_genes, precision=precision)
    return pd.concat(alpha_obj_df_list, axis=1)

# %% [markdown]
# 
# 
# # COC

# %%
@dataclass(kw_only=True)
class CocultureAlphaFinder(AlphaFinderConfig): # scale normal & knockout to 1-0
    model: List 
    search_alpha : float
    current_gene : str
    
    target_obj_val : float = 0.5 
    alpha_table : pd.DataFrame = None
    precision : int = 1 # precision of alpha
#     alpha_table : pd.DataFrame = alpha_table # init AFConfig
    trace_biomass : List = field(default_factory=list)
    trace_obj_val : List = field(default_factory=list) # TODO: combine all trace_XX & convert to dict, follow format of response_record in MonocultureAlphaFinder
    trace_alpha : List = field(default_factory=list)
    gr_ko : float = None
    n_dir : str = ''
    opt_alpha_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    nxt_alpha_table: pd.DataFrame = field(default_factory=pd.DataFrame)
    ko_intercepted : None
    acceptance_threshold_lower : float = .95
    acceptance_threshold_upper : float = 1.05
    exp_leap = None
    gr_Normal : float = None
    carbon_source_val : float = 5e-3
    add_nutrient_val : float = 10 # DO NOT set Met limiting in batch culture 
    is_growth_switch = False
    initial_pop : float = 1e-8
    p : None = None
    obj_style: str 
    
    # add_nutrient_val : float = field(default_factory=[.08])
    def __post_init__(self):
        super().__post_init__()
        print(self.exp_leap)
        # print(type(self.model[1]))

        self.is_monoculture = False
        self.eval_alpha_fun = self.evaluate_alpha_from_gr 
        self.get_alpha_use()
        self.E0 = [ele for ele in convert_arg_to_list(self.model) if ele.id == 'E0'][:1]
        # self.gr_Normal = 0.00064644123538125 # ?check effect of lcts conc changed in Div_COMETS
        self.alpha_pair : pd.DataFrame = None
        if not self.exp_leap:
            self.exp_leap = 1.3
        self.E0 = self.model[0]
        self.S0 = self.model[1]
        self.out_fun = self.out_opt_alpha_table
        self.calculate_gr_Normal()
#         _ = self.calculate_gr_ko() # moved to eval_alpha
    
    def calculate_gr_Normal(self):
        if not self.gr_Normal:
            if self.current_gene == 'Normal':
                self.gr_Normal = self.calculate_gr_ko()
            else: self.gr_Normal = 0.28672261
        return float(self.gr_Normal)

    def calculate_gr_ko(self):
        print('calculating ko gr')
        self.ko_intercepted = False
        
        if self.current_gene in ['folP', 'folA', 'dadX', 'pgk', 'acnB']:
            self.gr_ko = 0
            return self.gr_ko
        
        # alpha_in_table = self.alpha_table.loc[self.current_gene, self.E0.id] # E0 col
        full_df, *_ =  get_BM_df(self.current_gene,n_dir=self.n_dir,alpha_table=self.alpha_table,
                                E0=self.E0, S0=self.S0, mono=False, p=self.p, return_sim=True, ko=True, 
                                carbon_source_val=self.carbon_source_val, add_nutrient_val=self.add_nutrient_val,
                                initial_pop=self.initial_pop, obj_style=self.obj_style)
        
        coculture_biomass_df = full_df.iloc[:,[0]]
        if self.current_gene == 'Normal': 
            coculture_biomass_df.to_csv(f'./Data/Normal_biomass_{str(self.carbon_source_val)}.csv')
        
        coculture_biomass_df.columns = rename_columns(coculture_biomass_df)
        
        # desired_cycle = get_desired_cycle(coculture_biomass_df, scale_diff=0.05)
        # # return desired_cycle, full_df, coculture_biomass_df
        # growth_phase = desired_cycle.growth_phase[0]
        # gr_ko = (coculture_biomass_df.loc[growth_phase[1]] - coculture_biomass_df.loc[growth_phase[0]])/(growth_phase[1]-growth_phase[0])
        # gr_ko = float(gr_ko)
        
        gr_ko = self.cal_fitted_gr(coculture_biomass_df)
            
        self.cocdf = coculture_biomass_df 
        # not setting alpha, not upper bound
        self.gr_ko = float(gr_ko)
        
        if self.current_gene == 'Normal':
            print('set_normal')
            self.gr_Normal = gr_ko
        elif self.gr_ko > self.gr_Normal*.9:
            print('intercept')
            self.ko_intercepted = True
            self.opt_alpha_table = self.alpha_table.loc[[self.current_gene]]
        else: print('NNNOT intercept')
        
        print(self.opt_alpha_table)
        
        print(f'self.gr_ko > self.gr_Normal: {self.gr_ko},{self.gr_Normal*.85}')
        print(' gr_ko: ', gr_ko)
        self.trace_biomass.append(full_df.add_suffix(f'_ko'))
        return gr_ko
    
    def get_alpha_use(self):
        if self.current_gene == 'Normal':
            return 0
        Ealpha = self.alpha_table.loc[self.current_gene].iloc[0]
        Salpha = self.alpha_table.loc[self.current_gene].iloc[1]
        self.alpha_use = Ealpha if Ealpha<1e5 else Salpha
        if self.alpha_use > 1e4:
            self.alpha_use = 1.04
            self.exp_leap = 3
        
        if not self.search_alpha and (self.current_gene != 'Normal'):
                self.search_alpha = (self.alpha_use-1)/2+1 if self.alpha_use > 2 else self.alpha_use
        return None
    
    def generate_next_fixed_ratio_table(self):     
        multiplier = (self.search_alpha-1)/(self.alpha_use-1)
        nxt_alpha_table =  (self.alpha_table.loc[[self.current_gene]]-1)*multiplier + 1 # ratio offset by 1, minimum is alpha=1
         # trace alpha table change and corresponding gr
        # self.trace_alpha_table.append(nxt_alpha_table.to_dict())
        self.nxt_alpha_table = nxt_alpha_table
        return nxt_alpha_table
     
    def cal_fitted_gr(self, coculture_biomass_df):
        # def logistic(x, saturation, growth_rate, inflection_point, initial_pop):
        #     return (saturation / (1 + np.exp(-growth_rate * (x - inflection_point)))) + initial_pop
        
        def np_loglinear(s: pd.Series):
            np_gr, est_initial_pop = np.polyfit(s.index, np.log(s), 1) # [B,A] -- log(y) = A+Bx
            return np_gr
            
        def truncate_max_biomass(fit_col):
            cycles_to_fit = fit_col<fit_col.iloc[-1]*.9
            return pd.Series(fit_col.loc[cycles_to_fit])
        
        # x = coculture_biomass_df.index
        y = truncate_max_biomass(coculture_biomass_df.iloc[:, 0]) # E0 column
        
        print(self.current_gene, y.tail(1), y.head(1), 'to convert to float')
        if len(y) < 5 or (float(y.tail(1)) < float(y.head(1))+5e-8): #intercept minimal growth
            growth_rate = 0
            return growth_rate
            
        # desired_cycle = get_desired_cycle(coculture_biomass_df, scale_diff=0.05) 
        
        # growth_phase = desired_cycle.loc[self.current_gene, 'growth_phase']
        # bool_growing = desired_cycle.loc[self.current_gene, 'bool_growing']
        # c_max_gr = desired_cycle.loc[self.current_gene, 'c_max_gr']
        # saturation = float(coculture_biomass_df.loc[growth_phase[1]])
        # init_gr = .23
        # p0 = [saturation, init_gr, c_max_gr, self.initial_pop]

        # print('growing: ',bool_growing)
        # if (bool_growing == False) or (bool_growing <1):
        #     popt, pcov = curve_fit(logistic, x, y, p0=p0)
        #     print('logi')
            
        #     if self.current_gene =='Normal':
        #         print(popt)
        #     growth_rate = popt[1]
        # else:
        #     print('lolin')
        #     popt, pcov = np.polyfit(x, np.log(y), 1) # [B,A] -- log(y) = A+Bx
        #     growth_rate = popt
        # if self.current_gene == 'Normal':
        #     print('popt', popt)
        growth_rate = np_loglinear(y)
         
        return growth_rate
        
    def cal_std_gr(self, coculture_biomass_df, gr_Normal):
#         desired_biomass_df = pd.DataFrame(get_cycle_max_gr(coculture_biomass_df), index=['max_gr', 'start', 'end', 'bool_growing'])
        gr = self.cal_fitted_gr(coculture_biomass_df)
        standardized_gr = (float(gr) - self.gr_ko)/(gr_Normal-self.gr_ko) # offset scale self.gr_ko then scale
        if standardized_gr > 1:
            self.exp_leap = 3 # for coculture gr inherently > 1, require large leap
        
        # print('standardized gr: ', standardized_gr)
        self.trace_obj_val.append(standardized_gr)
        print('---GR, stdGR__:',gr, standardized_gr)
        return standardized_gr
    
    def evaluate_alpha_from_gr(self, alpha_overwrite=None): # first, search alppha = alpha_table alpha
        if self.ko_intercepted is None:
            _ = self.calculate_gr_ko() # also catch gr_ko ~ gr_Normal- Gene with nonessential reactions
        
        if alpha_overwrite:
            self.search_alpha = alpha_overwrite
        nxt_alpha_table = self.generate_next_fixed_ratio_table()
        # self.trace_alpha_table.append(nxt_alpha_table.to_dict())

        print(self.p.all_params['maxCycles'])
        full_df, out_dict, co_sim = get_BM_df(self.current_gene,n_dir=self.n_dir,alpha_table=self.nxt_alpha_table,
                                E0=self.E0, S0=self.S0, mono=False, p=self.p, return_sim=True, ko=False, 
                                carbon_source_val=self.carbon_source_val, add_nutrient_val=self.add_nutrient_val,
                                initial_pop=self.initial_pop, obj_style=self.obj_style)
        # get_BM_df(self.current_gene,n_dir=self.n_dir,alpha_table=nxt_alpha_table,mono=False, return_sim=True)
        self.co_sim = co_sim
        
        self.trace_biomass.append(full_df.add_suffix(f'_{self.i_iter}'))
        
        coculture_biomass_df = full_df.iloc[:,[0]]
        target_gr = self.target_obj_val
        std_gr = self.cal_std_gr(coculture_biomass_df, self.gr_Normal) # ?check this
#         std_gr = self.cal_std_gr(coculture_biomass_df, 1.3)

        self.is_new_ub = (std_gr < target_gr)# sharp bound, although lower estimate, no search for higher alpha if just meet
        self.converge_alpha = self.is_found(self.search_alpha, self.alpha_lb, self.alpha_ub, self.precision)
        
        # self.found_alpha = self.found_alpha or 
        #     (self.std_gr > self.target_obj_val*acceptance_threshold_upper)

        self.obj_found = ((std_gr > self.target_obj_val*self.acceptance_threshold_lower) and 
            (std_gr < self.target_obj_val*self.acceptance_threshold_upper)) 
        
        bool_gr_coculture_exceed_limit = (std_gr > 1)
        self.found_alpha = self.converge_alpha or self.obj_found or self.classify_growth_switch()

        if self.is_new_ub:
            self.nxt_alpha_table = nxt_alpha_table
        print('searc, found, new_un, lu')
        print(self.search_alpha, self.found_alpha, self.is_new_ub, self.alpha_lb, self.alpha_ub)
        print('convg, obj_f, exceed, gr_sw')
        print(self.converge_alpha, self.obj_found, bool_gr_coculture_exceed_limit, self.classify_growth_switch())
        # update optimal df
        
        obj_req = ((std_gr > self.target_obj_val*0.95) and (std_gr < self.target_obj_val*1.1)  # harsh lb
                    or self.is_growth_switch)
        
        temp_df = pd.DataFrame.from_dict(
            {self.current_gene :
                {'std_growth_rate' : std_gr,
                'lb_ub' : (self.alpha_lb, self.alpha_ub),
                'converge_alpha' : self.converge_alpha,
                'obj_found': self.obj_found,
                'bool_gr_coculture_exceed_limit': bool_gr_coculture_exceed_limit,
                'classify_growth_switch': self.is_growth_switch,
                'Non_essential': self.ko_intercepted
                }}, 
                orient='index')
        nxt_alpha_table = pd.concat([nxt_alpha_table, temp_df], axis=1)
        nxt_alpha_table.index.name='Gene_inhibition'
        self.trace_alpha.append(nxt_alpha_table.to_dict())
        
        if obj_req or (self.opt_df is None) or self.converge_alpha: # store only qualified alpha
            # if self.opt_df is not None: # do not append the first 
            #     pass
            #     temp_df = pd.DataFrame.from_dict(
            #         {self.current_gene :
            #             {'converge_alpha' : self.converge_alpha,
            #             'obj_found': self.obj_found,
            #             'bool_gr_coculture_exceed_limit': bool_gr_coculture_exceed_limit,
            #             'classify_growth_switch': self.is_growth_switch}}, 
            #             orient='index')
            # else: temp_df = pd.DataFrame()
            if self.opt_df is not None:
                nxt_alpha_table['obj_req_satisfied'] = True
            self.opt_df = [full_df, nxt_alpha_table]
            # self.opt_alpha_table = pd.concat([nxt_alpha_table, temp_df], axis=1)
            self.opt_alpha_table = nxt_alpha_table
        return full_df 
    
    def out_opt_alpha_table(self):
        result = {}
        for i in range(len(self.trace_obj_val)):
            result[f'iter_{i+1}'] = {'obj_val': self.trace_obj_val[i], 'alpha_table': self.trace_alpha[i]}
        
        biomass_df = pd.concat(self.trace_biomass, axis=1) if self.trace_biomass else pd.DataFrame()
        return biomass_df, self.opt_alpha_table, {self.current_gene : result}
#         return 'Temp Done' if (n_iter<=10 or n_iter ==2) else 'end at: '+ str(n_iter)

##--------Perform Coculture Secrch---------##

def coculture_search_job(**kwargs):     
    AF = CocultureAlphaFinder(**kwargs)
    # return pd.DataFrame([]), {0:{0:0}}
    AF.calculate_gr_ko()
    return AF.find_feasible_alpha()

def run_coculture_search_mp(potential_genes, filename, n_processor, **kwargs):
    def save_trace_biomass(trace_biomass):
        trace_biomass = pd.concat(trace_biomass,axis=1) # alpha_table
        trace_biomass.columns = rename_columns(trace_biomass)
        trace_biomass.to_csv(f'./Data/biomass_{filename}.csv')
     
    result_list = list()
     
    with concurrent.futures.ProcessPoolExecutor(n_processor) as executor:
        for i, current_gene in enumerate(potential_genes):
            future = executor.submit(coculture_search_job, current_gene=current_gene, n_dir=i, **kwargs)
            result_list.append(future)
    
    trace_biomass, opt_alpha_list, result_dict_list = zip(*[r.result() for r in result_list])
    opt_alpha = pd.concat(opt_alpha_list) # alpha_table
        
    opt_alpha.columns = rename_columns(opt_alpha)
    opt_alpha.to_csv(f'./Data/alpha_table_{filename}.csv')
    print(opt_alpha, result_dict_list)
    
    trace_dict = {k: v for d in result_dict_list for k, v in d.items()}
    with open(f"./Data/trace_record_{filename}.json", "w") as outfile: 
        json.dump(trace_dict, outfile) 
        
    
    save_trace_biomass(trace_biomass)
    print('finished alpha search')
    return opt_alpha_list, trace_dict