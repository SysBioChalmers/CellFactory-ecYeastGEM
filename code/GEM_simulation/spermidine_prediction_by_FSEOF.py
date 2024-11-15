import cobra
from cobra import Model, Reaction, Metabolite
import os
from os.path import join
import numpy as np
import pandas as pd
import copy
import re
#model = cobra.io.load_matlab_model(r"C:\Users\Administrator\Desktop\Augest\y8.mat")
# model = cobra.io.load_matlab_model("examples/siwei_product_FSEOF/y8.mat")
model=cobra.io.read_sbml_model(r'ModelFiles/xml/yeastGEM.xml')

Ysx = 0.122*180/1000

def simulateGrowth(model, alpha):
    tmpmodel = model.copy()
    tmpmodel.reactions.get_by_id('r_1714').bounds = -5,0
    # tmpmodel.reactions.get_by_id('r_1761').bounds = (0,1000)
    # model.reactions.get_by_id('r_1761_REV').bounds = (0,0)

    # max growth
    tmpmodel.objective = 'r_2111'
    sol = tmpmodel.optimize()

    # fix growth and max product
    tmpgrow = sol.objective_value * 0.999*alpha
    tmpmodel.reactions.get_by_id('r_2111').bounds = (tmpgrow, tmpgrow)
    # product_id='r_1589'
    product_id = 'r_2051'
    tmpmodel.objective = product_id
    product_max = tmpmodel.optimize()
    # print(product_max.objective_value)

    # fix product and perform pFBA
    tmpmodel.reactions.get_by_id(product_id).bounds = (product_max.objective_value * 0.999,
                                                  product_max.objective_value * 0.999)
    sol_pfba = cobra.flux_analysis.pfba(tmpmodel)
    return sol_pfba


def div_(start_, end_, num):
    l = []
    interval = (end_ - start_)/num
    co = 0
    while co < num + 1:
        l.append(start_ + co * interval)
        co += 1
    return l


def compare_substrate(model, Ysx):
    flux_WT = simulateGrowth(model, 1)
    wtproduct = Ysx / flux_WT.fluxes.loc['r_2111']
    alpha = div_(wtproduct / 2, wtproduct * 2, 15)
    v_matrix = []
    k_matrix = []
    f_gene = []
    for a in alpha:
        print(a)
        tmpflux = simulateGrowth(model, a)
        try:
            v_matrix.append(tmpflux.fluxes)
            k_matrix.append(np.asarray(tmpflux.fluxes)/
                            np.asarray(flux_WT.fluxes))
        except AttributeError:
            v_matrix.append(tmpflux)
            k_matrix.append(tmpflux)
        # print(a)
        # print(alpha.index(a))


    return v_matrix, k_matrix, f_gene, alpha

v_matrix, k_matrix, f_gene, alpha = compare_substrate(model, Ysx)


all_flux = pd.DataFrame(v_matrix)   # 生成dataframe
all_flux = all_flux.T


all_k = pd.DataFrame(k_matrix)      # 生成dataframe
all_k = all_k.T
all_k.index = all_flux.index

# delete rxn without genes
f1 = list()
for h in all_flux.index:
    if bool(model.reactions.get_by_id(h).gene_reaction_rule == ''):
        f1.append(h)
    else:
        pass
all_flux.drop(index = f1, inplace = True)
all_k.drop(index = f1, inplace = True)

# delete all nan
all_flux.dropna(axis=0, how='all', inplace=True)
all_k.dropna(axis=0, how='all', inplace=True)

# nan --> 1
all_flux.fillna(1, inplace=True)
all_k.fillna(1, inplace=True)

# inf --> 1000
all_k.replace([np.inf, -np.inf], 1000, inplace=True)

# delete inconsistent
f5 = []
for i in all_k.index:
    big = any(all_k.loc[i,:] > 1)
    small = any(all_k.loc[i,:] < 1)
    if big == small == True:
        f5.append(i)
all_k.drop(index = f5, inplace = True)

# order k(descend)
orderd_k = copy.deepcopy(all_k)
orderd_k['meank'] = orderd_k.mean(axis =1)
orderd_k = orderd_k.sort_values(by=["meank"], ascending=False)


cons_g = {}
for ge in model.genes:
    # 取基因的k值
    ge_rxn = [r.id for r in ge.reactions if r.id in orderd_k.index]
    ge_k = [orderd_k.loc[r.id, 'meank'] for r in ge.reactions if r.id in orderd_k.index]
    # 看一致性
    if ge_k != []:
        gk = np.asarray(ge_k)
        one = np.array(len(gk))
        if sum(gk >= 1) == one or sum(gk <= 1) == one:
            cons_g[ge.id] = np.mean(ge_k)


alpha_m = np.mean(alpha)
cons_g_f = {}
for k, v in cons_g.items():
    #if (v > 1e-3):
        if (v < alpha_m - 1e-3) or (v > 1 + 1e-3):
            cons_g_f[k] = v


cons_g_f = pd.Series(cons_g_f)
cons_g_f.sort_values(ascending=False, inplace=True)

# classification
# k_score > 1, overexpression
# 0.05 < k_score <= 0.5 down-regulation
# k_score <= 0.05 deletion
cons_g_f0 = pd.DataFrame(cons_g_f)
cons_g_f0.columns = ['k_scores']
actions = []
for i, x in cons_g_f0.iterrows():
    if x['k_scores'] > 1.1:
        actions.append('OE')
    elif x['k_scores'] <= 0.05:
        actions.append('KO')
    elif x['k_scores'] <= 0.5 and x['k_scores'] > 0.05:
        actions.append('KD')
    else:
        actions.append('Not_sure')

cons_g_f0['actions'] = actions
cons_g_f0['gene'] = list(cons_g_f0.index)

#cons_g_f.to_excel(r"C:\Users\Administrator\Desktop\product production enhance targets.xlsx")
cons_g_f0.to_excel(r"code/GEM_simulation/output/fseof_result.xlsx")



