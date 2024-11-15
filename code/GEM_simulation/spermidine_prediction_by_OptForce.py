import pandas as pd
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
import os

def run_FVA(model, target_rxn, target_const_rate, flux_constraints={}):
    with model as m:
        m.objective = target_rxn
        m.objective_direction = 'max'
        for rxn_id in flux_constraints:
            m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]
        target_max = m.slim_optimize()

    fva_result = {}
    for rxn in model.reactions:
        with model as m:
            tmp = target_const_rate * target_max
            m.reactions.get_by_id(target_rxn).bounds = (tmp, tmp)

            for rxn_id in flux_constraints:
                m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]

            m.objective = rxn.id

            m.objective_direction = 'min'
            min_flux = m.slim_optimize()
            m.objective_direction = 'max'
            max_flux = m.slim_optimize()
            fva_result[rxn.id] = (min_flux, max_flux)
    return fva_result


def compare_FV_range(bio_result, target_result):
    down_targets = {}
    up_targets = {}
    for rxn_id in target_result:
        target_lb, target_ub = target_result[rxn_id]
        bio_lb, bio_ub = bio_result[rxn_id]
        if target_lb > bio_ub:
            up_targets[rxn_id] = {'bio':(bio_lb, bio_ub),
                                  'target':(target_lb, target_ub)}
        elif target_ub < bio_lb:
            down_targets[rxn_id] = {'bio':(bio_lb, bio_ub),
                                    'target':(target_lb, target_ub)}
    return down_targets, up_targets


def recordResults(target_rxn, down_targets, up_targets):
    if not os.path.exists('./output'): # If the directory does not exist, make it
        os.makedirs('./output')
        
    with open(f'./output/OptForce_{target_rxn}.txt', 'w') as fp:
        fp.write(f'Target reaction: {target_rxn}\n')
        fp.write('\nDown targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in down_targets:
            bio_lb, bio_ub = down_targets[rxn_id]['bio']
            target_lb, target_ub = down_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')
        fp.write('\nUp targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in up_targets:
            bio_lb, bio_ub = up_targets[rxn_id]['bio']
            target_lb, target_ub = up_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')

def rxns_to_genes(model,rxnIDlist):
    geneIDlist=list()
    for rxnID in rxnIDlist:
        for g in model.reactions.get_by_id(rxnID).genes:
            geneIDlist.append(g.id)

    geneIDlist=list(set(geneIDlist))
    return geneIDlist

        
if __name__ == '__main__':
    # os.chdir(r'code/GEM_simulation')
    model_dir = '../../ModelFiles/xml/yeastGEM.xml'
    biomass_rxn = 'r_2111'
    target_rxn = 'r_2051'

    biomass_const = 0.9
    target_const = 0.9

    model = read_sbml_model(model_dir)
    with model as m:
        bio_max = m.slim_optimize()

    bio_fva_result = run_FVA(model, biomass_rxn, biomass_const)

    ###
    flux_constraint = {biomass_rxn:(0.1*bio_max, 0.1*bio_max)}
    ###

    target_fva_result = run_FVA(model, target_rxn, target_const, flux_constraints=flux_constraint)

    down_targets, up_targets = compare_FV_range(bio_fva_result, target_fva_result)
    recordResults(target_rxn, down_targets, up_targets)

    # get gene targets
    down_targets_idlist=down_targets.keys()
    up_targets_idlist=up_targets.keys()
    down_gene_idlist=rxns_to_genes(model,down_targets_idlist)
    up_gene_idlist=rxns_to_genes(model,up_targets_idlist)

    # remove conflicting targets
    common_gene_idlist=list(set(down_gene_idlist)&set(up_gene_idlist))
    if len(common_gene_idlist)>0:
        down_gene_idlist=list(set(down_gene_idlist)-set(common_gene_idlist))
        up_gene_idlist=list(set(up_gene_idlist)-set(common_gene_idlist))

    df_gene_result=pd.DataFrame(down_gene_idlist+up_gene_idlist,columns=['geneID'])
    df_gene_result.loc[df_gene_result['geneID'].isin(down_gene_idlist),'target']='KD'
    df_gene_result.loc[df_gene_result['geneID'].isin(up_gene_idlist),'target']='OE'
    df_gene_result.to_csv(f'./output/OptForce_{target_rxn}_gene.csv',index=False)


