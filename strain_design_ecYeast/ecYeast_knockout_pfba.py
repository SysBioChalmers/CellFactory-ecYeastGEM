###################################
#  application for ecYeast
###################################
from cobra.io import load_matlab_model, read_sbml_model
import os
import sys
import matplotlib.pyplot as plt # For python 3.6: sudo apt-get install python3.6-tk   TO INSTALL Tkinter
import pandas as pd
from cobra.flux_analysis import pfba
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

from cobra.manipulation import remove_genes
from cobra.manipulation import *
from cobra.flux_analysis.moma import moma


sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")
sys.path.append(r"/Users/luho/PycharmProjects/model/strain_design/code")
os.chdir('/Users/luho/PycharmProjects/model/strain_design/code')



# import self function
from mainFunction import *
from in_silico_strain_design_scripts import *


# load the ecYeast model
ecYeast = load_matlab_model("../data/ec3HP.mat")


"""knock-out"""
# pFBA method is used
# run the gene deletion with pFBA method, minimization of proteins pool
# in this simulation, 3HP production is as the objective function
# qs is unconstrait
def getModelWithRemoveGene (model0, gene_remove0):
    model1 = model0.copy()
    remove_genes(model1, gene_remove0, remove_reactions=True) # this script can't be used
    return model1

def solveEcModel(model, obj_id, target_id, glucose_uptake ='r_1714_REV'):
    model.objective = obj_id  # we put it in the mutant model
    solution = model.optimize()  # we put it in the mutant model
    # mutModel.summary() # we put it in the mutant model
    solution_value0 = solution.objective_value  # we put it in the mutant model
    model.reactions.get_by_id(obj_id).bounds = (0.999 * solution_value0, solution_value0)  # we put it in the mutant model
    model.objective = {model.reactions.prot_pool_exchange: -1}
    solution2 = model.optimize()
    yield_p = solution2.fluxes[target_id]/solution2.fluxes[glucose_uptake]
    return  yield_p

# pipeline

# prepare the geneID
geneID = []
for x in ecYeast.genes:
    print(x.id)
    geneID.append(x.id)

# prepare ecModel for calculation
model2 = ecYeast.copy()
model2 = ecYeastMinimalMedia(model2)
model2.reactions.get_by_id("r_1714_REV").bounds = (0,1000)  # D-glucose exchange (reversible)

# using pfba
yield_p_ref0 =[]
yield_p_mutant0 =[]
for i, x in enumerate(geneID):
    print(i)
    gene_remove = []
    gene_remove.append(x)
    mutModel = getModelWithRemoveGene(model0=model2, gene_remove0=gene_remove)

    try:
        yield_p_ref = solveEcModel(model=model2, obj_id = 'newRxn3', target_id = 'newRxn3')
        yield_p_mutant = solveEcModel(model=mutModel, obj_id='newRxn3', target_id='newRxn3')
        print(yield_p_ref)
        print(yield_p_mutant)
        yield_p_ref0.append(yield_p_ref)
        yield_p_mutant0.append(yield_p_mutant)
    except:
        yield_p_ref = None
        yield_p_mutant = None
        yield_p_ref0.append(yield_p_ref)
        yield_p_mutant0.append(yield_p_mutant)
pfba_result = pd.DataFrame({'gene':geneID, 'yield_p_ref0': yield_p_ref0, 'yield_p_mutant0': yield_p_mutant0 })
saveExcel(pfba_result, "../result/pfba_result.xlsx")