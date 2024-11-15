# ecYeast cell-factory

# Instructions

* All scripts are refered relative to their location inside of the subfolder `cellFactory-yeast-GEM/code`.

1. Run MATLAB script `eval_fluxDist.m` to evaluate biomass yield, product yield and product rate for each of the +100 chemicals using ecYeastGEM (wild-type). This script also produces a matrix of stoichiometric distances between optimal distributions for maximum production of each chemical.

2. Run MATLAB script `run_predictions.m` to run the `ecFactory` method on all the models available in this repository. This will generate lists of gene target predictions for enhanced production of +100 chemicals in *S. cerevisiae*.

3. Run MATLAB scripts `summarize_results.m` and `getTargetsMatrix.m` to obtain a summary of the number of predicted targets for OE, KD and KO for each product, and a representation of this in a discrete matrix form, respectively.

4. Run R script `analyse_results.R`. This script will generate:
 	- Boxplots, showing distribution of number of targets, per type, for all products and levels of the multi-step robust ecFSEOF (1. FSEOF, 2. mechanisctically validated individual targets, 3. optimal combination of targets).

5. Run R script `geneCentricAnalysis.R` to perform analysis on the gene target profiles on a gene-centric manner. This will generate:
 - Bar plots showing the top10 most represented targets (for KO, KD and OE) across +100 chemicals
 - Bar plots showing the representation level of the top10 targets (for KO, KD and OE) within each the chemical families.
 - Histograms showing the target specificity (how many products are related to each gene target)
 
 6. Run R script `stoich_distance_analysis.R`, to perform a correlation analysis between pairwise distances between gene target profile vectors and optimal flux distribution vectors.

7. Run R script `geneTargets_by_prodCluster.R` for identification of common gene targets for all chemical within each cluster of products identified by t-SNE analysis.
