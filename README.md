# CellFactory-yeast-GEM: In silico strain design of _Saccharomyces cerevisiae_ in large scale



* Brief Model Description:

This repository contains the data and script for _in silico_ strain with protein constrainted genome-scale metabolic model of _Saccharomyces cerevisiae_. 

* Model KeyWords:

**GEM Category:** Species; **Utilisation:** predictive simulation, multi-omics integrative analysis, _in silico_ strain design; **Field:** cell factory _in silico_ design; **Model Source:** [Yeast 8.3](https://github.com/SysBioChalmers/yeast-GEM); **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** General Metabolism; **Condition:** aerobic, glucose-limited, defined media, maximization of growth.

* Last update: 2018-12-01

* Main Model Descriptors:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 8.3](https://github.com/SysBioChalmers/yeast-GEM)|3962|2691|1150|



## Installation

### Required Software - User:

* Matlab user:
  * A functional Matlab installation (MATLAB 7.3 or higher).
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
* Python user:
  * Python 2.7, 3.4, 3.5 or 3.6
  * [cobrapy](https://github.com/opencobra/cobrapy)

### Required Software - Contributor:

* Both of the previous Matlab requirements.
* The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN).
* A [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.

### Dependencies - Recommended Software:
* For Matlab, the [libSBML MATLAB API](https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface/) (version 5.17.0 is recommended).
* [Gurobi Optimizer](http://www.gurobi.com/registration/download-reg) for any simulations.

### Installation Instructions
* For users: Clone it from [`master`](https://github.com/SysBioChalmers/yeast-GEM) in the Github repo, or just download [the latest release](https://github.com/SysBioChalmers/yeast-GEM/releases).
* For contributors: Fork it to your Github account, and create a new branch from [`devel`](https://github.com/SysBioChalmers/yeast-GEM/tree/devel).

## Usage

Make sure to load/save the model with the corresponding wrapper functions!
* In Matlab:
  * Loading: `complementaryScripts/loadYeastModel.m`
  * Saving: `complementaryScripts/saveYeastModel.m`
* In Python:
  * Loading: `complementaryScripts/loadYeastModel.py`
  * Saving: currently unavailable

## Model Files

The model is available in `.xml`, `.txt`, `.yml`, `.mat` and `.xlsx` (the last 2 extensions only in `master`). Additionally, the following 2 files are available:
* `dependencies.txt`: Tracks versions of toolboxes & SBML used for saving the model.
* `boundaryMets.txt`: Contains a list of all boundary metabolites in model, listing the id and name. 

### Complementary Scripts

* `missingFields`: Folder with functions for adding missing fields to the model.
* `modelCuration`: Folder with curation functions.
* `otherChanges`: Folder with other types of changes.
* `increaseVersion.m`: Updates the version of the model in `version.txt` and as metaid in the `.xml` file. Saves the model as `.mat` and as `.xlsx`
* `loadYeastModel.m`: Loads the yeast model from the `.xml` file for Matlab.
* `loadYeastModel.py`: Loads the yeast model from the `.xml` file for Python.
* `saveYeastModel.m`: Saves yeast model as a `.xml`, `.yml` and `.txt` file, and updates `boundaryMets.txt` and `dependencies.txt`.

### Complementary Data

* `databases`: Yeast data from different databases (KEGG, SGD, swissprot, etc).
* `modelCuration`: Data files used for performing curations to the model. Mostly lists of new rxns, mets or genes added (or fixed) in the model.
* `physiology`: Data on yeast growth under different conditions, biomass composition, gene essentiality experiments, etc.

## Contributing

Contributions are always welcome! Please read the [contributions guideline](https://github.com/SysBioChalmers/yeast-GEM/blob/master/.github/CONTRIBUTING.md) to get started.
  
## Contributors

* [Eduard J. Kerkhoven](https://www.sysbio.se/people/eduard-kerkhoven/) ([@edkerk](https://github.com/edkerk)), Chalmers University of Technology, Sweden
* [Feiran Li](https://www.sysbio.se/people/feiran-li) 
([@feiranl](https://github.com/feiranl)), Chalmers University of Technology, Sweden
* [Hongzhong Lu](https://www.sysbio.se/people/hongzhong-lu) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Sweden
* [Yao Lu](hhttps://www.sysbio.se/people/yao-lu), Chalmers University of Technology, Sweden
