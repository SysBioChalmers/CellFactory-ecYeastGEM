% The theoretical yield analysis of Suc

% Set constraints (glyoxylate cycle) in yeastGEM model
cd ../ComplementaryScripts
model = loadYeastModel;
model=changeRxnBounds(model,'r_1714',-15.2,'l');
model=changeRxnBounds(model,'r_1889',-1,'l');
model=changeObjective(model,'r_2056');
% For wild type strain
model=changeRxnBounds(model,'r_2111',0.25,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh2?sdh1?idh1 strain
model=changeRxnBounds(model,'r_1021',0,'u');
model=changeRxnBounds(model,'r_0658',0,'u');
model=changeRxnBounds(model,'r_2111',0.22,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh2?sdh1?idh1?idp1 strain
model=changeRxnBounds(model,'r_2131',0,'u');
model=changeRxnBounds(model,'r_2111',0.2,'l');
FBAsolution=optimizeCbModel(model)

% Set constraints (glyoxylate cycle) in ecYeast model
cd ../ModelFiles/mat
load('ecYeastGEM_batch.mat');
model = ecModel_batch;
model=changeRxnBounds(model,'r_1714_REV',15.2,'u');
model=changeRxnBounds(model,'r_1889_REV',1,'u');
model=changeObjective(model,'r_2056');
% For wild type strain
model=changeRxnBounds(model,'r_2111',0.25,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh2?sdh1?idh1 strain
model = removeGenes(model,'YLL041C');
model = removeGenes(model,'YKL148C');
model = removeGenes(model,'YNL037C');
model=changeRxnBounds(model,'r_2111',0.22,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh2?sdh1?idh1?idp1 strain
model = removeGenes(model,'YDL066W');
model=changeRxnBounds(model,'r_2111',0.2,'l');
FBAsolution=optimizeCbModel(model)

% Set constraints (pathway link cell growth) in yeastGEM model
cd ../../ComplementaryScripts
model = loadYeastModel;
model=changeRxnBounds(model,'r_1714',-15.2,'l');
model=changeRxnBounds(model,'r_1810',-1,'l');
model=changeObjective(model,'r_2056');
% For wild type strain
model=changeRxnBounds(model,'r_2111',0.33,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh3?ser3?ser33 strain
model=changeRxnBounds(model,'r_1021',11.36,'u');% partially inactive
model=changeRxnBounds(model,'r_0891',0,'u');
model=changeRxnBounds(model,'r_2111',0.22,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh3?ser3?ser33ICL¡ü strain
model=changeRxnBounds(model,'r_1021',11.36,'u');% partially inactive
model=changeRxnBounds(model,'r_0891',0,'u');
model=changeRxnBounds(model,'r_0662',12.71,'l');% The flux of reaction contains ICL doubled
model=changeRxnBounds(model,'r_2111',0.12,'l');
FBAsolution=optimizeCbModel(model)

% Set constraints (pathway link cell growth) in ecYeast model
cd ../ModelFiles/mat
load('ecYeastGEM_batch.mat');
model = ecModel_batch;
model=changeRxnBounds(model,'r_1714_REV',15.2,'u');
model=changeRxnBounds(model,'r_1810_REV',1,'u');
model=changeObjective(model,'r_2056');
% For wild type strain
model=changeRxnBounds(model,'r_2111',0.33,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh3?ser3?ser33 strain
model = removeGenes(model,'YER081W');
model = removeGenes(model,'YIL074C');
model=changeRxnBounds(model,'r_2111',0.22,'l');
FBAsolution=optimizeCbModel(model)
% For ?sdh3?ser3?ser33ICL¡ü strain
model=changeRxnBounds(model,'r_0662',0.22,'l');% The flux of reaction contains ICL doubled
model=changeRxnBounds(model,'r_2111',0.12,'l');
FBAsolution=optimizeCbModel(model)

% Pathway construction of Suc in ecYeast model

% Add new reactions of reTCA pathway to the model
load('ecYeastGEM_batch.mat');
model = ecModel_batch;

% Add rxns of gene FumC
Kcat1=1150*3600;
MW1=50.489;
model = addReaction(model,'newRxn1','metaboliteList',{'s_0725','s_0803','prot_A0A024L3D9','s_0066'},'stoichCoeffList',[-1 -1 -1/Kcat1 1],'reversible',false);
Kcat2=11.2*3600;
model = addReaction(model,'newRxn2','metaboliteList',{'s_0066','prot_A0A024L3D9','s_0725','s_0803'},'stoichCoeffList',[-1 -1/Kcat2 1 1],'reversible',false);
model = addReaction(model,'newRxn3','metaboliteList',{'prot_pool','prot_A0A024L3D9'},'stoichCoeffList',[-MW1 1],'reversible',false);

% Add gene rules to the reaction
model=changeGeneAssociation(model,'newRxn1','FumC');
model=changeGeneAssociation(model,'newRxn2','FumC');

% Normalization of geneShortNames, enzymes, and enzGenes
model.geneShortNames(1128)={'FumC'};
model.enzymes(964)={'A0A024L3D9'};
model.enzGenes(964)={'FumC'};

% deletion PDC1,5,6 and FUM
model = removeGenes(model,'YGR087C');
model = removeGenes(model,'YLR044C');
model = removeGenes(model,'YLR134W');
model = removeGenes(model,'YPL262W');
% deletion ICL and MAS
model = removeGenes(model,'YER065C');
model = removeGenes(model,'YNL117W');

% Set constraints of reTCA pathway
model=changeRxnBounds(model,'r_1714_REV',15.2,'u');
model=changeRxnBounds(model,'r_2111',0.25,'l');
model=changeRxnBounds(model,'r_1671_REV',1,'u');
model=changeObjective(model,'r_2056');
FBAsolution=optimizeCbModel(model)

% Save model:
cd ../../result_ecYeast/Suc/Models
save ecSuc-reTCA.mat model

% KcatSensitivities analysis
cd ../../../strain_design_ecYeast
load('ecSuc-reTCA.mat')
lyEnzymeSensitivityAnalysis_Suc

% Prediction of ME targets in ecSuc-reTCA model
cd ../../strain_design_ecYeast/scripts_metEngFinders
load('ecSuc-reTCA.mat');
metList = {'Succinate'};
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ecSuc_reTCA_OE.txt');
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','deletion','ecSuc_reTCA_deletion.txt');

% Prediction of ME targets in ecSuc model
cd ../strain_design_ecYeast/scripts_metEngFinders
load('ecYeastGEM_batch.mat');
model = ecModel_batch;
% Set constraints of general Suc pathway
model=changeRxnBounds(model,'r_1714_REV',15.2,'u');
model=changeRxnBounds(model,'r_2111',0.25,'l');
model=changeObjective(model,'r_2056');
FBAsolution=optimizeCbModel(model)
metList = {'Succinate'};
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ecSuc_OE.txt');
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','deletion','ecSuc_deletion.txt');
