% Pathway construction of 3HP in yeastGEM model

% Add new reactions to the model
cd ../ComplementaryScripts
model = loadYeastModel;

% Load stoichiometry data:
fid = fopen('../ComplementaryData/3HP/newpathway_newRxnMatrix.tsv');
newreaction = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
matrix.rxnIDs  = newreaction{1};
matrix.metcoef = cellfun(@str2num, newreaction{2});
matrix.metIDs  = newreaction{3};
matrix.mettype = newreaction{4};
matrix.metcompartments = newreaction{5};
fclose(fid);

% Load rxn properties data:
fid  = fopen('../ComplementaryData/3HP/newpathway_newRxnProp.tsv','r');
rev = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newrxn.ID  = rev{1};
newrxn.Rev = cellfun(@str2num, rev{2});
newrxn.GPR = rev{3};
newrxn.rxnNames      = rev{4};
newrxn.rxnECNumbers  = rev{5};
newrxn.rxnKEGGID     = rev{6};
newrxn.rxnNotes      = rev{7};
newrxn.rxnMetaNetXID = newrxn.ID;
for i = 1:length(newrxn.rxnMetaNetXID)
    if ~startsWith(newrxn.rxnMetaNetXID{i},'MNXR')
        newrxn.rxnMetaNetXID{i} = '';
    end
end
newrxn.rxnMetaNetXID = regexprep(newrxn.rxnMetaNetXID,'_cv','');
fclose(fid);

% Change coefficients for reactants:
for i=1:length(matrix.rxnIDs)
    if strcmp(matrix.mettype(i),'reactant')
        matrix.metcoef(i) = matrix.metcoef(i)*-1;
    end
end

% Change compartments:
CONValldata = cat(2,model.compNames,model.comps);
lbracket    = ' [' ;
llbracket   = '[';
rbrackets   = ']';
space       = ' ';
[m, n]      = size(CONValldata);
for i = 1:m
    aa = CONValldata(i,1);
    aa = char(aa);
    for j=1:length(matrix.rxnIDs)
        bb = matrix.metcompartments(j,1);
        bb = char(bb);
        if strcmp(bb,aa)
            matrix.Newcomps(j,1) = CONValldata(i,2);
        end
    end
end
for i=1:length(matrix.rxnIDs)
    matrix.metnames(i) = strcat(matrix.metIDs(i),lbracket,matrix.metcompartments(i),rbrackets);
    matrix.Newcomps(i) = strcat(llbracket,matrix.Newcomps(i),rbrackets);
end

% Map mets to model.metnames, get s_index for new mets:
cd otherChanges
for j = 1:length(matrix.metnames)
    [~,metindex] = ismember(matrix.metnames(j),model.metNames);
    if metindex ~= 0
        matrix.mets(j) = model.mets(metindex);
    elseif metindex == 0
        newID = getNewIndex(model.mets);
        matrix.mets(j) = strcat('s_',newID,matrix.Newcomps(j));
        model = addMetabolite(model,char(matrix.mets(j)), ...
            'metName',matrix.metnames(j));
    end
end

% Add metabolite data:
fid = fopen('../../ComplementaryData/3HP/newpathway_newRxnMetAnnotation.tsv');
newmet_annot         = textscan(fid,'%s %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
newmet.metNames      = newmet_annot{1};
newmet.metFormulas   = newmet_annot{2};
newmet.metCharges    = cellfun(@str2num, newmet_annot{3});
newmet.metKEGGID     = newmet_annot{5};
newmet.metChEBIID    = newmet_annot{6};
newmet.metMetaNetXID = newmet_annot{7};

fclose(fid);
for i = 1:length(newmet.metNames)
    [~,metID] = ismember(newmet.metNames(i),model.metNames);
    if metID ~= 0
        model.metFormulas{metID}   = newmet.metFormulas{i};
        model.metCharges(metID)    = newmet.metCharges(i);
        model.metKEGGID{metID}     = newmet.metKEGGID{i};
        model.metChEBIID{metID}    = newmet.metChEBIID{i};
        model.metMetaNetXID{metID} = newmet.metMetaNetXID{i};
        model.metNotes{metID}      = 'NOTES: added after the heterologous update';
    end
end

% Add new reactions according to rev ID: Met Coef needs to be a column, not
% a row. Coef should be a double, which was converted at the import section
EnergyResults     = {};
MassChargeresults = {};
RedoxResults      = {};
if ~isfield(model,'rxnMetaNetXID')
    model.rxnMetaNetXID = cell(size(model.rxns));
end
for i = 1:length(newrxn.ID)
    cd ../otherChanges
    newID = getNewIndex(model.rxns);
    j     = find(strcmp(matrix.rxnIDs,newrxn.ID{i}));
    Met   = matrix.mets(j);
    Coef  = transpose(matrix.metcoef(j));
    [model,rxnIndex] = addReaction(model, ['r_' newID],...
        'reactionName', newrxn.ID{i},...
        'metaboliteList',Met,...
        'stoichCoeffList',Coef,...
        'reversible',newrxn.Rev(i,1),...
        'geneRule',newrxn.GPR{i},...
        'checkDuplicate',1);
    cd ../modelexpansion/
    [EnergyResults,RedoxResults] = CheckEnergyProduction(model,{['r_' newID]},EnergyResults,RedoxResults);
    [MassChargeresults] = CheckBalanceforSce(model,{['r_' newID]},MassChargeresults);
    if isempty(rxnIndex)
        rxnIndex = strcmp(model.rxns,['r_' newID]);
    end
    % Add rxn annotation:
    model.rxnNames{rxnIndex}      = newrxn.rxnNames{i};
    model.rxnECNumbers(rxnIndex)  = newrxn.rxnECNumbers(i);
    model.rxnKEGGID(rxnIndex)     = newrxn.rxnKEGGID(i);
    model.rxnMetaNetXID(rxnIndex) = newrxn.rxnMetaNetXID(i);
    model.rxnConfidenceScores(rxnIndex) = 1;   %reactions without gene but needed for modelling
    model.rxnNotes{rxnIndex} = ['NOTES: heterologous pathway; ',newrxn.rxnNotes{i}];
end

% Save model:
cd ../../result_ecYeast/3HP/Models
save 3HP_GEM.mat model

% The theoretical yield analysis of different pathways of 3HP
% Set constraints of asp pathway
model=changeRxnBounds(model,'r_4601',0,'u');
model=changeRxnBounds(model,'r_4607',0,'u'); % Block other added pathways by setting the 'ub' of 0 when analysing a specific pathway
model=changeRxnBounds(model,'r_1714',-10,'l'); % minimal glucose uptake rate
model=changeRxnBounds(model,'r_2111',0.1,'l'); % maximum growth rate
model=changeObjective(model,'r_4604'); % RxnID of new product's exchange reaction
FBAsolution=optimizeCbModel(model)

% Set constraints of malcoa pathway
model=changeRxnBounds(model,'r_4607',1000,'u');
model=changeRxnBounds(model,'r_4601',0,'u');
model=changeRxnBounds(model,'r_4602',0,'u');
model=changeRxnBounds(model,'r_4605',0,'u');
model=changeRxnBounds(model,'r_4606',0,'u'); % Block other added pathways by setting the 'ub' of 0 when analysing a specific pathway
model=changeRxnBounds(model,'r_1714',-10,'l'); % minimal glucose uptake rate
model=changeRxnBounds(model,'r_2111',0.1,'l'); % maximum growth rate
model=changeObjective(model,'r_4604'); % RxnID of new product's exchange reaction
FBAsolution=optimizeCbModel(model)

% Pathway construction of 3HP in ecYeast model

% Add new reactions of asp pathway to the model
cd ../../../ModelFiles/mat
load('ecYeastGEM_batch.mat');
model = ecModel_batch;

% Add rxns of gene fdyG
Kcat1=115*3600;
MW1=27.249;
model = addReaction(model,'newRxn1','metaboliteList',{'s_0794','s_1212','s_4184','prot_P39831','s_1207','s_4207'},'stoichCoeffList',[-1 -1 -1 -1/Kcat1 1 1],'reversible',false);
model = addReaction(model,'newRxn7','metaboliteList',{'prot_pool','prot_P39831'},'stoichCoeffList',[-MW1 1],'reversible',false);
% Add rxns of gene bcere0029_32090
model = addReaction(model,'newRxn4','metaboliteList',{'s_0441','s_1399','s_0955','s_4184'},'stoichCoeffList',[-1 -1 1 1],'reversible',false);
% Add rxns of gene A7U8C7
Kcat2=7.03*3600;
MW2=61.24;
model = addReaction(model,'newRxn5','metaboliteList',{'s_0794','s_0973','prot_A7U8C7','s_0441','s_0456'},'stoichCoeffList',[-1 -1 -1/Kcat2 1 1],'reversible',false);
model = addReaction(model,'newRxn8','metaboliteList',{'prot_pool','prot_A7U8C7'},'stoichCoeffList',[-MW2 1],'reversible',false);
% Add rxns of gene YGR019Wly
Kcat3=0.1324*3600;
MW3=52.946;
model = addReaction(model,'r_4572','metaboliteList',{'s_0180','s_0441','prot_P17649ly','s_0991','s_4184'},'stoichCoeffList',[-1 -1 -1/Kcat3 1 1],'reversible',false);
model = addReaction(model,'newRxn9','metaboliteList',{'prot_pool','prot_P17649ly'},'stoichCoeffList',[-MW3 1],'reversible',false);
% Add transport and exchange rxns of 3HP
model = addReaction(model,'newRxn2','metaboliteList',{'s_4207','s_4208'},'stoichCoeffList',[-1 1],'reversible',false);
model = addReaction(model,'newRxn3','metaboliteList',{'s_4208'},'stoichCoeffList',[-1],'reversible',false);

% Add gene rules to the reaction
model=changeGeneAssociation(model,'newRxn1','fdyG');
model=changeGeneAssociation(model,'newRxn4','bcere0029_32090');
model=changeGeneAssociation(model,'newRxn5','A7U8C7');
model=changeGeneAssociation(model,'r_4572','YGR019Wly');

% Normalization of geneShortNames, metComps, enzymes, and enzGenes
model.geneShortNames(1128)={'fdyG'};
model.geneShortNames(1129)={'bcere0029_32090'};
model.geneShortNames(1130)={'A7U8C7'};
model.geneShortNames(1131)={'YGR019Wly'};

model.metComps(4147)=1;
model.metComps(4148)=1;
model.metComps(4149)=1;
model.metComps(4150)=1;
model.metComps(4151)=3;

model.enzymes(964)={'P39831'};
model.enzymes(965)={'A7U8C7'};
model.enzymes(966)={'P17649ly'};

model.enzGenes(964)={'fdyG'};
model.enzGenes(965)={'A7U8C7'};
model.enzGenes(966)={'YGR019Wly'};

% The theoretical yield analysis of asp pathway of 3HP
% Set constraints of asp pathway
model=changeRxnBounds(model,'r_1714_REV',10,'u');
model=changeRxnBounds(model,'r_2111',0.1,'l');
model=changeObjective(model,'newRxn3');
FBAsolution=optimizeCbModel(model)

% Save model:
cd ../../result_ecYeast/3HP/Models
save ec3HP-asp.mat model

% Add new reactions of malcoa pathway to the model
cd ../../../ModelFiles/mat
load('ecYeastGEM_batch.mat');
model = ecModel_batch;

% Add rxns of gene MCRca
Kcat1=50*3600;
MW1=300;
model = addReaction(model,'newRxn10','metaboliteList',{'s_0794','s_1101','s_1212','prot_Q6QQP7','s_0529','s_1207','s_4207'},'stoichCoeffList',[-2 -1 -2 -1/Kcat1 1 2 1],'reversible',false)
model = addReaction(model,'newRxn11','metaboliteList',{'prot_pool','prot_Q6QQP7'},'stoichCoeffList',[-MW1 1],'reversible',false);

% Add transport and exchange rxns of 3HP
model = addReaction(model,'newRxn2','metaboliteList',{'s_4207','s_4208'},'stoichCoeffList',[-1 1],'reversible',false);
model = addReaction(model,'newRxn3','metaboliteList',{'s_4208'},'stoichCoeffList',[-1],'reversible',false);

% Add gene rules to the reaction
model=changeGeneAssociation(model,'newRxn10','MCRca');

% Normalization of geneShortNames, metComps, enzymes, and enzGenes
model.geneShortNames(1128)={'MCRca'};
model.metComps(4147)=1;
model.metComps(4148)=3;
model.metComps(4149)=1;
model.enzymes(964)={'Q6QQP7'};
model.enzGenes(964)={'MCRca'};

% The theoretical yield analysis of malcoa pathway of 3HP
% Set constraints of malcoa pathway
model=changeRxnBounds(model,'r_1714_REV',10,'u');
model=changeRxnBounds(model,'r_2111',0.1,'l');
model=changeObjective(model,'newRxn3');
FBAsolution=optimizeCbModel(model)

% Save model:
cd ../../result_ecYeast/3HP/Models
save ec3HP-malcoa.mat model

% KcatSensitivities analysis
cd ../../../strain_design_ecYeast
load('ec3HP-asp.mat')
lyEnzymeSensitivityAnalysis_3HP

% Prediction of ME targets in ec3HP-asp model
cd ../../strain_design_ecYeast/scripts_metEngFinders
load('ec3HP-asp.mat');
metList = {'s_4208'};
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ec3HP_asp_OE.txt');
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','deletion','ec3HP_asp_deletion.txt');

% Prediction of ME targets in ec3HP-malcoa model
cd ../strain_design_ecYeast/scripts_metEngFinders
load('ec3HP-malcoa.mat');
metList = {'s_4208'};
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ec3HP_malcoa_OE.txt');
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','deletion','ec3HP_malcoa_deletion.txt');