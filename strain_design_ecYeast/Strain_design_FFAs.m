% Pathway construction of FFAs in yeastGEM model

% Add new reactions to the model
cd ../ComplementaryScripts
model = loadYeastModel;

% Load stoichiometry data:
fid = fopen('../ComplementaryData/FFAs/newpathway_newRxnMatrix.tsv');
newreaction = textscan(fid,'%s %s %s %s %s','Delimiter','\t','HeaderLines',1);
matrix.rxnIDs  = newreaction{1};
matrix.metcoef = cellfun(@str2num, newreaction{2});
matrix.metIDs  = newreaction{3};
matrix.mettype = newreaction{4};
matrix.metcompartments = newreaction{5};
fclose(fid);

% Load rxn properties data:
fid  = fopen('../ComplementaryData/FFAs/newpathway_newRxnProp.tsv','r');
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
fid = fopen('../../ComplementaryData/FFAs/newpathway_newRxnMetAnnotation.tsv');
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
cd ../../result_ecYeast/FFAs/Models
save FFAs_GEM.mat model

% The theoretical yield analysis of FFAs in yeastGEM model
% Set constraints
model=changeRxnBounds(model,'r_1714',-10,'l'); % minimal glucose uptake rate
model=changeRxnBounds(model,'r_2111',0.03,'l'); % maximum growth rate
model=changeObjective(model,'r_4602'); % RxnID of new product's exchange reaction
FBAsolution=optimizeCbModel(model)

% Pathway construction of FFAs in ecYeast model

% Add new reactions of asp pathway to the model
cd ../../../ModelFiles/mat
load('ecYeastGEM_batch.mat');
model = ecModel_batch;
model=changeRxnBounds(model,'r_1714_REV',10,'u');
model=changeRxnBounds(model,'r_2111',0.03,'l');
model = addReaction(model,'newRxn1','metaboliteList',{'s_0597','s_1067','s_1163','s_1288','s_1295','s_1450','s_2826','s_4209'},'stoichCoeffList',[-1 -1 -1 -1 -1 -1 -1 1],'reversible',false);
model = addReaction(model,'newRxn2','metaboliteList',{'s_4209'},'stoichCoeffList',[-1],'reversible',false);
model=changeObjective(model,'newRxn2');
FBAsolution=optimizeCbModel(model)

% add rxns of tesA
Kcat1=10.13*3600
MW1=23.622
model = addReaction(model,'newRxn10','metaboliteList',{'s_1302','s_0803','prot_P0ADA1','s_0529','s_1286','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn11','metaboliteList',{'s_1454','s_0803','prot_P0ADA1','s_0529','s_1449','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn12','metaboliteList',{'s_1262','s_0803','prot_P0ADA1','s_0529','s_1260','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn13','metaboliteList',{'s_4210','s_0803','prot_P0ADA1','s_0529','s_0595','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn14','metaboliteList',{'s_1073','s_0803','prot_P0ADA1','s_0529','s_1065','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn15','metaboliteList',{'s_1176','s_0803','prot_P0ADA1','s_0529','s_1161','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn16','metaboliteList',{'s_2877','s_0803','prot_P0ADA1','s_0529','s_1293','s_0794'},'stoichCoeffList',[-1 -1 -1/Kcat1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn20','metaboliteList',{'prot_pool','prot_P0ADA1'},'stoichCoeffList',[-MW1 1],'reversible',false)

% add rxn of ME
model = addReaction(model,'newRxn3','metaboliteList',{'s_0068','s_1200','s_0460','s_1205','s_1401'},'stoichCoeffList',[-1 -1 1 1 1],'reversible',false)

% add rxn of ACL
Kcat2=27.1383*3600
MW2=119.728
model = addReaction(model,'newRxn4','metaboliteList',{'s_0434','s_0522','s_0529','prot_Q91V92','s_0373','s_0394','s_1271','s_1322'},'stoichCoeffList',[-1 -1 -1 -1/Kcat2 1 1 1 1],'reversible',false)
model = addReaction(model,'newRxn17','metaboliteList',{'prot_pool','prot_Q91V92'},'stoichCoeffList',[-MW2 1],'reversible',false)

% add rxn of FAS
Kcat3=4.16*3600
MW3=318.115
MW4=137.582
model = addReaction(model,'newRxn5','metaboliteList',{'s_0373','s_0794','s_1101','s_1212','prot_A0A2Z6EZ44','prot_A0A2Z6EYW1','s_0456','s_0529','s_0803','s_1207','s_1302'},'stoichCoeffList',[-1 -21 -7 -14 -1/Kcat3 -1/Kcat3 7 7 7 14 1],'reversible',false)
model = addReaction(model,'newRxn6','metaboliteList',{'s_0373','s_0794','s_1101','s_1212','prot_A0A2Z6EZ44','prot_A0A2Z6EYW1','s_0456','s_0529','s_0803','s_1207','s_1454'},'stoichCoeffList',[-1 -24 -8 -16 -1/Kcat3 -1/Kcat3 8 8 8 16 1],'reversible',false)
model = addReaction(model,'newRxn7','metaboliteList',{'s_0373','s_0794','s_1101','s_1212','prot_A0A2Z6EZ44','prot_A0A2Z6EYW1','s_0456','s_0529','s_0803','s_1207','s_1176'},'stoichCoeffList',[-1 -18 -6 -12 -1/Kcat3 -1/Kcat3 6 6 6 12 1],'reversible',false)
model = addReaction(model,'newRxn8','metaboliteList',{'s_0373','s_0794','s_1101','s_1212','prot_A0A2Z6EZ44','prot_A0A2Z6EYW1','s_0456','s_0529','s_0803','s_1207','s_1073'},'stoichCoeffList',[-1 -15 -5 -10 -1/Kcat3 -1/Kcat3 5 5 5 10 1],'reversible',false)
model = addReaction(model,'newRxn9','metaboliteList',{'s_0373','s_0794','s_1101','s_1212','prot_A0A2Z6EZ44','prot_A0A2Z6EYW1','s_0456','s_0529','s_0803','s_1207','s_4210'},'stoichCoeffList',[-1 -12 -4 -8 -1/Kcat3 -1/Kcat3 4 4 4 8 1],'reversible',false)
model = addReaction(model,'newRxn18','metaboliteList',{'prot_pool','prot_A0A2Z6EZ44'},'stoichCoeffList',[-MW3 1],'reversible',false)
model = addReaction(model,'newRxn19','metaboliteList',{'prot_pool','prot_A0A2Z6EYW1'},'stoichCoeffList',[-MW4 1],'reversible',false)

% deletion FAA1,4
model = removeGenes(model,'YOR317W')
model = removeGenes(model,'YMR246W')
% deletion POX1
model = removeGenes(model,'YGL205W')
% deletion HFD1
model = removeGenes(model,'YMR110C')
% deletion of MDH3
model = removeGenes(model,'YDL078C')

% GeneAssociation
model=changeGeneAssociation(model,'newRxn3','ME')
model=changeGeneAssociation(model,'newRxn4','ACL')
model=changeGeneAssociation(model,'newRxn5','FAS1 and FAS2')
model=changeGeneAssociation(model,'newRxn6','FAS1 and FAS2')
model=changeGeneAssociation(model,'newRxn7','FAS1 and FAS2')
model=changeGeneAssociation(model,'newRxn8','FAS1 and FAS2')
model=changeGeneAssociation(model,'newRxn9','FAS1 and FAS2')
model=changeGeneAssociation(model,'newRxn10','tesA')
model=changeGeneAssociation(model,'newRxn11','tesA')
model=changeGeneAssociation(model,'newRxn12','tesA')
model=changeGeneAssociation(model,'newRxn13','tesA')
model=changeGeneAssociation(model,'newRxn14','tesA')
model=changeGeneAssociation(model,'newRxn15','tesA')
model=changeGeneAssociation(model,'newRxn16','tesA')

% Normalization of geneShortNames, metComps, enzymes, and enzGenes
model.geneShortNames(1123)={'ME'};
model.geneShortNames(1124)={'ACL'};
model.geneShortNames(1125)={'FAS1'};
model.geneShortNames(1126)={'FAS2'};
model.geneShortNames(1127)={'tesA'};

model.metComps(4147)=3;
model.metComps(4148)=1;
model.metComps(4149)=1;
model.metComps(4150)=1;
model.metComps(4151)=1;
model.metComps(4152)=1;

model.enzymes(964)={'P0ADA1'};
model.enzymes(965)={'Q91V92'};
model.enzymes(966)={'A0A2Z6EZ44'};
model.enzymes(967)={'A0A2Z6EYW1'};

model.enzGenes(964)={'tesA'};
model.enzGenes(965)={'ACL'};
model.enzGenes(966)={'FAS1'};
model.enzGenes(967)={'FAS2'};

% Save model:
cd ../../result_ecYeast/FFAs/Models
save ecFFAs.mat model

% KcatSensitivities analysis
load('ecFFAs.mat')
cd ../../../strain_design_ecYeast
lyEnzymeSensitivityAnalysis_FFAs

% Prediction of ME targets in ecFFAs model
cd ../../strain_design_ecYeast/scripts_metEngFinders
metList = {'s_4209'};
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ecFFAs_OE.txt');
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','OE','ecFFAs_downregulation.txt',0.05);
ResultsTable = lymetEng_TargetsFinder(model,'D-glucose exchange (reversible)',metList,'production','extracellular','deletion','ecFFAs_deletion.txt');
