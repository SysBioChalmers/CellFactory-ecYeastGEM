function [mutantStrain,filtered,step] = robust_ecFSEOF(model,rxnTarget,expYield,CS_MW,resultsFolder)
mkdir(resultsFolder)
current      = pwd;
tol          = 1E-12;
OE           = 2;
thresholds   = [0.5 1];
delLimit     = 0.05;
step         = 0;
essential = readtable('../../ComplementaryData/Essential_ORFs.txt','Delimiter','\t');
essential = strtrim(essential.Ids);

cd GECKO
%Get model parameters
cd geckomat
parameters = getModelParameters;
bioRXN     = parameters.bioRxn;
c_source   = parameters.c_source;
%Parameters for FSEOF method
Nsteps     = 16;
alphaLims  = [0.5*expYield 2*expYield];
%output files for genes and rxn targets
file1   = 'results/genesResults_ecFSEOF.txt';
file2   = 'results/rxnsResults_ecFSEOF.txt';
% 1.- Run FSEOF to find gene candidates
cd utilities/ecFSEOF
mkdir('results')
step = 1;
disp([num2str(step) '.-  **** Running ecFSEOF method (from GECKO utilities) ****'])
results = run_ecFSEOF(model,rxnTarget,c_source,alphaLims,Nsteps,file1,file2);
genes   = results.genes;
disp(' ')
disp(['ecFSEOF yielded ' num2str(length(genes)) ' targets'])
disp(' ')
%Format results table
geneShorts = results.geneNames;
actions    = results.k_genes;
actions(actions<0.5) = 0;
actions(actions>1)   = 1;
MWeigths             = [];
%Identify candidate genes in model enzymes
disp('Extracting enzymatic information for target genes')
[~,iB]     = ismember(genes,model.enzGenes);
candidates = {};
for i=1:numel(iB)
    if iB(i)>0
        candidates = [candidates; model.enzymes(iB(i))];
        MWeigths   = [MWeigths; model.MWs(iB(i))];
    else
        candidates = [candidates; {''}];
        MWeigths   = [MWeigths; nan];
    end
end
disp(' ')
%Get results files structures
candidates = table(genes,candidates,geneShorts,MWeigths,actions,results.k_genes,'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'actions' 'k_scores'});
% Keep top results
disp(['Removing targets ' num2str(thresholds(1)) ' < K_score < ' num2str(thresholds(2))])
toKeep     = find((candidates.k_scores>=thresholds(2)|candidates.k_scores<=thresholds(1)));
candidates = candidates(toKeep,:);
disp([num2str(height(candidates)) ' gene targets remain'])
disp(' ')

% 2.- discard essential genes from deletion targets
step = step+1;
disp([num2str(step) '.-  **** Removing essential genes from KD and KO targets list ****'])
[~,iB]    = ismember(candidates.genes,essential);
toRemove  = iB & candidates.k_scores<=0.05;
candidates(toRemove,:) = [];
cd (current)
disp([num2str(height(candidates)) ' gene targets remain'])
disp(' ')
writetable(candidates,[resultsFolder '/candidates_ecFSEOF.txt'],'Delimiter','\t','QuoteStrings',false);

% 3.- Enzyme usage variability analysis (EVA)
step = step+1;
disp([num2str(step) '.-  **** Running enzyme usage variability analysis ****'])
tempModel = model;
%Get relevant rxn indexes
targetIndx  = find(strcmpi(tempModel.rxns,rxnTarget));
CUR_indx    = find(strcmpi(tempModel.rxnNames,c_source));
growth_indx = find(strcmpi(tempModel.rxns,bioRXN));
prot_indx = find(contains(tempModel.rxns,'prot_pool'));
%Fix suboptimal experimental biomass yield conditions
V_bio = expYield*CS_MW;
tempModel.lb(growth_indx) = V_bio;
%Fix unit C source uptake
tempModel.lb(CUR_indx) = (1-tol)*1;
tempModel.ub(CUR_indx) = (1+tol)*1;
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', targetIndx, +1);
sol       = solveLP(tempModel,1);
WT_prod   = sol.x(targetIndx);
WT_CUR    = sol.x(CUR_indx);
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = (1+tol)*WT_prod;
%Run FVA for all enzyme usages subject to fixed CUR and Grates
FVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);
%sort results
candidateUsages = FVAtable.pU;
%Classify enzyme variability ranges
disp(' ')
disp('Classifying enzyme usage variability ranges')
candidates.EV_type = cell(height(candidates),1);
idxs = find(FVAtable.minU<=tol & FVAtable.maxU<=tol);
candidates.EV_type(idxs) = {'unusable'};
idxs = find(FVAtable.minU>tol);
candidates.EV_type(idxs) = {'essential'};
idxs = find(FVAtable.minU<=tol & FVAtable.maxU>tol);
candidates.EV_type(idxs) = {'totally_variable'};
idxs = find((FVAtable.minU./FVAtable.maxU)>=0.99 & FVAtable.minU>tol);
candidates.EV_type(idxs) = {'tightly_constrained'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & (candidateUsages./FVAtable.maxU)>=0.99);
candidates.EV_type(idxs) = {'variable_optimal'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages<=tol);
candidates.EV_type(idxs) = {'variable_NOToptimal'};
%Classify overexpression types
for i=1:length(candidates.enzymes)
    if FVAtable.maxU(i)~=0 && candidates.actions(i)>0
        candidates.actions(i) = 1;
        %Enzymes that are more tightly constrained are classified as candidates
        %for overexpression by modification on Kcats
        if FVAtable.maxU(i)< OE*candidateUsages(i)
            candidates.actions(i) = 2;
        end 
    end
end
candidates.OE(candidates.actions>0)  = OE;
candidates.OE(candidates.actions==0) = 0;
candidates.minUsage = FVAtable.minU;
candidates.maxUsage = FVAtable.maxU;
candidates.pUsage   = candidateUsages;
%Discard enzymes whose usage LB = UB = 0
disp(' ')
disp('Discard OE targets with lb=ub=0')
toRemove = strcmpi(candidates.EV_type,'unusable') & candidates.actions>0;
candidates(toRemove,:) = [];
disp([num2str(height(candidates)) ' gene targets remain'])
disp(' ')
disp('Discard essential enzymes from deletion targets')
toRemove = (strcmpi(candidates.EV_type,'tightly_constrained') | strcmpi(candidates.EV_type,'essential')) & ...
       (candidates.k_scores<=delLimit);
candidates(toRemove,:) = [];       
disp([num2str(height(candidates)) ' gene targets remain'])
disp(' ')
%Generate table with FVA results
writetable(candidates,[resultsFolder '/candidates_enzUsageFVA.txt'],'Delimiter','\t','QuoteStrings',false);

% 4.- Mechanistic validations of FSEOF results
step = step+1;
disp([num2str(step) '.-  **** Mechanistic validation of results ****'])
%Relevant rxn indexes
relIndexes = [CUR_indx, targetIndx, growth_indx];
%relax target rxn bounds
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = 1000;
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',targetIndx,+1);
%Calculate WT yields
WT_prod_yield = WT_prod/WT_CUR;
%Run mechanistic validation of targets
[FCs_y,FCs_p,validated]  = testAllmutants(candidates,tempModel,relIndexes,WT_prod_yield);
%Discard genes with a negative impact on production yield
candidates.foldChange_yield = FCs_y; 
candidates.foldChange_pRate = FCs_p; 
candidates            = candidates(validated,:);
disp('Discard gene modifications with a negative impact on product yield or rate')
disp([num2str(height(candidates)) ' gene targets remain'])
disp(' ')
writetable(candidates,[resultsFolder '/candidates_mech_validated.txt'],'Delimiter','\t','QuoteStrings',false);
% 4.- Assess genes redundancy
% 
%Get Genes-metabolites network
[GeneMetMatrix,~,Gconect] = getGeneMetMatrix(tempModel,candidates.genes);
%Get independent genes from GeneMetMatrix
[indGenes,G2Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%Append algebraic results to candidates table
candidates.unique = indGenes;
candidates.conectivity = Gconect.mets_number;
[~,groups]        = getGenesGroups(G2Gmatrix,candidates.genes);
candidates.groups = groups;
% Rank candidates by priority
disp('Ranking gene targets by priority level:')
disp('  1.- Unique genes candidates for OE with pUsage>0')
disp('  1.- Unique genes candidates for del with pUsage=0 and maxUsage>0')
disp('  2.- Isoenzymes candidates for OE with the lowest MW for a metabolic rxn')
disp('  2.- Groups of isoenzymes candidates for deletion with all pUsage=0')
disp('  3.- The heaviest enzymes in a group of isoenzymes candidates for deletion and at least one pUsage>0')
disp(' ')
priority = zeros(height(candidates),1);
%%% 1st. unique=1 OEs with both min and pUsage>0 & Deletions with pUsage=0
%unique Enzymes that are necesarily used
cond1 = (candidates.actions>0 & candidates.minUsage>0);
%unique enzymes that are not used in a parsimonious simulation
cond2   = (candidates.actions==0 & candidates.pUsage==0);
indexes = (candidates.unique==1 & (cond2 | cond1));
priority(indexes) = 1; 
%%% 2nd. unique=0, for OEs pick the enzyme with the lowest MW for each group, these 
% are usually isoenzymes, therefore the lowest protein burden impact is
% desired. For deletions assign 2 to groups of isoenzymes in which none of
% them is used. For groups of isoenzymes candidates for deletions in which 
% there are used enzymes in the parsimonious distribution then delete the
% non-used ones.
for i=1:max(candidates.groups)
    %Find group genes
    groupIndxs = find(candidates.groups==i);
    groupTable = candidates(groupIndxs,:); 
    if sum(candidates.actions(groupIndxs))==0
        nonZeros = (candidates.pUsage(groupIndxs)>0);
        if sum(nonZeros)==0
            priority(groupIndxs) = 2;
        else
            priority(groupIndxs(~nonZeros)) = 3;
        end
    else
        if all(candidates.actions(groupIndxs)>0)
            %Select those with pUsage>0  and minUsage>0 and the lowest MW
            groupTable = groupTable((groupTable.pUsage>0),:);
            [~,minMW]  = min(groupTable{:,'MWs'});
            groupGenes = groupTable.genes(minMW);
            if ~isempty(groupGenes)                
                prtyIndx   = (strcmpi(candidates.genes,groupGenes));
                priority(prtyIndx) = 2;
            end
        end
    end
end
candidates.priority = priority;
%Keep priority genes and sort them accordingly
candidates = candidates(priority>0,:);
disp('Discard genes with priority level = 0')
disp([num2str(height(candidates)) ' gene targets remain'])


% 5.- Find compatible combinations
step = 5;
disp(num2str(step))
disp(' ')
candidates = sortrows(candidates,'priority','ascend');
writetable(candidates,[resultsFolder '/candidates_priority.txt'],'Delimiter','\t','QuoteStrings',false);
% get optimal strain according to priority candidates
disp('Constructing optimal strain')
[mutantStrain,filtered] = getOptimalStrain(tempModel,candidates,[targetIndx CUR_indx prot_indx],false);
step = 6;
disp(num2str(step))
%[mutantStrain,filtered,] = getOptimalStrain(tempModel,filtered,[targetIndx CUR_indx prot_indx],false);
%step = 7;
cd (current)
actions = cell(height(filtered),1);
actions(filtered.actions==0)= {'deletion'};
actions(filtered.actions>0) = {'OE'};
filtered.actions = actions;
disp([num2str(height(filtered)) ' gene targets remain'])
disp(' ')
writetable(filtered,[resultsFolder '/compatible_genes_results.txt'],'Delimiter','\t','QuoteStrings',false);
origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
copyfile(origin,resultsFolder)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FChanges_y,FChanges_p,validated] = testAllmutants(candidates,tempModel,indexes,tol)
if nargin<4
    tol = 0.001;
end
FChanges_y = [];
FChanges_p = [];
CUR_indx    = indexes(1);
targetIndx  = indexes(2);
growth_indx = indexes(3);
%Index to minimize (bi-level optimization)
minIndex = find(contains(tempModel.rxnNames,'prot_pool'));
%Unconstrain CUR and biomass formation
tempModel = setParam(tempModel,'ub',CUR_indx,1000);
tempModel = setParam(tempModel,'lb',CUR_indx,0);
tempModel = setParam(tempModel,'ub',growth_indx,1000);
tempModel = setParam(tempModel,'lb',growth_indx,0);
tempModel = setParam(tempModel,'obj',growth_indx,1);
%Get WT solutions
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',minIndex);
maxGrowth = WTsol(growth_indx);
tempModel = setParam(tempModel,'lb',growth_indx,0.9*maxGrowth);
tempModel = setParam(tempModel,'obj',targetIndx,1);
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',minIndex);
prodWT    = WTsol(targetIndx);
WTval     = prodWT/WTsol(CUR_indx);

for i=1:height(candidates)
    gene   = candidates.genes{i};
    enzyme = candidates.enzymes{i};
    short  = candidates.shortNames{i};
    action = candidates.actions(i);
    if action ==1
        action = 2;
    end
    OEf    = candidates.OE(i);
    modifications = {gene action OEf};
    if ~isempty(enzyme)
        index  = find(contains(tempModel.rxnNames,['draw_prot_' enzyme]));
        pUsage = WTsol(index);
    else 
        pUsage = [];
    end
    if action ==0 & isempty(pUsage)
        pUsage = 1E-9;
    end
    mutantModel     = getMutant(tempModel,modifications,pUsage);
    [mutSolution,~] = solveECmodel(mutantModel,mutantModel,'pFBA',minIndex);
    if ~isempty(mutSolution)
        yield = mutSolution(targetIndx)/mutSolution(CUR_indx);
        FC_y  = yield/WTval;
        FC_p  = mutSolution(targetIndx)/prodWT;
    else
        FC_y = 0;
        FC_p = 0;
    end
    FChanges_y = [FChanges_y; FC_y];
    FChanges_p = [FChanges_p; FC_p];
    disp(['Ready with genetic modification #' num2str(i) '[' short ': ' num2str(action) '] FC: ' num2str(FC_y)])
end
validated  = mean([FChanges_y,FChanges_p])>=1-tol;
[maxVal,I] = max(FChanges_y);
if ~(maxVal>=1)
    maxVal = [];
    gene   = [];
else 
    TOPgene = candidates.genes{I};
    FC_y    = FChanges_y(I);
    disp(['candidate gene: ' short ' FC: ' num2str(FC_y)])
end
end

