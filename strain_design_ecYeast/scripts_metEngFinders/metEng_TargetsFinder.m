function ResultsTable = metEng_TargetsFinder(model,c_sourceID,metList,direction,compartment,action,filename)
% metEng_TargetsFinder
%
% Function that computes the effect of all the possible gene deletions or 
% overexpressions in the total production or consumption yields for a given
% set of internal metabolites in the network.
% 
% model         a model structure
% c_sourceID    Rxn ID for the carbon source uptake reaction
% metList       Cell array containing the metName for all the metabolites to
%               analyze.
% direction     String 'production' when the total flux that produces a given
%               metabolite is required, 'consumption' when the total flux that 
%               consumes a metabolite is required instead.
% compartment   String that indicates the model compartment in which the total 
%               flux calculations should be done
% action        'deletion' for deleting all the genes in the model (one at a
%               time) and 'OE' for overexpressions with a default OE factor of 20.
% filename      String indicating the name of the file (with extension) in
%               which the numeric results are saved.
%
% ResultsTable  Matrix [# of genes x number of tracked metabolites+1].
%               Contains the fold-changes (1-FC) for the total production of 
%               the metabolites of interest for each simulated mutant
%               strain.
%
% Created.  Ivan Domenzain 2018-10-15
%

nMets = length(metList);
%Identify target metabolites and the rxns that produce them on the desired
%compartment
metIndexes = cell(0,nMets);
stoich     = cell(0,nMets);
for i=1:nMets
    [~,metIndexes{i},stoich{i}] = getMetProdIndexes(model,metList{i},direction,compartment);
    disp([metList{i} ' ' direction ' reactions:'])
    disp(model.rxnNames(metIndexes{i}))
end

%preallocate results matrixes
m = length(model.genes)+1;
resultsMat = zeros(m,nMets+1);
%Set culture media
[model,~] = changeMedia_batch(model,c_sourceID,'Min');
%Get all protein usage indexes (if model is enzyme-constrained)
if isfield(model,'enzymes') && isfield(model,'enzGenes')
    prots = find(contains(model.rxnNames,'prot'));
    %Exclude the protein pool
    prots = prots(1:end-1);
else 
    prots = [];
end
%Find the glucose uptake rate
glucUptkIndx = find(strcmpi(model.rxnNames,c_sourceID));
%Objective index
objIndex   = find(model.c);
%Get Wild type solution 
[base_sol,~] = solveMutant(model,[],'pFBA',prots);
glucUptake   = base_sol(glucUptkIndx)
WTgrowth     = base_sol(objIndex)

%Calculate WT yield for selected metabolites
indexes = {};
WT_yields = zeros(1,nMets);
for i=1:nMets
    %get the total flux towards/from the i-th metabolite
    totalFlux    = sum(base_sol(metIndexes{i}));
    indexes      = [indexes;metIndexes(i)];
    WT_yields(i) = totalFlux/glucUptake;
    disp(['WT ' metList{i} ' ' direction ' yield: ' num2str(WT_yields(i)) ' [mol/mol glucose]'])
end
WT_yields(WT_yields==0) = 1E-6;
indexes     = [indexes;{glucUptkIndx};{objIndex}];
growthYield = WTgrowth/(glucUptake*0.180);
disp(['WT Growth yield: ' num2str(growthYield) ' [gBiomass/gGlucose]'])
fprintf('\n')
%Loop through all the original genes
for j=1:length(model.genes)    
    gene = model.genes(j);
    %Get met production yields for every  mutant
    [Fchanges,successWT,mutGrowth] = getResultsForGene(model,model,gene,'pFBA',WT_yields,indexes,stoich,action,prots);
    %save results
    gRateFC         = mutGrowth/WTgrowth;
    resultsMat(j,:) = [Fchanges,gRateFC];
    disp(['Ready with gene ' num2str(j) ' feasible: ' num2str(successWT) ' average FC: ' num2str(Fchanges(end))])
end
%Rearrange results WT yields in the first row and each mutant in the rest
%of the rows
resultsMat(2:end,:) = resultsMat(1:end-1,:);
resultsMat(1,:)     = [WT_yields,growthYield];
resultsMat = num2cell(resultsMat);
metList    = strrep(metList,'-','_');
metList    = strrep(metList,' ','_');
metList    = [metList,'Growth'];
Rows       = ['WT_yields';model.genes];
ResultsTable = cell2table(resultsMat,'VariableNames',metList,'RowNames',Rows);
writetable(ResultsTable,filename,'WriteVariableNames',true,'WriteRowNames',true,'Delimiter','\t')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [results,success,gRate] = getResultsForGene(strainModel,model,gene,method,WT_yields,indexes,stoich,action,prots)
success = 0;
nMets   = length(WT_yields);
objIndex  = indexes{end};
glucUptkIndx  = indexes{end-1};

results = zeros(1,nMets); 
gRate   = 0;
%find gene in model's genes
if ~isempty(find(strcmpi(strainModel.genes,gene), 1)) || strcmpi(gene,'') 
    if strcmpi(action,'deletion')
        %Get mutant (single deletion)
         mutant = removeGenes(model,gene);     
    elseif strcmpi(action,'OE')
        %Get mutant (single OE)
        mutant = getMutant(strainModel,{gene,1,20});
    end
    %optimize
    [solution,flag] = solveMutant(mutant,model,method,prots);   
    %if the optimization was feasible, then get results
    if flag ==1 
        %glucose uptake rate
        mutUptake = solution(glucUptkIndx);
        if mutUptake>0
            success = 1;
            %get yield fold-changes for each precursor
            for i=1:nMets
                results(1,i) = getYieldFoldChange(indexes{i},solution,mutUptake,WT_yields(i),stoich{i});
            end
            gRate = solution(objIndex);
        end
    end
else
    [solution,~] = solveMutant(model,model,method,prots);
    results      = ones(1,nMets);
    gRate        = solution(objIndex);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function YieldFC = getYieldFoldChange(prodIndex,mutSolution,Cuptake,WT_Yield,stoich)
%Multiply by the stoichiometric coefficients
product =mutSolution(prodIndex).*stoich;
mutProdYield = sum(product)/Cuptake;
%Get metabolite production yields ratio betweeen WT and
%double mutants
YieldFC = mutProdYield/WT_Yield;
end



            


