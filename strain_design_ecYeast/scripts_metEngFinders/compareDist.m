%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = compareDist(model)
model = getOptimalStrain(model,7,20);
%model = modelModifications(model);
results = FSEOF(model,'FA');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = FSEOF(model,target)
if strcmpi(target,'AA')
    precursors = {'phosphoenolpyruvate','D-erythrose 4-phosphate','malonyl-CoA','NADPH'};
    results.PEP    = getResults(model,precursors(1));
    results.E4P    = getResults(model,precursors(2));
    results.MalCoA = getResults(model,precursors(3));
    results.NADPH  = getResults(model,precursors(4));
    candidates = intersect(results.PEP.genes,results.E4P.genes);
    candidates = intersect(candidates,results.MalCoA.genes);
    candidates = intersect(candidates,results.NADPH.genes);
    results.candidates = candidates;
    [~,~,indexes] = intersect(candidates,results.PEP.genes);
    results.candidatesMat(:,1) = num2cell(results.PEP.k_genes(indexes));
    [~,~,indexes] = intersect(candidates,results.E4P.genes);
    results.candidatesMat(:,2) = num2cell(results.E4P.k_genes(indexes));
    [~,~,indexes] = intersect(candidates,results.MalCoA.genes);
    results.candidatesMat(:,3) = num2cell(results.MalCoA.k_genes(indexes));
    [~,~,indexes] = intersect(candidates,results.NADPH.genes);
    results.candidatesMat(:,4) = num2cell(results.NADPH.k_genes(indexes));
    
elseif strcmpi(target,'FA')
    precursors = {'Acetyl-CoA','malonyl-CoA','NADPH'};
    results.AcCoA    = getResults(model,precursors(1));
    results.MalCoA = getResults(model,precursors(2));
    results.NADPH  = getResults(model,precursors(3));
    candidates = intersect(results.AcCoA.genes,results.AcCoA.genes);
    candidates = intersect(candidates,results.MalCoA.genes);
    candidates = intersect(candidates,results.NADPH.genes);
    results.candidates = candidates;
    [~,~,indexes] = intersect(candidates,results.AcCoA.genes);
    results.candidatesMat(:,1) = num2cell(results.AcCoA.k_genes(indexes));
    [~,~,indexes] = intersect(candidates,results.MalCoA.genes);
    results.candidatesMat(:,2) = num2cell(results.MalCoA.k_genes(indexes));
    [~,~,indexes] = intersect(candidates,results.NADPH.genes);
    results.candidatesMat(:,3) = num2cell(results.NADPH.k_genes(indexes));
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modelModifications(model)
model.ub(strcmp(model.rxns,'r_0670'))     = 0; %H2O + L-kynurenine -> anthranilate + H+ + L-alanine
model.ub(strcmp(model.rxns,'r_1040'))     = 0; %L-threonine  -> acetaldehyde + L-glycine
model.ub(strcmp(model.rxns,'r_1936'))     = 0; %dihydroxyacetone phosphate  -> methylglyoxal + phosphate
model.ub(strcmp(model.rxns,'r_1632_REV')) = 0; %To mitochondrion: acetaldehyde  -> acetaldehyde
model.ub(strcmp(model.rxns,'r_0159')) = 0;  %acetyl-CoA + ethanol    -> coenzyme A + ethyl acetate
model.ub(strcmp(model.rxns,'r_0160')) = 0;  %acetyl-CoA + isoamylol  -> coenzyme A + isoamyl acetate
model.ub(strcmp(model.rxns,'r_0161')) = 0;  %acetyl-CoA + isobutanol -> coenzyme A + isobutyl acetate
[model,~] = addReaction(model,'acetylUsage',{'s_0373','s_0529'},[-1,1],false,0,1000);
[model,~] = addReaction(model,'malonylUsage',{'s_1101','s_0529'},[-1,1],false,0,1000);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = getResults(model,precursor)
positions          = getPrecursorsProdPos(model,precursor);
results            = compare_substrate(model,positions);
geneTable          = getGeneTable(results);
results.geneTable  = geneTable;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function geneTable = getGeneTable(results)
geneTable      = cell(length(results.genes),2);
geneTable(:,1) = results.genes;
geneTable(:,2) = num2cell(results.k_genes);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FC = compare_substrate(model,product)
%Simulate WT (100% growth) and forced (X% growth and the rest towards product):
FC.flux_WT = simulateGrowth(model,product,1);
alpha      = 0.3:0.05:0.8;
v_matrix   = zeros(length(model.rxns),length(alpha));
k_matrix   = zeros(length(model.rxns),length(alpha));
for i = 1:length(alpha)
    FC.flux_MAX   = simulateGrowth(model,product,alpha(i));
    v_matrix(:,i) = FC.flux_MAX;
    k_matrix(:,i) = FC.flux_MAX./FC.flux_WT;
end

%Generate rxn equations:
rxnEqs = printRxnFormula(model,model.rxns,true,true,true);

%Take out rxns with no grRule:
withGR   = ~cellfun(@isempty,model.grRules);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
gene_rxn = model.rxnGeneMat(withGR,:);
FC.rxns  = [model.rxns(withGR) model.rxnNames(withGR) model.grRules(withGR) rxnEqs(withGR)];

%Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
gene_rxn = gene_rxn(non_nan,:);
FC.rxns  = FC.rxns(non_nan,:);

%Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

%Replace any Inf value with 100 (maximum value is ~70):
k_matrix(isinf(k_matrix)) = 100;

%Filter out values that are inconsistent at different alphas:
always_down  = sum(k_matrix <= 1,2) == length(alpha);
always_up    = sum(k_matrix >= 1,2) == length(alpha);
incons_rxns  = always_down + always_up == 0;
incons_genes = sum(gene_rxn(incons_rxns,:),1) > 0;
incons_rxns  = sum(gene_rxn(:,incons_genes),2) > 0;
v_matrix     = v_matrix(~incons_rxns,:);
k_matrix     = k_matrix(~incons_rxns,:);
gene_rxn     = gene_rxn(~incons_rxns,:);
FC.rxns      = FC.rxns(~incons_rxns,:);

%Order from highest to lowest k:
FC.k_rxns   = mean(k_matrix,2);
[~,order]   = sort(FC.k_rxns,'descend');
FC.k_rxns   = FC.k_rxns(order,:);
FC.v_matrix = v_matrix(order,:);
FC.k_matrix = k_matrix(order,:);
gene_rxn    = gene_rxn(order,:);
FC.rxns     = FC.rxns(order,:);

%Create list of remaining genes and filter out any inconsistent score:
FC.genes   = model.genes(sum(gene_rxn,1) > 0);
FC.k_genes = zeros(size(FC.genes));
gene_rxn   = gene_rxn(:,sum(gene_rxn,1) > 0);
cons_genes = false(size(FC.genes));
for i = 1:length(FC.genes)
    k_set         = FC.k_rxns(gene_rxn(:,i) > 0);
    always_down   = sum(k_set <= 1) == length(k_set);
    always_up     = sum(k_set >= 1) == length(k_set);
    cons_genes(i) = always_down + always_up == 1;
    FC.k_genes(i) = mean(k_set);
end
FC.genes   = FC.genes(cons_genes);
FC.k_genes = FC.k_genes(cons_genes);

%Filter any value between mean(alpha) and 1:
unchanged  = (FC.k_genes >= mean(alpha) - 1e-3) + (FC.k_genes <= 1 + 1e-3) == 2;
FC.genes   = FC.genes(~unchanged);
FC.k_genes = FC.k_genes(~unchanged);

%Order from highest to lowest k:
[~,order]  = sort(FC.k_genes,'descend');
FC.genes   = FC.genes(order,:);
FC.k_genes = FC.k_genes(order,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function geneTable = addGenes(geneTable,k_pos,new_set)

N = length(geneTable(:,1));
for i = 1:length(new_set.genes)
    pos = strcmp(geneTable(:,1),new_set.genes{i});
    if sum(pos) > 0
        geneTable{pos,k_pos} = new_set.k_genes(i);
    else
        N = N + 1;
        geneTable{N,1}     = new_set.genes{i};
        geneTable{N,k_pos} = new_set.k_genes(i);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function positions = getPrecursorsProdPos(model,precursors)
positions = [];
for i=1:length(precursors)
    met = precursors(i);
    [~,indexes] = getMetProdIndexes(model,met,true);
    positions   = [positions,indexes];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [metindx,metProduction] = getMetProdIndexes(model,metName,cytFlag)
metProduction  = [];
metindx        = find(strcmpi(model.metNames,metName));
%metindx        = metindx(find(model.metComps(metindx)==cytIndex));
if ~cytFlag
    for i=1:length(metindx)
        row = model.S(metindx(i),:);
        %metProduction  = find(row>0);
        metProduction  = [metProduction,find(row>0)];
    end
else
    cytIndex       = find(strcmpi(model.compNames,'cytoplasm'));
    metindx        = metindx(find(model.metComps(metindx)==cytIndex));
    row            = model.S(metindx,:);
    metProduction  = find(row>0);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%