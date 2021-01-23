function mutantModel = getMutant(model,modifications,message)
% getMutant
%   Get an in-Silico mutant strain based on the provided modifications
%   list. Multiple and combinatorial gene deletions, overexpressions and
%   heterologous enzyme expressions are allowed. Overexpressions are
%   applied directly on enzyme levels (upper bounds) for enzyme-constrained
%   GEMs and on stoichiometric coefficients for non-enzyme related
%   reactions.
%
%   model           Original model (either GEM or ecGEM)
%   modifications   [nx3] Cell array with the desired modifications. Column
%                   (1) contains strings with each of the individual gene IDs
%                   to modify. Column (2) indicates the modification
%                   action; 0 for deletions, 1 for overexpressions and 2
%                   for enzyme activity modification (substitution
%                   endogenous<-heterologous enzyme). Column (3) Over-expression 
%                    factor, 0 for deletions, a given number for overexpressions 
%                    and the new Kcat for heterologous expression (1/s)
%   message         Flag, true if a success message should be displayed
%
%   mutantModel     Mutant model with new constraints
%
%   Usage: mutantModel = getMutant(model,modifications,message)
%
%   Ivan Domenzain.     Last edited 2018/06/27
%           

if nargin<3
    message = false;
end

mutantModel = model;
genes2mod = modifications(:,1);
actions = modifications(:,2);
OE = modifications(:,3);

for i=1:length(genes2mod)
    gene = genes2mod{i};
    action = actions{i};
    OEFactor = OE{i};
    %Deletion mutants
    switch action
        case 0
            mutantModel = removeGenes(mutantModel,gene);
        %Overexpression mutants
        case 1
            gene2modIndex = find(strcmpi(mutantModel.enzGenes,gene));
            %If the gene to overexpress has an associated enzyme then
            %"overexpress" the enzyme usage pseudoRxn Upper bound
            if ~isempty(gene2modIndex)
                enzyme   = mutantModel.enzymes(gene2modIndex);
                %find enzyme exchange or draw reaction
                enzUsage = find(contains(mutantModel.rxnNames,enzyme));
                %If the enzUsage is bounded "overexpress" the UB
                if contains(mutantModel.rxnNames(enzUsage),'exchange')
                    mutantModel.ub(enzUsage) = mutantModel.ub(enzUsage)*OEFactor;
                %For non measured proteins multiply the kinetic coeffs by the OE
                else
                    enzName    = ['prot_' enzyme{1}];
                    enzMetIndx = find(strcmpi(mutantModel.metNames,enzName));
                    enzKcats   = find(mutantModel.S(enzMetIndx,:));
                    enzKcats   = enzKcats(1:end-1);
                    mutantModel.S(enzMetIndx,enzKcats) = mutantModel.S(enzMetIndx,enzKcats)./OEFactor;
                end           
            %If not, Overexpress all the stoichiometric coefficient for all
            %metabolites in the rxns encoded by the gene
            else
                geneRxns = find(contains(mutantModel.grRules,gene));
                if ~isempty(geneRxns)
                    mutantModel.S(:,geneRxns) = mutantModel.S(:,geneRxns).*OEFactor;
                end
            end
        %Change endogenous for Heterologous enzymes
        case 2
            gene2modIndex = find(strcmpi(mutantModel.enzGenes,gene));
            if ~isempty(gene2modIndex)
                enzyme   = mutantModel.enzymes(gene2modIndex);
                enzName    = ['prot_' enzyme{1}];
                enzMetIndx = find(strcmpi(mutantModel.metNames,enzName));
                enzKcats   = find(mutantModel.S(enzMetIndx,:));
                %Avoid enzyme usage reaction
                enzKcats   = enzKcats(1:end-1);
                %Substitute kinetic  coefficients
                mutantModel.S(enzMetIndx,enzKcats) = -1/(3600*OEFactor);%mutantModel.S(enzMetIndx,enzKcats)./OEFactor;
            end
        otherwise
            mutantModel = model;
    end
end
sol = solveLP(mutantModel,1);
if ~isempty(sol.x) && message
    disp(sol.f)
    disp('Mutant strain successfully constructed')
end
end