function [optStrain,optGenes,FChanges,iB] = getOptimalStrain(model,candidates,rxnIndxs,protFlag)
%Get WT yield
targetIndex = rxnIndxs(1);
GURindex    = rxnIndxs(2);
prot_index  = rxnIndxs(3);
[sol,flag]  = solveECmodel(model,model,'pFBA',prot_index);
if flag
    if protFlag
        minIndex = prot_index;
    else
        minIndex = GURindex;
    end
    WTyield = sol(targetIndex)/(sol(minIndex));
end
medianUsage = (candidates.maxUsage-candidates.minUsage)/2.001; 
%Create mutants iteratively
optStrain  = model;
FChanges   = [];
genesFC    = [];
counter    = 0;
previousFC = 1;
for i=[1 2 3]
    levelCandidates = candidates(candidates.priority==i,:);
    levelCandidates = sortrows(levelCandidates,'foldChange','descend');
    for j=1:length(levelCandidates.genes)
        gene   = levelCandidates.genes{j};
        enzyme = levelCandidates.enzymes{j};
        enzRxn = find(contains(optStrain.rxnNames,['prot_' enzyme]),1);
        if ~isempty(enzRxn)
            iterationModel = optStrain;
            iterationModel = setParam(iterationModel,'obj',enzRxn,1);
            tempSolution   = solveLP(iterationModel);         
            fluxE  = tempSolution.x(enzRxn)>1E-18;
            short  = levelCandidates.shortNames{j};
            action = levelCandidates.actions(j);
            maxUse = levelCandidates.maxUsage(j);
            OEf    = levelCandidates.OE(j);
            %Avoid including enzymes that cannot carry any flux
            if ~(~fluxE & maxUse>0)
                modifications = {gene action OEf};
                if action == 0
                    enzUsage = candidates.maxUsage(i);
                else
                    enzUsage = medianUsage(i);
                end
                tempMutant = getMutant(optStrain,modifications,enzUsage);
                [mutSol,~] = solveECmodel(tempMutant,model,'pFBA',prot_index);
                
                if ~isempty(mutSol)
                    yield = mutSol(targetIndex)/(mutSol(minIndex));
                    FC    = yield/WTyield;
                    %Just keep those genes that don't affect the production phenotype
                    if FC >= (previousFC-1E-12)
                        FChanges   = [FChanges; FC];
                        genesFC    = [genesFC;{gene}];
                        optStrain  = tempMutant;
                        previousFC = FC;
                        counter = counter+1;
                        disp(['Ready with gene #' num2str(counter) ' (' short ')' '  FC:' num2str(FC)])
                    end
                end
            end
        end
    end
end
[~,iB]   = ismember(genesFC,candidates.genes);
iB       = sort(iB,'ascend');
FChanges = table(genesFC,FChanges,'VariableNames',{'genes' 'FC'});
optGenes = candidates(iB,:);
end