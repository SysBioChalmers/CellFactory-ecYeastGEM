current = pwd;
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
[presence,iB] = ismember(ecModel_batch.genes,ecModel_batch.enzGenes);
genesTable = table(ecModel_batch.genes,ecModel_batch.geneShortNames,'VariableNames',{'genes' 'shortNames'});
genesTable.enzymes = cell(height(genesTable),1);
genesTable.subSystems = cell(height(genesTable),1);
subSystems = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
%Avoid empty subsystems
newCol = strrep(subSystems,' // ','');
emptyCells = cellfun(@isempty,newCol);
subSystems(emptyCells) = {''};
genesTable.enzymes(presence) = ecModel_batch.enzymes(iB(iB~=0));
genesTable.subSystems(presence) = subSystems(iB(iB~=0));
%retrieve chemical classes info
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
comp_classes   = unique(chemicals_info.class);
%Iterate through all chemical compounds
d = dir('../results');
isub = [d(:).isdir]; %# returns logical vector
nameFolders = {d(isub).name}';
for i=1:length(nameFolders)
    folder = nameFolders{i};
    if contains(folder,'_targets')        
        chemical = strrep(folder,'_targets','');
        chemical = regexprep(chemical,'[^a-zA-Z]','');
        model_idx = find(strcmpi(chemicals_info.ecModel,['ec' strrep(folder,'_targets','') '.mat']));
        if ~isempty(model_idx)
            class = find(strcmpi(comp_classes,chemicals_info.class(model_idx)));
        else
            class = 0;
        end
        idx    = find(strcmpi(comp_classes,class));
        newStr = [chemical '_del_' comp_classes{class}];
        newStr = strrep(newStr,' ','_');
        disp(newStr)
        %create new column for chemical
        eval(['genesTable.' newStr '=zeros(height(genesTable),1);'])
        %Open targets file
        try
            candidates = readtable(['../results/' folder '/candidates_mech_validated.txt'],'Delimiter','\t');
            OEs =candidates.genes(strcmpi(candidates.actions,'OE'));
            OEs =candidates.genes(candidates.actions>0);

            [~,iB]=ismember(OEs,genesTable.genes);
            eval(['genesTable.' newStr '(iB(iB>0)) = 3;'])
            dels=candidates.genes(candidates.k_scores<=0.05);
            [~,iA]=ismember(dels,genesTable.genes);
            eval(['genesTable.' newStr '(iA(iA>0)) = 1;'])
            dRegs=candidates.genes(candidates.k_scores>0.05 & candidates.k_scores<=0.5);
            [~,iA]=ismember(dRegs,genesTable.genes);
            eval(['genesTable.' newStr '(iA(iA>0)) = 2;'])
        catch
            %disp(chemical)
        end
    end
end
writetable(genesTable,'../results/targetsMatrix_mech_validated.txt','delimiter','\t','QuoteStrings',false)
