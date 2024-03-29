current = pwd;
load('../ModelFiles/ecYeastGEM_batch.mat')
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
origin = table();
%retrieve chemical classes info
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');
comp_classes   = unique(chemicals_info.class);
class_short    = {'alc' 'alk' 'AAs' 'aro' 'bio' 'FAL' 'fla' 'oAc' 'stb' 'ter'};
%Iterate through all chemical compounds
d = dir('../results/production_targets');
isub = [d(:).isdir]; %# returns logical vector
nameFolders = {d(isub).name}';
for i=1:length(nameFolders)
    folder = nameFolders{i};
    if contains(folder,'_targets')        
        chemical = strrep(folder,'_targets','');
        %chemical = regexprep(chemical,'[^a-zA-Z]','');
        str = strrep(folder,'_targets','');
        model_idx = find(strcmpi(chemicals_info.internal_ids,str));
        if ~isempty(model_idx)
            class = find(strcmpi(comp_classes,chemicals_info.class(model_idx)));
            native = chemicals_info.Group(model_idx);
        %else
        %    disp(['No class: ' chemical])
        %end
        idx    = find(strcmpi(comp_classes,class));
        %if ~isempty(idx)
        newStr = [lower(chemical) '_fam_' class_short{class}];
        %newStr = lower(chemical);
        newStr = strrep(newStr,' ','_');
        newStr = strrep(newStr,'-','_');
        newStr = strrep(newStr,',','');
        newStr = strrep(newStr,'(','');
        newStr = strrep(newStr,')','');
        newStr = regexprep(newStr,'[0-9]','');
        newStr = strrep(newStr,'__','');
        while startsWith(newStr,'_')
        	newStr = newStr(2:end);
        end
        
        %disp(newStr)
        %create new column for chemical
        eval(['genesTable.' newStr '=ones(height(genesTable),1);'])
        %Open targets file
      %  try
            candidates = readtable(['../results/production_targets/' folder '/candidates_L3.txt'],'Delimiter','\t');
                        
            OEs=candidates.genes(candidates.k_scores>1);
            [~,iB]=ismember(OEs,genesTable.genes);
            eval(['genesTable.' newStr '(iB(iB>0)) = 4;'])
            dels=candidates.genes(candidates.k_scores<=0.05);
            [~,iA]=ismember(dels,genesTable.genes);
            eval(['genesTable.' newStr '(iA(iA>0)) = 0;'])
            dRegs=candidates.genes(candidates.k_scores>0.05 & candidates.k_scores<=0.5);
            [~,iA]=ismember(dRegs,genesTable.genes);
            eval(['genesTable.' newStr '(iA(iA>0)) = 0.25;'])
            
%             [iA,iB]=ismember(candidates.genes,genesTable.genes);
% %            values = candidates.maxUsage(find(iA))./(1E-15+candidates.maxUsageBio(find(iA)));
%             %values(values<=0.05) = 0;
%             %eval(['genesTable.' newStr '(iB(iB>0)) = candidates.k_scores(find(iA));'])
%             eval(['genesTable.' newStr '(iB(iB>0)) = values;'])
            origin = [origin;[chemical {class} native]];
%         catch
%            disp(chemical)
%         end
        end
    end
end
%writetable(genesTable,'../results/production_targets/targetsMatrix_mech_validated.txt','delimiter','\t','QuoteStrings',false)
%writetable(genesTable,'../results/production_targets/targetsMatrix_compatible.txt','delimiter','\t','QuoteStrings',false)
writetable(genesTable,'../results/production_targets/targetsMatrix_discrete_EV.txt','delimiter','\t','QuoteStrings',false)
writetable(origin,'../results/production_targets/chemicals_origin.txt','delimiter','\t','QuoteStrings',false)
