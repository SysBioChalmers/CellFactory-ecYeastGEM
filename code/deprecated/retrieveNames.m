clear
fileNames  = dir('../ModelFiles/production_ecModels');
chemTable  = readtable('../data/chemicals_info.txt','Delimiter','\t');
targetsMat = readtable('../results/processed_results/targetsMatrix_mech_validated.txt','Delimiter','\t');
chemNames  = targetsMat.Properties.VariableNames;
conditions = readtable('../data/strain_conditions.txt','Delimiter','\t');


clc
internalIds = cell(height(chemTable),1);
for i=1:length(fileNames)
    file = fileNames(i).name;
    if contains(file,'.mat') & startsWith(file,'ec')
        modelName = file;%file(1:end-4);
        modelName = strrep(modelName,'_WBG','');
        %modelName = lower(modelName);
        chemical = lower(modelName(3:(end-4)));
        index = find(strcmpi(conditions.Chemicals,chemical));
        if isempty(index)
            disp(['No hay conditions: ' chemical])
        end
        chemical = strrep(chemical,' ','_');
        chemical = strrep(chemical,'-','_');
        chemical = strrep(chemical,',','');
        chemical = strrep(chemical,'(','');
        chemical = strrep(chemical,')','');
        chemical = regexprep(chemical,'[0-9]','');
        chemical = strrep(chemical,'__','');
        while startsWith(chemical,'_')
            chemical = chemical(2:end);
        end
        %chemical = regexprep(chemical,'[0-9]','');
        %disp([modelName ':  #' num2str(i)])
        index = find(strcmp(chemTable.ecModel,modelName));
        if isempty(index)
            disp(['No hay model: ' chemical])
        else
            internalIds{index} = chemical;
            if ~isfolder(['../results/production_targets/' chemical '_targets'])
                disp(['No hay targets: ' chemical])
            end
            
            index = find(strcmp(chemNames,chemical));
            if isempty(index)
                disp(['No hay matrix: ' chemical])
            end
        end
        
    end
end
chemTable.internal_ids = internalIds;
%writetable(chemTable,'../data/chemicals_info.txt','Delimiter','\t','QuoteStrings',false);
