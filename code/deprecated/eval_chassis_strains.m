current = pwd;
chemicals_info    = readtable('../data/chemicals_info.txt','Delimiter','\t');
strain_conditions = readtable('../data/genetic_background.txt','Delimiter','\t');

comp_classes = unique(chemicals_info.class);
class_short  = {'alc' 'alk' 'AAs' 'aro' 'bio' 'FAL' 'fla' 'oAc' 'stb' 'ter'};
%load chassis strain files
files       = dir('../results');
chassis_mod = readtable('../results/chassis_strains_modifications.txt','Delimiter','\t');
chemicals_info.Name = strtrim(chemicals_info.Name);
maxMod = 6;
cd  (current)
fileStr = '../results/chassis_strain_7_chemicals.txt';
results = table();
chemicals = readtable(fileStr,'Delimiter','\t');
chemicals.product = strrep(chemicals.product,'coumaricacid','p-coumaric acid');
chemicals.product = strrep(chemicals.product,'taxadienalphaylacetate', 'taxadien_5alpha_yl_acetate');
chemicals.product = strrep(chemicals.product,'oleanolate','oleate');
chemicals.product = strrep(chemicals.product,'glycyrrhetinicacid','glycyrrhetinic_acid');
chemicals.product = strrep(chemicals.product,'betacarotene','_-carotene');
chemicals.product = strrep(chemicals.product,'betaionone','_-ionone');
chemicals.product = strrep(chemicals.product,'dha','DHA');
chemicals.product = strrep(chemicals.product,'nicotianamine','Nicotianamine');
chemicals.product = strrep(chemicals.product,'sreticuline','(S)-Reticuline');
chemicals.product = strrep(chemicals.product,'epa','EPA');
%iterate through all levels of modification
%load each model
mkdir('../results/chassis_strain')
for j=1:height(chemicals)  
    figure
    for i=1:maxMod
        chem  = chemicals.product(j);
        index = find(contains(chemicals_info.Name,chem),1);
        if ~isempty(index)
            modelStr = chemicals_info.ecModel{index};
            MWeight  = chemicals_info.MW(index);
            model = [];
            try
                %disp(modelStr)
                load(['../ModelFiles/production_ecModels/' modelStr])
            catch
                %disp(modelStr)
                modelStr = strrep(modelStr,'.mat','_WBG.mat');
                load(['../ModelFiles/production_ecModels/' modelStr])
            end
            targets   = chassis_mod(1:i,:);
            genes = [];
            for k = 1:length(targets.geneTargets)
                idx     = find(strcmpi(model.geneShortNames,targets.geneTargets{k}));
                genes   = [genes; model.genes(idx)];
                actions = strcmpi(targets.modifications,'OE')*2;
                actions = num2cell(actions);
            end
            cd method
            modifications = [genes actions actions];
            mutant = getMutant(model,modifications,[],false);
            cd ..
            index = find(contains(strain_conditions.ecModel,modelStr));
            cSource = 'D-glucose exchange (reversible)';
            model   = changeMedia_batch(model,cSource,'Min',false);
            mutant  = changeMedia_batch(mutant,cSource,'Min',false);
            objFlux = haveFlux(model,1E-6,find(model.c));
            objFluxM = haveFlux(mutant,1E-6,find(model.c));
            if ~isempty(index) & ~objFlux & ~objFluxM
                media = strain_conditions.media{index};
                cSource = 'D-glucose exchange (reversible)';
                model  = changeMedia_batch(model,cSource,media,false);
                mutant = changeMedia_batch(mutant,cSource,media,false);           
            end

            MWeight = MWeight/1000;
            %
            objIndx = find(model.c);
            if i==1
                disp(modelStr)
                printModel(model,objIndx)
                [BioYield_WT,yield_WT,oxygen,EtOH] = getYieldPlot(model,objIndx,1,MWeight);
                hold on
                xStr = 'Biomass yield [g_{biomass}/g_{glucose}]';
                yStr = 'Product yield [g_{product}/g_{glucose}]';
                plot(BioYield_WT,yield_WT,'-','LineWidth',5)%,'Color','red')
                xlabel(xStr)
                ylabel(yStr)
                xlim([0 1])
                %ylim([0 1])
                legendStr = {'WT'};
                
            end
            [BioYield_i,yield_i] = getYieldPlot(mutant,objIndx,1,MWeight);
            plot(BioYield_i,yield_i,'-.','LineWidth',2)%,'Color','blue')'
            %legend({'ecModel low' 'ecModel high' 'GEM low' 'GEM high'})
        else
            disp([num2str(i) ': ' chem])
        end
        legendStr = [legendStr;{['Level ' num2str(i)]}];
        title(chem)
        disp(' ')
    end
    legend(legendStr)
    set(gca,'FontSize',22)
    modelStr = strrep(modelStr,'.mat','');
    saveas(gcf,['../results/plots/chassis_strain/'  modelStr '_yieldPlot.jpg'])
    hold off
    close all
    %hold off
    %pause
    %close all
end
results.Properties.VariableNames = {'chemicals' 'mod_number' 'FC_bY' 'FC_pY' 'FC_pR'};
%writetable(results,['../results/chassis_strain_' num2str(maxMod) '_prod_capabilities.txt'],'delimiter','\t','QuoteStrings',false)