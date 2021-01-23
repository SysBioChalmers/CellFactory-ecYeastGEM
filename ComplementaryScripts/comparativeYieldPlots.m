clear
current = pwd;
%subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
d         = dir('../results');
isub      = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
load('../ModelFiles/mat/yeastGEM.mat')
model = ravenCobraWrapper(model);
subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
models = [];
chemicals_info.Name = strtrim(chemicals_info.Name);
maxYield = table([],[],[],[]);
%Loop through all maxYield folders
mkdir(['../results/yieldPlots'])
for i=1:length(nameFolds)
    folder = nameFolds{i};
    cd (current)
    if contains(folder,'_targets')
        compound  = strrep(folder,'_targets','');
        models    = [models; {['ec' compound]}];
        model_idx = find(strcmpi(chemicals_info.ecModel,['ec' strrep(folder,'_targets','') '.mat']));
        if ~isempty(model_idx)
            if strcmpi(chemicals_info.Group{model_idx},'native')
                native = true;
            else
                native = false;
            end
            met     = chemicals_info.Name{model_idx};
            MWeight = chemicals_info.MW(model_idx)/1000;
            index   = find(strcmpi(ecModel_batch.metNames,met));
            %Check presence of exchange reaction
            str      = [met ' exchange'];
            index_ec = find(strcmpi(ecModel_batch.rxnNames,str));
            index    = find(strcmpi(model.rxnNames,str));
            if ~isempty(index_ec)
                c_source = 'D-glucose exchange';
                ecModel  = setParam(ecModel_batch,'obj',index_ec,1);
                ecModel  = changeMedia_batch(ecModel,[c_source ' (reversible)'],'Min',false,1);
                sol1 = solveLP(ecModel);
                obj1 = -sol1.f;
                %Get biomass and product yields for low glucose
                %consumption
                [BioYield_ec_l,yield_ec_l] = getYieldPlot(ecModel,index_ec,1,MWeight);
                [BioYield_ec_h,yield_ec_h] = getYieldPlot(ecModel,index_ec,100,MWeight);
                %Plot results for ecModel
                figure
                hold on
                xStr = 'Biomass yield [g_{biomass}/g_{glucose}]';
                yStr = 'Product yield [g_{product}/g_{glucose}]';
                plot(BioYield_ec_l,yield_ec_l,'-.','LineWidth',5)%,'Color','red')
                plot(BioYield_ec_h,yield_ec_h,'-.','LineWidth',5)%,'Color','purple')
                xlabel(xStr)
                ylabel(yStr)
                xlim([0 1])
                ylim([0 1])
                    
                
                index    = find(strcmpi(model.rxnNames,str));
                if ~isempty(index)
                    model = setParam(model,'obj',index,1);
                    model = changeMedia_Original(model,1,1,1);
                    sol2  = solveLP(model);
                    obj2  = -sol2.f;
                    FC    = obj2/obj1;
                    maxYield = [maxYield;{met} {-sol1.f} {-sol2.f} {FC}];
                    [BioYield_l,yield_l] = getYieldPlot(model,index,1,MWeight);
                    [BioYield_h,yield_h] = getYieldPlot(model,index,20,MWeight);
                    plot(BioYield_l,yield_l,'-','LineWidth',5)%,'Color','blue')'
                    plot(BioYield_h,yield_h,'-','LineWidth',5)%,'Color','yellow')
                    legend({'ecModel low' 'ecModel high' 'GEM low' 'GEM high'})
                else
                    legend({'ecModel low' 'ecModel high'})
                end
                set(gca,'FontSize',22)
                saveas(gcf,['../results/yieldPlots/'  met '_yieldPlot.jpg'])
                hold off
                close all 
            end
        end
    else
        disp(folder)
        class = {''};
    end
end