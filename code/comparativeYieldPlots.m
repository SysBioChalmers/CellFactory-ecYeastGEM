function comparativeYieldPlots
clear
current = pwd;
%subSystems_GEM = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');
d         = dir('../results');
isub      = [d(:).isdir]; %# returns logical vector
load('../ModelFiles/ecYeastGEM_batch.mat')
load('../ModelFiles/yeastGEM.mat')
GEM  = ravenCobraWrapper(model);
chemicals_info.Name = strtrim(chemicals_info.Name);
maxRate = table([],[],[],[]);
failed = [];
%Loop through all maxRate folders
mkdir('../results/production_capabilities/yieldPlots')
for i=1:height(chemicals_info)
    compound = chemicals_info.Name{i};
    model    = [];
    %try to load GEM
    modelStr = ['../ModelFiles/production_ecModels/' chemicals_info.ecModel{i}];
    try
        load(modelStr);
    catch
        modelStr = strrep(modelStr,'.mat','_WBG.mat');
        load(modelStr);
    end
    
    if ~isempty(model)
        ecModel = model;
        MWeight = chemicals_info.MW(i)/1000;
        %Check presence of objective reaction
        index_ec = find(model.c);
        objRxn   = model.rxns(index_ec);
        indexGEM = find(strcmpi(GEM.rxns,objRxn));
        if ~isempty(index_ec)
            %c_source = 'D-glucose exchange';
            if startsWith(compound,'L-') | contains(compound,'lutathione') 
                AA = true;
            else
                AA = false;
            end
            ecModel = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','Min',AA);
            %unconstrain growth
            ecModel = setParam(ecModel,'lb',find(strcmpi(ecModel.rxnNames,'growth')),0);
            %Check growth capabilities
            temp = setParam(ecModel,'obj',find(strcmpi(ecModel.rxnNames,'growth')),1);
            sol1 = solveLP(temp);
            gRate = -sol1.f>0;
            %check flux
            sol1    = solveLP(ecModel);
            obj1    = -sol1.f;
            flux_ec = obj1>=1E-6;

            if ~flux_ec | ~gRate
                pause
                ecModel = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','YEP',AA);
                %check flux
                sol1    = solveLP(ecModel,1);
                flux_ec = (-sol1.f)>0;
            end
            
            if flux_ec
                disp(['Generating yield plots for: ' compound])
                %Get biomass and product yields for low glucose
                %consumption
                temp = setParam(ecModel,'obj','r_2111',1);
                [BioYield_ec_l,yield_ec_l] = getYieldPlot(ecModel,index_ec,1,MWeight);
                [BioYield_ec_h,yield_ec_h] = getYieldPlot(ecModel,index_ec,10,MWeight);
                if any(yield_ec_l>1) | any(yield_ec_h>1)
                    disp(['Yield > 1 for : ' compound ' (ecModel)'])
                end
                %Plot results for ecModel
                ylimit = 1;
                if max(yield_ec_l)<= 0.1
                    ylimit = 0.1;
                    if max(yield_ec_l)<= 0.01
                        ylimit = 0.01;
                        if max(yield_ec_l)<= 0.001
                            ylimit = 0.001;
                            if max(yield_ec_l)<= 0.0001
                                ylimit = 0.0001;
                                if max(yield_ec_l)<= 0.00001
                                    ylimit = 0.00001;
                                end
                            end
                        end
                    end
                end  
                figure
                hold on
                xStr = 'Biomass yield [g_{biomass}/g_{glucose}]';
                yStr = 'Product yield [g_{product}/g_{glucose}]';
                plot(BioYield_ec_l,yield_ec_l,'-','LineWidth',4)%,'Color','red')
                plot(BioYield_ec_h,yield_ec_h,'-','LineWidth',4)%,'Color','purple')
                xlabel(xStr)
                ylabel(yStr)
                xlim([0 1])
                ylim([0 ylimit])
                %If objective reaction is also present in GEM, then add GEM
                %yields to plot
                if ~isempty(indexGEM)
                    GEM = changeMedia_Original(GEM,AA,1,1000);
                    GEM = setParam(GEM,'obj',indexGEM,1);
                    sol2  = solveLP(GEM);
                    obj2  = -sol2.f;
                    FC    = obj2/obj1;
                    if obj2>0
                        maxRate = [maxRate;{compound} {-sol1.f} {-sol2.f} {FC}];
                        [BioYield_l,yield_l] = getYieldPlot(GEM,indexGEM,1,MWeight);
                        [BioYield_h,yield_h] = getYieldPlot(GEM,indexGEM,10,MWeight);
                        plot(BioYield_l,yield_l,'-.','LineWidth',4)%,'Color','blue')'
                        plot(BioYield_h,yield_h,'-.','LineWidth',4)%,'Color','yellow')
                        legend({'ecModel low' 'ecModel high' 'GEM low' 'GEM high'})
                        if any(yield_l>1) | any(yield_h>1)
                            disp(['Yield > 1 for : ' compound ' (GEM)'])
                        end
                        %redefine plot limits (if needed)
                        if ylimit<max(yield_l)
                            ylimit = ceil(log10(max(yield_l)));
                            ylimit = 10^(ylimit);
                            ylim([0 ylimit])
                        end
                    end
                else
                    legend({'ecModel low' 'ecModel high'})
                end
                %save plot file a jpg
                set(gca,'FontSize',22)
                saveas(gcf,['../results/production_capabilities/yieldPlots/'  compound '_yieldPlot.jpg'])
                hold off
                close all
            end
        end
    else
        sprintf(['GEM file for: ' compound ' not found.\n'])
        failed = [failed; {compound}];
    end
end
end