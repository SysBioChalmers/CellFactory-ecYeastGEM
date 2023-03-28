clear
current = pwd;
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');
clusters_info  = readtable('../results/cluster_strains/product_clusters.txt','Delimiter','\t');

d    = dir('../results');
isub = [d(:).isdir]; %# returns logical vector
%load model files
load('../ModelFiles/ecYeastGEM_batch.mat');
%initialize variables
prod_capabilities  = table();
rxns               = ecModel_batch.rxns;
fluxes             = table(rxns);
families           = [];
biomass_prod       = true;
protFree           = false;
mainPrec           = false;
Prot_cost          = []; 

%list of 12  main metabolic precursors
if mainPrec
    precursors = {'D-glucose 6-phosphate' 'D-fructose 6-phosphate' ...
                'ribose-5-phosphate' 'D-erythrose 4-phosphate' ...
             	'glyceraldehyde 3-phosphate' '3-phosphonato-D-glycerate(3-)' ...
                'phosphoenolpyruvate' 'pyruvate' ...
                'acetyl-CoA' '2-oxoglutarate' ...
                'succinyl-CoA' 'oxaloacetate'};   
%list of other relevant intracellular molecules to track 
else
    precursors = {'ATP' 'NADH' 'NADPH' 'coenzyme A' 'oxygen'}; 
end
compartments = [1 1 1 1 1 1 1 1 1 1 9 1];
prec_mTO   = table(precursors');
%create directory for generation of yield plots
mkdir('../results/production_capabilities/yieldPlots')
ecModelWT = changeMedia_batch(ecModel_batch,'D-glucose exchange (reversible)','Min',1);
growthPos = find(strcmpi(ecModelWT.rxnNames,'growth'));
CS_index  = find(strcmpi(ecModelWT.rxnNames,'D-glucose exchange (reversible)'));
bio       = calculate_potential(ecModelWT,growthPos,growthPos,CS_index,0.180,true);
metTOsWT  = [];

for k=1:length(precursors)
    midx = find(strcmpi(ecModel_batch.metNames,precursors{k}));
    %midx = midx(find(ecModel_batch.metComps(midx) == compartments(k)));
    metTO = 0;
    for i=1:length(midx)
        
        if ~isempty(midx)
            rxns = find(ecModel_batch.S(midx(i),:));
            [iA,rxns2] = ismember(ecModel_batch.rxns(rxns),bio.fluxDist.rxns);
            rxns = rxns(iA);
            rxns2 = rxns2(find(rxns2));
            %compute turnover numbers
            metTO = metTO + 0.5*sum(abs(ecModel_batch.S(midx(i),rxns))*bio.fluxDist.flux(rxns2));
        else
            disp(precursors{k})
            pause
        end
    end
    metTOsWT = [metTOsWT;metTO];
end


for i=1:height(chemicals_info)
    compound = chemicals_info.Name{i};
    intID    = chemicals_info.internal_ids{i};
    model    = [];
    %try to load GEM
    modelStr = ['../ModelFiles/production_ecModels/' lower(chemicals_info.ecModel{i})];
    try
        modelStr = strrep(modelStr,'.mat','_WBG.mat');
        load(modelStr);
    catch
        modelStr = strrep(modelStr,'_WBG.mat','.mat');
        try
            load(modelStr);
        catch
        end
    end
    
    if ~isempty(model)
        ecModel = model;
        %Check presence of objective reaction
        index_ec = find(model.c);
        indexProt= find(contains(model.rxnNames,'prot_pool_exchange'));
        objRxn   = model.rxns(index_ec);        
        if ~isempty(index_ec)
            %c_source = 'D-glucose exchange';
            if startsWith(compound,'L-')
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
                ecModel = changeMedia_batch(ecModel,'D-glucose exchange (reversible)','YEP',AA);
                %check flux
                sol1    = solveLP(ecModel,1);
                flux_ec = (-sol1.f)>0;
            end
            
            if flux_ec
                disp(['Performing pFBA for: ' compound])
                %Get biomass and product yields for low glucose
                %consumption
                CS_index  = find(strcmpi(ecModel.rxnNames,'D-glucose exchange (reversible)'));
                growthPos = find(strcmpi(ecModel.rxnNames,'growth'));
                %ecModel   = setParam(ecModel,'ub',indexProt,1000);
                ecM = calculate_potential(ecModel,growthPos,index_ec,CS_index,0.180,biomass_prod);
                ecM.fluxDist_1 = ecM.fluxDist;
                enzUsages  = ecM.fluxDist(contains(ecM.fluxDist.rxns,'prot_'),:);
                %Get burden of CCM enzymes
                CCMratio = get_CCM_enzBurden(ecM.fluxDist,ecModel);
                %Get flux distribution
                cost          = ecM.fluxDist.flux(strcmpi(ecM.fluxDist.rxns,ecModel.rxns(indexProt)));
                [presence,iB] = ismember(fluxes.rxns,ecM.fluxDist.rxns);
                fluxes        = fluxes(presence,:);
                ecM.fluxDist  = ecM.fluxDist(iB(presence),:);
                eval(['fluxes.' intID '=ecM.fluxDist.flux;'])    
                %Get metabolic precursors turnover
                metTOs = [];
                for k=1:length(precursors)
                    midx = find(strcmpi(ecModel_batch.metNames,precursors{k}));
                    %midx = midx(find(ecModel_batch.metComps(midx) == compartments(k)));
                    metTO = 0;
                    for j=1:numel(midx)
                    if ~isempty(midx)
                        rxns = find(ecModel_batch.S(midx(j),:));
                        [iA,rxns2] = ismember(ecModel_batch.rxns(rxns),ecM.fluxDist_1.rxns);
                        rxns = rxns(iA);
                        rxns2 = rxns2(find(rxns2));
                        %compute turnover numbers
                        metTO = metTO + 0.5*sum(abs(ecModel_batch.S(midx(j),rxns))*ecM.fluxDist_1.flux(rxns2));                        
                    else
                        disp(precursors{k})
                    end 
                    end
                    metTOs = [metTOs;metTO];
                end
                metTOs = metTOs./metTOsWT;
                eval(['prec_mTO.' intID '=metTOs;'])
                %now with GEM
                %no prot case
                compound = chemicals_info.internal_ids{i};
                UEC = calculate_potential(ecModel,growthPos,index_ec,CS_index,0.180,biomass_prod,protFree);
                newRow            = [{compound},chemicals_info.Group(i),chemicals_info.class(i),chemicals_info.MW(i),UEC.bioYield,UEC.prod_yield,UEC.prod_rate,ecM.bioYield,ecM.prod_yield,ecM.prod_rate,ecM.cFlux_l,ecM.cFlux_h,cost,CCMratio,ecM.prot_yield];
                prod_capabilities = [prod_capabilities; newRow];
                families = [families;chemicals_info.class(i)];
                Prot_cost =[Prot_cost;cost]; 
            end
        end
    else
        sprintf(['ecModel file for: ' compound ' not found.\n'])
    end
end
%
mkdir('../results/production_capabilities')
if mainPrec
    file5 = '../results/met_turnover/met_precursors_turnovers_allChemicals';
else
    file5 = '../results/met_turnover/met_cofactors_turnovers_allChemicals';
end

if biomass_prod
    file1 = '../results/production_capabilities/prodCapabilities_allChemicals_wBio.txt';
    file2 = '../results/production_capabilities/proteinLimitations_allChemicals_wBio.txt';
    file3 = '../results/production_capabilities/fluxDist_distance_allChemicals_wBio.txt';
    file4 = '../results/production_capabilities/fluxDistributions_allChemicals_wBio.txt';
    file5 = [file5 '_wBio.txt'];

elseif biomass_prod==false
    file1 = '../results/production_capabilities/prodCapabilities_allChemicals_noBio.txt';
    file2 = '../results/production_capabilities/proteinLimitations_allChemicals_noBio.txt';
    file3 = '../results/production_capabilities/fluxDist_distance_allChemicals_noBio.txt';
    file4 = '../results/production_capabilities/fluxDistributions_allChemicals_noBio.txt';
    file5 = [file5 '_noBio.txt'];

elseif protFree
    file1 = '../results/production_capabilities/prodCapabilities_allChemicals_noProt.txt';
    file2 = '../results/production_capabilities/proteinLimitations_allChemicals_noProt.txt';
    file3 = '../results/production_capabilities/fluxDist_distance_allChemicals_noProt.txt';
    file4 = '../results/production_capabilities/fluxDistributions_allChemicals_noProt.txt';
    file5 = [file5 '_noProt.txt']

end
    
prod_capabilities.Properties.VariableNames = {'compound' 'type' 'family' 'MW' 'bioYield_gem' 'prodYield_gem' 'prodRate_gem' 'bioYield_ec' 'prodYield_ec' 'prod_rate_ec' 'cFlux_l_ec' 'cFlux_h_ec' 'Pburden' 'CCMratio' 'protScaled'};
writetable(prod_capabilities,file1,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
%calculate euclidean distance matrix (flux distributions)
[m,n]   = size(fluxes);
distMat = zeros(n-1,n-1);
for i=2:n
    for j=2:n
        D = norm(table2array(fluxes(:,i)) - table2array(fluxes(:,j)));
        distMat((i-1),(j-1)) = D;
    end
end
%save distance matrix
distMat = array2table(distMat);
distMat.Properties.VariableNames = fluxes.Properties.VariableNames(2:end);
distMat.Properties.RowNames = fluxes.Properties.VariableNames(2:end);
writetable(fluxes,file4,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
writetable(prec_mTO,file5,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
writetable(distMat,file3,'QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true,'Delimiter','\t')
    