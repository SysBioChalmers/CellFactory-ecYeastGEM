function cost_trajectory
clear
current = pwd;
%subSystems_GEM = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../data/chemicals_info.txt','Delimiter','\t');

chemicals_info.Name = strtrim(chemicals_info.Name);
prod_capabilities   = table();
biomass_prod = false;

%Loop through all maxRate folders
mkdir('../results/production_capabilities/yieldPlots')
products = [{'Psilocybin'};{'Miltiradiene'};{'Valencene'}];
MW = [0.28425, 0.2725, 0.20435];

for i = 1:length(products)
    
    prod = products{i};
    modelStr  = ['../ModelFiles/production_ecModels/ec' prod '_WBG.mat'];
    prod = lower(prod);
    load('../ModelFiles/ecYeastGEM_batch.mat');
    ecModel = ecModel_batch;
    load(modelStr);
    candidates    = readtable(['../results/production_targets/' prod '_targets/compatible_genes_results.txt'],'Delimiter','\t');
    modifications = [candidates.genes(1),{2},{2}];

    OEfactors     = [1, 2, 5,10,20,50,100];
    mutantModel   = getMutant(model,modifications);
    
    index_mut     = find(model.c);
    indexProt_mut = find(contains(model.rxnNames,'prot_pool_exchange'));
    objRxn_mut    = model.rxns(index_mut);
    %If objective reaction is also present in GEM, then add GEM
    CS_index_mut  = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
    ox_index_mut  = find(strcmpi(model.rxnNames,'oxygen exchange (reversible)'));

    growthPos_mut = find(strcmpi(model.rxnNames,'growth'));
    if i==1
        figure
    end
    
    for j=1:length(OEfactors)
        OEf = OEfactors(j);
        modifications = [candidates.genes(1),{2},{OEfactors(j)}];
        mutantModel   = getMutant(model,modifications);
        %unconstrain growth
        mutantModel = setParam(mutantModel,'lb',find(strcmpi(mutantModel.rxnNames,'growth')),0);
        %Check growth capabilities
        mutantModel = setParam(mutantModel,'obj',find(strcmpi(mutantModel.rxnNames,'growth')),1);
        mutantModel = changeMedia_batch(mutantModel,'D-glucose exchange (reversible)','Min',false);
        [bioY_mut,proY_mut,rate_mut,fluxDist_mut,cFlux_l_mut,cFlux_h_mut] = calculate_potential(mutantModel,growthPos_mut,index_mut,CS_index_mut,0.180,biomass_prod);
        %Get flux distribution
        Pcost_mut = fluxDist_mut.flux(strcmpi(fluxDist_mut.rxns,mutantModel.rxns(indexProt_mut)))/(rate_mut*MW(i));
        Ox_in_mut = fluxDist_mut.flux(ox_index_mut);
        Ccost_mut = fluxDist_mut.flux(strcmpi(fluxDist_mut.rxns,mutantModel.rxns(CS_index_mut)))*0.18/(rate_mut*MW(i));
        newRow = [{prod}, {OEf}, {Pcost_mut},{Ccost_mut},{fluxDist_mut.flux(CS_index_mut)} {rate_mut} {Ox_in_mut}];
        prod_capabilities = [prod_capabilities; newRow];
        if i==1
            [BioYield_ec_l,yield_ec_l] = getYieldPlot(mutantModel,index_mut,1,MW(i));
            plot(BioYield_ec_l,yield_ec_l,'LineWidth',4)%,'Color','blue')'
            hold on
        end
        
    end
    
    if i==1
        legend({'1' '2' '5' '10' '20' '50' '100'})
        %save plot file a jpg
        set(gca,'FontSize',22)
        saveas(gcf,['../results/production_capabilities/plots/psilocybin_mutantCost.jpg'])
        hold off
        %close all
    end

end
prod_capabilities.Properties.VariableNames = {'product' 'OEf' 'Pcost' 'Ccost' 'c_uptake' 'p_rate' 'O2_in'};
writetable(prod_capabilities,'../results/production_capabilities/mutant_costs.txt','Delimiter','\t','QuoteStrings',false)
end