cand_1_total = [];
cand_1_OE    = [];
cand_1_del   = [];
cand_1_dR    = [];

cand_2_total = [];
cand_2_OE    = [];
cand_2_del   = [];
cand_2_dR    = [];

cand_3_total = [];
cand_3_OE    = [];
cand_3_del   = [];
cand_3_dR    = [];

del_targets = [];
OE_targets  = [];
dR_targets  = [];

del_genes = [];
OE_genes  = [];
dR_genes  = [];

subSystems_del = [];
subSystems_OE  = [];
subSystems_dR  = [];

chem_class_del = [];
chem_class_OE  = [];
chem_class_dR  = [];

chem_comp_del = [];
chem_comp_OE  = [];
chem_comp_dR  = [];

models = [];
chemClass = [];
current = pwd;
%fileNames  = dir();
d = dir('../results');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
load('../ModelFiles/mat/ecYeastGEM_batch.mat')
subSystems_model = mapEnzymeSubSystems(ecModel_batch.enzymes,ecModel_batch);
chemicals_info = readtable('../ComplementaryData/chemicals_info.txt','Delimiter','\t');
for i=1:length(nameFolds)
    folder = nameFolds{i};
    cd (current)
    if contains(folder,'_targets')
        compound  = strrep(folder,'_targets','');
        models    = [models; {['ec' compound]}];
        model_idx = find(strcmpi(chemicals_info.ecModel,['ec' strrep(folder,'_targets','') '.mat']));
        if ~isempty(model_idx)
            class = chemicals_info.class(model_idx);
        else
            class = {''};
        end
        chemClass = [chemClass; class];
        try
            cand_1       = readtable(['../results/' folder '/candidates_ecFSEOF.txt'],'Delimiter','\t');
            cand_1_del   = [cand_1_del;sum(cand_1.k_scores<=0.05)];
            cand_1_dR    = [cand_1_dR;sum(cand_1.k_scores>0.05 & cand_1.k_scores<=0.5)];
            cand_1_OE    = [cand_1_OE;sum(cand_1.k_scores>=2)];
            cand_1_total = [cand_1_total;height(cand_1)];

        catch
            cand_1_total = [cand_1_total;0];
            cand_1_del   = [cand_1_del;0];
            cand_1_OE    = [cand_1_OE;0];
            cand_1_dR    = [cand_1_dR;0];
        end
        
        try
            cand_2 = readtable(['../results/' folder '/candidates_mech_validated.txt'],'Delimiter','\t');
            cand_2_total = [cand_2_total;height(cand_2)];
            cand_2_del   = [cand_2_del;sum(cand_2.k_scores<=0.05)];
            cand_2_dR    = [cand_2_dR;sum(cand_2.k_scores>0.05 & cand_2.k_scores<=0.5)];
            cand_2_OE    = [cand_2_OE;sum(cand_2.k_scores>=2)];
        catch
            cand_2_total = [cand_2_total;0];
            cand_2_dR    = [cand_2_dR;0];
            cand_2_del   = [cand_2_del;0];
            cand_2_OE    = [cand_2_OE;0];
        end
              
        try
            cand_3 = readtable(['../results/' folder '/compatible_genes_results.txt'],'Delimiter','\t');
            %cand_3 = cand_3(cand_3.k_scores>=2 | cand_3.k_scores<=0.05,:);
            %cand_3 = cand_3(cand_3.priority==1,:);
            
            cand_3_total   = [cand_3_total;height(cand_3)];
            idxs           = find(contains(cand_3.actions,'deletion') & cand_3.k_scores<=0.05);
            cand_3_del     = [cand_3_del;length(idxs)];
            deletions      = cand_3.enzymes(idxs);
            subSystems_del = [subSystems_del;mapEnzymeSubSystems(deletions,ecModel_batch)];
            del_genes      = [del_genes;cand_3.genes(idxs)];
            del_targets    = [del_targets;cand_3.shortNames(idxs)];
            temp           = repelem(class,length(idxs),1);
            tempComp       = repelem({compound},length(idxs),1);
            chem_class_del = [chem_class_del;temp];
            chem_comp_del  = [chem_comp_del;tempComp];
            
            idxs           = find(cand_3.k_scores>0.05 & cand_3.k_scores<=0.5);
            cand_3_dR      = [cand_3_dR;length(idxs)];
            dRegs          = cand_3.enzymes(idxs);
            subSystems_dR  = [subSystems_dR;mapEnzymeSubSystems(dRegs,ecModel_batch)];
            dR_genes       = [dR_genes;cand_3.genes(idxs)];
            dR_targets     = [dR_targets;cand_3.shortNames(idxs)];
            temp           = repelem(class,length(idxs),1);
            tempComp       = repelem({compound},length(idxs),1);
            chem_class_dR  = [chem_class_dR;temp];
            chem_comp_dR   = [chem_comp_dR;tempComp];
            
            cand_3_OE     = [cand_3_OE;sum(contains(cand_3.actions,'OE'))];
            idxs          = find(contains(cand_3.actions,'OE') & cand_3.k_scores>=2);
            OE_genes      = [OE_genes;cand_3.genes(idxs)];
            OE_targets    = [OE_targets;cand_3.shortNames(idxs)];
            temp          = repelem(class,length(idxs),1);
            tempComp      = repelem({compound},length(idxs),1);
            chem_class_OE = [chem_class_OE;temp];
            OEs           = cand_3.enzymes(idxs);
            subSystems_OE = [subSystems_OE;mapEnzymeSubSystems(OEs,ecModel_batch)];
            chem_comp_OE  = [chem_comp_OE;tempComp];

        catch
            cand_3_total = [cand_3_total;0];
            cand_3_del   = [cand_3_del;0];
            cand_3_dR   = [cand_3_dR;0];
            cand_3_OE    = [cand_3_OE;0];
        end    
    end
    
end
models = strrep(models,'_targets','');
t = table(models,chemClass,cand_1_total,cand_1_OE,cand_1_dR,cand_1_del,cand_2_total,cand_2_OE,cand_2_dR,cand_2_del,cand_3_total,cand_3_OE,cand_3_dR,cand_3_del);
writetable(t,'../results/targets_summary.txt','delimiter','\t');

t = table(del_genes,del_targets,subSystems_del,chem_class_del,chem_comp_del);
writetable(t,'../results/all_deletions.txt','delimiter','\t','QuoteStrings',false);

t = table(dR_genes,dR_targets,subSystems_dR,chem_class_dR,chem_comp_dR);
writetable(t,'../results/all_downRegs.txt','delimiter','\t','QuoteStrings',false);

expression = '(\w+)''';
OE_targets = regexprep(OE_targets,expression,'_');
OE_targets = strrep(OE_targets,'4"OMT','4_OMT');
OE_targets = strrep(OE_targets,'F3"H','F3_H');
OE_genes = strrep(OE_genes,'4"OMT','4_OMT');
OE_genes = strrep(OE_genes,'F3"H','F3_H');

t = table(OE_genes,OE_targets,subSystems_OE,chem_class_OE,chem_comp_OE);
writetable(t,'../results/all_OEs.txt','delimiter','\t','QuoteStrings',false);
