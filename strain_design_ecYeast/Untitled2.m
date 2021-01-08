Kcat1=115*3600;
MW1=27.249;
model = addReaction(model,'3HP_newRxn1','metaboliteList',{'s_0794','s_1212','s_4184','prot_P39831','s_1207','s_4207'},'stoichCoeffList',[-1 -1 -1 -1/Kcat1 1 1],'reversible',false);
model = addReaction(model,'3HP_newRxn7','metaboliteList',{'prot_pool','prot_P39831'},'stoichCoeffList',[-MW1 1],'reversible',false);
% Add rxns of gene bcere0029_32090
model = addReaction(model,'3HP_newRxn4','metaboliteList',{'s_0441','s_1399','s_0955','s_4184'},'stoichCoeffList',[-1 -1 1 1],'reversible',false);
% Add rxns of gene A7U8C7
Kcat2=7.03*3600;
MW2=61.24;
model = addReaction(model,'3HP_newRxn5','metaboliteList',{'s_0794','s_0973','prot_A7U8C7','s_0441','s_0456'},'stoichCoeffList',[-1 -1 -1/Kcat2 1 1],'reversible',false);
model = addReaction(model,'3HP_newRxn8','metaboliteList',{'prot_pool','prot_A7U8C7'},'stoichCoeffList',[-MW2 1],'reversible',false);
% Add rxns of gene YGR019Wly
Kcat3=0.1324*3600;
MW3=52.946;
model = addReaction(model,'3HP_newRxn6','metaboliteList',{'s_0180','s_0441','prot_P17649ly','s_0991','s_4184'},'stoichCoeffList',[-1 -1 -1/Kcat3 1 1],'reversible',false);
model = addReaction(model,'3HP_newRxn9','metaboliteList',{'prot_pool','prot_P17649ly'},'stoichCoeffList',[-MW3 1],'reversible',false);
% Add transport and exchange rxns of 3HP
model = addReaction(model,'3HP_newRxn2','metaboliteList',{'s_4207','s_4208'},'stoichCoeffList',[-1 1],'reversible',false);
model = addReaction(model,'3HP_newRxn3','metaboliteList',{'s_4208'},'stoichCoeffList',[-1],'reversible',false);

% Add gene rules to the reaction
model=changeGeneAssociation(model,'3HP_newRxn1','fdyG');
model=changeGeneAssociation(model,'3HP_newRxn4','bcere0029_32090');
model=changeGeneAssociation(model,'3HP_newRxn5','A7U8C7');
model=changeGeneAssociation(model,'3HP_newRxn6','YGR019Wly');

% Normalization of geneShortNames, metComps, enzymes, and enzGenes
model.geneShortNames(1128)={'fdyG'};
model.geneShortNames(1129)={'bcere0029_32090'};
model.geneShortNames(1130)={'A7U8C7'};
model.geneShortNames(1131)={'YGR019Wly'};

model.metComps(4147)=1;
model.metComps(4148)=1;
model.metComps(4149)=1;
model.metComps(4150)=1;
model.metComps(4151)=3;

model.enzymes(964)={'P39831'};
model.enzymes(965)={'A7U8C7'};
model.enzymes(966)={'P17649ly'};

model.enzGenes(964)={'fdyG'};
model.enzGenes(965)={'A7U8C7'};
model.enzGenes(966)={'YGR019Wly'};
