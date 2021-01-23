%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = changeModel(model)

%Correct glucan coefficients in biomass reaction:
model.S(strcmp(model.mets,'s_0002'),strcmp(model.rxns,'r_4041')) = 0;
model.S(strcmp(model.mets,'s_0001'),strcmp(model.rxns,'r_4041')) = -0.8506;
model.S(strcmp(model.mets,'s_0004'),strcmp(model.rxns,'r_4041')) = -0.2842;

%Correctly represent proton balance inside cell:
model.lb(strcmp(model.rxns,'r_1824')) = 0;  %Block free H+ export
model.ub(strcmp(model.rxns,'r_1250')) = 0;  %Block free putrescine export
model.ub(strcmp(model.rxns,'r_1259')) = 0;  %Block free spermidine export

%COMPLEX V: For 1 ATP 3 H+ are needed, not 4:
model.S(strcmp(model.mets,'s_0799'),strcmp(model.rxns,'r_0226')) = +2;
model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_0226')) = -3;

%Delete blocked rxns (LB = UB = 0):
to_remove = boolean((model.lb == 0).*(model.ub == 0));
model     = removeRxns(model,model.rxns(to_remove));

%Correct rev vector: true if LB < 0 & UB > 0 (and it's not an exchange reaction):
model.rev = false(size(model.rxns));
for i = 1:length(model.rxns)
    if (model.lb(i) < 0 && model.ub(i) > 0) || ~isempty(strfind(model.rxnNames{i},'exchange'))
        model.rev(i) = true;
    end
end

%Additional changes:
model = standardizeModel(model,'COBRA');
model = convertToIrreversibleModel(model);
model = fixedModifications(model);

%Manual changes: block weird rxns:
model.ub(strcmp(model.rxns,'r_2045_REV')) = 0; %L-serine transport from mit to cyt
model.ub(strcmp(model.rxns,'r_0659'))     = 0; %isocit + NADP(+) -> 2-oxoglut + CO2 + NADPH 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%