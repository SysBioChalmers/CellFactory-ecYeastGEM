function [bioYield,prod_yield, prod_rate,fluxDist,cFlux_l,cFlux_h,prot_yield] = calculate_potential(model,growthPos,targetIndex,CS_index,CS_MW,biomass,prot)
if nargin<7
    prot = false;
if nargin<6
    biomass = true;
end
end
protIndex = find(contains(model.rxns,'prot_pool_exchange'));
%Unconstrain growth
tempModel = setParam(model,'lb',growthPos,0);
tempModel = setParam(tempModel,'ub',growthPos,1000);
%unconstrain production
tempModel = setParam(tempModel,'ub',targetIndex,1000);
tempModel = setParam(tempModel,'lb',targetIndex,0);
%fix unit glucose uptake rate
tempModel = setParam(tempModel,'obj',growthPos,1);
if ~isfield(model,'enzymes')
    tempModel = setParam(tempModel,'ub',CS_index,0);
    tempModel = setParam(tempModel,'lb',CS_index,-1.000001);
else
    tempModel = setParam(tempModel,'ub',CS_index,1.000001);
    tempModel = setParam(tempModel,'lb',CS_index,0);    
end
%Get biomass yield for a unit glucose uptake rate
solution  = solveLP(tempModel);
bioYield  = abs(solution.x(growthPos)/(solution.x(CS_index)*CS_MW));
%disp(['The maximum biomass yield is ' num2str(bioYield) '[g biomass/g carbon source]']);
%fix suboptimal growth
if biomass
    tempModel = setParam(tempModel,'lb',growthPos,0.5*solution.x(growthPos));
end
if prot
    tempModel = setParam(tempModel,'ub',find(contains(model.rxns,'prot_pool_exchange')),1000);
end
%get max. product yield
tempModel  = setParam(tempModel,'obj',targetIndex,1);
solution   = solveLP(tempModel,1);
fluxDist   = table(tempModel.rxns,solution.x,'VariableNames',{'rxns' 'flux'});
prod_yield = abs(solution.x(targetIndex)/(solution.x(CS_index)));
prot_yield = abs(solution.x(protIndex)/(solution.x(CS_index)));
cFlux_l    = solution.x(CS_index);
%disp(['The maximum product yield is ' num2str(prod_yield) '[mol product/mol carbon source]']);
%Get product max. rate for a unconstrained glucose uptake rate
if isfield(model,'enzymes')
    tempModel = setParam(tempModel,'ub',CS_index,1000);
    tempModel = setParam(tempModel,'lb',CS_index,0);
else
    tempModel = setParam(tempModel,'ub',CS_index,0);
    tempModel = setParam(tempModel,'lb',CS_index,-1000);
end
solution  = solveLP(tempModel);
prod_rate = solution.x(targetIndex);
cFlux_h     = solution.x(CS_index);
%disp(['The maximum production rate is ' num2str(prod_rate) '[mmol/gDW h]']);
end