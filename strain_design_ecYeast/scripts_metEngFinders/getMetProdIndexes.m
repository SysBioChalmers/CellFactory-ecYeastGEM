function [metindx,rxnIndexes,stoich] = getMetProdIndexes(model,metName,direction,compName)
% getMetProdIndexes
%
% Function that get the indexes for all the reactions that produce or
% consume a given metabolite in the indicated compartment.
% 
% model         a GEM matlab structure
% metName       metabolite name
% direction     'production' if all the reactions that have the metabolite
%               as a product are requested, 'consumption' if the reactions 
%               that consume the metabolite are required instead
% compName      Name of the metabolite in which the reactions for the 
%               metabolite should be located
%
% metindx       Index of the metabolite in the model (accounting for
%               compartmentalization)
% rxnIndexes    Production/consumption rxns indexes
% stoich        Stoichiometric coefficients of the metabolite in each of
%               the found reactions.
% 
% Created.  Ivan Domenzain 2018-10-15
%
rxnIndexes  = [];
stoich      = [];
metindx     = find(strcmpi(model.metNames,metName));

if strcmpi(direction,'production')
    coeff = 1;
elseif strcmpi(direction,'consumption')
    coeff = -1;
end

if nargin<4
    for i=1:length(metindx)
        row = model.S(metindx(i),:);
        rxnIndexes  = [rxnIndexes,find(coeff*row>0)];
    end
else
    compIndex   = find(strcmpi(model.compNames,compName));
    metindx     = metindx(find(model.metComps(metindx)==compIndex));
    row         = model.S(metindx,:);
    rxnIndexes  = find(coeff*row>0);
    %Save the stoichiometric coefficient for the metabolite in each of its
    %reactions for calculating total production rates
    stoich = full(model.S(metindx,rxnIndexes))';
end
end