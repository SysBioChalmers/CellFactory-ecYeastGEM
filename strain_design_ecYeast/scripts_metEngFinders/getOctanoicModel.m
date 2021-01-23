function strainModel = getOctanoicModel(model,kma,growth,gRate)

glucIndex = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
model.lb(glucIndex) = 0;
model.ub(glucIndex) = 1000;
%Fix minimal growth rate
gIndex = find(strcmpi(model.rxnNames,'growth'));
sol    = solveLP(model);
model.lb(gIndex) = 0.5*sol.x(gIndex);
model.c(:) = 0;

if ~kma
    %Mets IDs
    metsToAdd.mets = {'OctAcid_[p]';'OctAcid_[c]';'OctAcid_[e]'};
    %metsToAdd.mets = {'fatty acid[c]';'fatty acid[e]'};
    %metsToAdd.mets = {'fatty acid[e]'};
    
    %Mets names
    metsToAdd.metNames = {'octanoic acid';'octanoic acid';'octanoic acid'};
    %metsToAdd.metNames = {'fatty acid';'fatty acid'};
    %metsToAdd.metNames = {'fatty acid'};
    
    %Compartments
    metsToAdd.compartments  = {'p';'c';'e'};
    %metsToAdd.compartments  = {'c';'e'};
    %metsToAdd.compartments  = {'e'};
    
    % Reactions
    rxnsToAdd.rxns = {'octanoic acid synthesis';'octanoic acid transport P';...
                      'octanoic acid transport c';'octanoic acid exchange'};
    %rxnsToAdd.rxns = {'fatty acid transport';'fatty acid exchange'};              
    
    %rxnsToAdd.rxnNames  = rxnsToAdd.rxns;
    
    rxnsToAdd.equations = {'H2O[p] + octanoyl-CoA[p] => octanoic acid[p] + coenzyme A[p]';... %synthesis
                          'octanoic acid[p] => octanoic acid[c]';...%Transport
                          'octanoic acid[c] => octanoic acid[e]';...
                          'octanoic acid[e] => '}; %exchange                  
    %rxnsToAdd.equations = {'fatty acid[c] => fatty acid[e]';'fatty acid[e] => '}; %exchange                   
    
    rxnsToAdd.c  = [0 0 0 1];
    %rxnsToAdd.c  = [0 1];
    
    rxnsToAdd.lb = [0 0 0 0];
    %rxnsToAdd.lb = [0 0];
    
    rxnsToAdd.ub = [1000 1000 1000 1000]; 
    %rxnsToAdd.ub = [1000 1000];

else
    metsToAdd.mets          = {'OctAcid_[c]';'OctAcid_[e]'};
    metsToAdd.metNames      = {'octanoic acid';'octanoic acid'};
    metsToAdd.compartments  = {'c';'e'};
    %metsToAdd.unconstrained = [0 0];
    % Reactions
    rxnsToAdd.rxns      = {'octanoic acid synthesis';'octanoic acid transport c';...
                           'octanoic acid exchange'};
    rxnsToAdd.rxnNames  = rxnsToAdd.rxns;
    rxnsToAdd.equations = {'H2O[c] + octanoyl-CoA[c] => octanoic acid[c] + coenzyme A[c]';... %synthesis
                          'octanoic acid[c] => octanoic acid[e]';...
                          'octanoic acid[e] => '}; %exchange
    rxnsToAdd.c  = [0 0 1];
    rxnsToAdd.lb = [0 0 0];
    rxnsToAdd.ub = [1000 1000 1000]; 
end


strainModel = addMets(model,metsToAdd);
strainModel = addRxns(strainModel,rxnsToAdd,3);
%Block competing rxns
%strainModel = blockProduction(strainModel,'peroxisomal acyl-CoA thioesterase','octanoyl-CoA');
%strainModel = blockProduction(strainModel,'Octanoyl-CoA:oxygen 2-oxidoreductase','octanoyl-CoA');
%strainModel = blockProduction(strainModel,'trans-Oct-2-enoyl-CoA reductase','octanoyl-CoA');
%strainModel = blockProduction(strainModel,'fatty-acid--CoA ligase (octanoate) (reversible)','octanoyl-CoA');

end

function model = blockProduction(model,name,product)
Rxns = find(strcmpi(model.rxnNames,name));

for i=1:length(Rxns)
    index = Rxns(i);
    substrates = model.metNames(find(model.S(:,index)<0));
    if any(contains(substrates,product)) 
        disp(name)
        model.ub(index) = 0;
    end
end
end