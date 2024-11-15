% Prediction of metabolic engineering targets for increased production of spermidine in Saccharomyces cerevisiae using ecYeastGEM and ecFactory

% 1.Get and load ecYeastGEM (Note: this model is an ecModel for S. cerevisiae. If the user wants to predict targets for other host organisms, they can load other ecModels as needed.)
load('../ModelFiles/ecYeastGEM_batch.mat');
ecModel = ecModel_batch;
current = pwd;

% 2. add GECKO2 and ecFactory code to the path:D:\code\github\ecFactory
ecFactory_path = 'D:\code\github\ecFactory\code/';
addpath(genpath(ecFactory_path));



% 3.- Set case-specific parameters
% Create a structure with the following model-specific parameters:
% modelParam.CSrxn = preferred carbon source uptake reaction ID
% modelParam.CS_MW = Molecular weight of the preferred carbon source
% modelParam.growthRxn = Growth reaction name
% modelParam.rxnTarget = ID for the production target reaction

%Find carbon source uptake reaction 
CSname = 'D-glucose exchange (reversible)';
modelParam.CSrxn  = ecModel.rxns{strcmpi(ecModel.rxnNames,CSname)};
%Indicate carbon source mol. weight in grams/mmol
modelParam.CS_MW  = 0.18015;
%Indicate the name of the growth or biomass pseudoreaction
modelParam.growthRxn = ecModel.rxns{strcmpi(ecModel.rxnNames,'biomass pseudoreaction')};
%identify the production target reaction
product_name         = 'spermidine';
target_pos           = find(strcmpi(ecModel.rxnNames,[product_name ' exchange']));
modelParam.rxnTarget = ecModel.rxns(target_pos);
disp(['*The production reaction for ' product_name ' has been found under the ID: ' modelParam.rxnTarget{1}]);

% Provide a directory path for the results files
%results folder name
% results_folder = 'results/production_targets/spermidine_targets';
results_folder_base = strcat(current,'/../results/production_targets/spermidine_targets');
mkdir(results_folder_base);


% 4. Constrain ecModel
%Set media conditions for batch growth
const_ecModel = changeMedia_batch(ecModel,CSname,'Min'); %model-specific script

CS_index      = find(strcmpi(const_ecModel.rxns,modelParam.CSrxn));
growthPos     = find(strcmpi(const_ecModel.rxns,modelParam.growthRxn));
%Enable cellular growth
const_ecModel = setParam(const_ecModel,'lb',growthPos,0);
const_ecModel = setParam(const_ecModel,'ub',growthPos,1000);
%block artificial growth
revIndex = find(strcmpi(const_ecModel.rxnNames,'growth (reversible)'));
const_ecModel.lb(revIndex) = 0;
const_ecModel.ub(revIndex) = 0;
%Set a fixed unit glucose uptake rate
const_ecModel = setParam(const_ecModel,'ub',CS_index,1);


% 5. Get a feasible biomass yield range to be scanned by the ecFactory method
%Get biomass yield for a unit glucose uptake rate
const_ecModel = setParam(const_ecModel,'obj',growthPos,1);
solution      = solveLP(const_ecModel,1);
WT_yield      = solution.x(growthPos)/(solution.x(CS_index)*modelParam.CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);

%Obtain a suboptimal yield value to run ecFactory
expYield = 0.49*WT_yield;
disp('* The ecFactory method will scan flux distributions spanning from')

disp(['a suboptimal biomass yield of: ' num2str(0.5*expYield) ' to: ' num2str(2*expYield) ' [g biomass/g carbon source]']);


% 7. Run ecFactory method
% try
%     % swich workdir to ecFactory code folder
%     cd(ecFactory_path);
%     [optStrain,candidates,step] = run_ecFactory(const_ecModel,modelParam,expYield,results_folder,true);
% catch
%     disp('The model is not suitable for the ecFactory method')
% end

for i = 1:100
    % Set random seed
    rng(i);
    
    % Create a folder for this iteration's results
    results_folder = strcat(results_folder_base, '/run_', num2str(i));

    % Check if the results folder already exists
    if exist(results_folder, 'dir')
        disp(['Results folder for run ' num2str(i) ' already exists. Skipping this run.']);
        continue;
    end
    
    mkdir(results_folder);
    
    try
        % Switch workdir to ecFactory code folder
        cd(ecFactory_path);
        [optStrain,candidates,step] = run_ecFactory(const_ecModel,modelParam,expYield,results_folder,true);
    catch
        disp(['The model is not suitable for the ecFactory method in run ' num2str(i)]);
    end
end

% return workdir
cd(current);
