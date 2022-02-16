%% simulationCrabtree

% Timing: ~ 700 s
tic;
mkdir('crabtree')
cd('crabtree')
% load model
load('pcSecYeast.mat');
load('enzymedata.mat');
load('enzymedataSEC.mat');
load('enzymedataDummyER.mat');
load('enzymedataMachine.mat')
%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1810',0,'b');% glycine production
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;

clear tot_protein f_modeled_protein;
factor_k = 1;
%% Solve LPs
mu_list = [0.1 0.2 0.3 0.35 0.37 0.38 0.39 0.4];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
%     f_carbon = 5.244714732847007e-01*mu;
%     f_carbon = f_carbon/0.38067
%     model_tmp = changeBiomass(model,f_carbon,'r_4041','s_3718[c]');
     model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,num2str(mu*100));
    allName{i} = fileName;
end

% 
writeclusterfileLP(allName(:),'sub_1')
cd ../
