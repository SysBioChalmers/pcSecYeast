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

%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
f_erm = 0.0083;
factor_k = 1;
%% Solve LPs
mu_list = [0.1,0.15,0.2,0.25,0.28,0.3,0.32,0.34,0.35,0.36,0.37,0.38];
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);

 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
%     f_carbon = 5.244714732847007e-01*mu;
%     f_carbon = f_carbon/0.38067
%     model_tmp = changeBiomass(model,f_carbon,'r_4041','s_3718[c]');
     model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    %fileName = writeLP(model_tmp,mu,f,f_mito,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,num2str(mu*100));
    fileName = writeLP(model_tmp,mu,f,f_mito,f_unmodelER,f_erm,osenseStr,rxnID,enzymedata_all,factor_k,num2str(mu*100));

    allName{i} = fileName;
end

% 
writeclusterfileLP(allName(:),'sub_1')
cd ../

%% Crabtree withouSec %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
mkdir('crabtreeWithoutSec')
cd('crabtreeWithoutSec')
% load model
load('pcSecYeastWithoutSEC.mat');
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

%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = 0.046;
f_mito = 0.1;
clear tot_protein f_modeled_protein;

factor_k = 1;
%% Solve LPs
mu_list = [0.1,0.15,0.2,0.25,0.27,0.28,0.3,0.32,0.34,0.35,0.36,0.37,0.38];

fluxes = zeros(length(model.rxns),length(mu_list));
 allName = cell(0,1);
for i = 1:length(mu_list)
    mu = mu_list(i);
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    fileName = writeLPWithoutSec(model_tmp,mu,f,f_mito,f_unmodelER,osenseStr,rxnID,enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER,factor_k,num2str(mu*100));
    allName{i} = fileName;
end

% 
writeclusterfileLP(allName(:),'sub_1')
