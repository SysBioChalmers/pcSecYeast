%% SimulateTP
% max production rate under various growth rates on minimal media
TP = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'HGCSF';'Insulin2'}; % Insulin2 is the Insulin without Disulfide bonds

mulist = [0.02:0.02:0.25,0.25:0.005:0.3,0.32 0.34 0.36 0.38 0.4];

mkdir('SimulateTP_SDAA');
cd SimulateTP_SDAA/;
allname = cell(0,1);

%% Solve LPs
for i = 1:length(TP)
    % load model and param
    load('enzymedata.mat')
    load('enzymedataMachine.mat')
    load('enzymedataSEC.mat')
    load('pcSecYeast.mat')
    load('enzymedataDummyER.mat');
    
    %% Set model
    % set medium
    model = setMedia(model,4);% SDAA
    % set carbon source
    model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
    % set oxygen
    model = changeRxnBounds(model,'r_1992',-1000,'l');
    % block reactions
    model = blockRxns(model);
    model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
    model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
    model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
    rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
    model.ub(rxn) = 0;
    %% Set optimization
    
    tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
    f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
    % r_4041 is pseudo_biomass_rxn_id in the GEM
    % s_3717[c] is protein id
    
    f = tot_protein * f_modeled_protein;
    f_mito = 0.1;
    f_unmodelER = 0.046;
    f_erm = 0.083;
    clear tot_protein f_modeled_protein;
    
    factor_k = 1; % doesn't matter this is to tune the saturation
    rxnID = strcat(TP{i},' exchange');
    osenseStr = 'Maximize';
    [model_tmp,enzymedataTP] = addTargetProtein(model,TP(i));
    model = model_tmp;
    save(['model',TP{i},'.mat'],'model')
    [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);
    enzymedata_new = enzymedata;
    enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
    enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
    enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
    enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
    enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
    enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
    enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
    enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
    enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    
    for j = 1:length(mulist)
        disp(num2str(j))
        mu = mulist(j);
        
        model_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
        name = [TP{i} '_',num2str(mu*100)];
        fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);% no extra constraint
        allname = [allname;{fileName}];
    end
end

%% reference condition to facilitate further identification of engenieering
% targets
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');
model = setMedia(model,4);% SDAA
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_2033',0,'b');% acetaldehyde production

rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c');
model.ub(rxn) = 0;

% Set optimization
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

enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
model.ub(contains(model.rxns,'dilution_misfolding')) = 0; % block the accumulation in the model;

for j = 1:length(mulist)
    disp(num2str(j))
    mu = mulist(j);
    
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    name = ['ref_',num2str(mu*100)];
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,name);% no extra constraint
    allname = [allname;{fileName}];
end
% write cluster file
writeclusterfileLP(allname,'sub_TP1')
cd ../


