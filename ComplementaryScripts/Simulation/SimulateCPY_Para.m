function SimulateCPY_Para(TPrate,misP,n,accstate,extraconstrait)
%% Sismulate CPY Para
% This function is for Supplementary Figure to show the parameter
% importance of Retrotranlocation. 
initcluster
%extraconstrait = [0,4,5,6,7]; % 4 means total er constraint, 5 means total erm constraint, 6 means total retro-tranloc constraint, 7 means total ERAD constraint
% misP = [0,0.1,0.3,0.5,0.7,0.9]
TP = 'YMR297W';

mkdir('SimulateCPY_Para');
cd SimulateCPY_Para/;

load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');

model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

% close other protein dilution misfolding 
 rxn = contains(model.rxns,'_dilution_misfolding_m')|contains(model.rxns,'_dilution_misfolding_c')|contains(model.rxns,'_dilution_misfolding_er');
 model.ub(rxn) = 0;

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
%% Write LPs in three condition
allName = cell(0,1);

% add TP
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';
[model,enzymedataTP] = addTargetProtein(model,{TP});
% close refold rxns
refold_list = model.rxns(contains(model.rxns,'refolding_'));
model = changeRxnBounds(model,refold_list,0,'b');

% change the mannose number in the CPY according to the literature median
% number for mannose in the N-glycan is around 15 
rxn = find(contains(model.rxns,'YMR297W_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex'));
metname = {'GDP-alpha-D-mannose [Golgi]';'GDP [Golgi]'};
[~,idx] = ismember(metname,model.metNames);
model.S(idx,rxn) = 1/9*model.S(idx,rxn);
rxn = find(contains(model.rxns,'YMR297W_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex'));
metname = {'GDP-alpha-D-mannose [Golgi]';'GDP [Golgi]'};
[~,idx] = ismember(metname,model.metNames);
model.S(idx,rxn) = 1/30*model.S(idx,rxn);

% change the loop
rxn = find(contains(model.rxns,'cycle_accumulation_sec_pdi1p_ero1p_complex'));
metname = {'oxygen [endoplasmic reticulum]';'glutathione [endoplasmic reticulum]';'glutathione disulfide [endoplasmic reticulum]';'hydrogen peroxide [endoplasmic reticulum]'};
[~,idx] = ismember(metname,model.metNames);
model.S(idx,rxn) = n/10.*model.S(idx,rxn);
rxn = find(contains(model.rxns,'cycle_accumulation_sec_acc_Kar2p_complex'));
metname = {'H+ [endoplasmic reticulum]';'H2O [endoplasmic reticulum]';'phosphate [endoplasmic reticulum]';'ATP [endoplasmic reticulum]';'ADP [endoplasmic reticulum]'};
[~,idx] = ismember(metname,model.metNames);
model.S(idx,rxn) = n/10.*model.S(idx,rxn);


[enzymedataTP] = SimulateRxnKcatCoef(model,enzymedataSEC,enzymedataTP);
% change the coef for the enzyme usage
rxn = find(contains(enzymedataTP.rxns,'YMR297W_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex'));
enzymedataTP.rxnscoef(rxn) = 1/9*enzymedataTP.rxnscoef(rxn);
rxn = find(contains(enzymedataTP.rxns,'YMR297W_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex'));
enzymedataTP.rxnscoef(rxn) = 1/30*enzymedataTP.rxnscoef(rxn);
enzymedataTP.kdeg = misP;
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

x = find(contains(enzymedata_new.rxns,'cycle_accumulation_sec_acc_Kar2p_complex'));
enzymedata_new.rxnscoef(x) = enzymedata_new.rxnscoef(x).*n/10;
x = find(contains(enzymedata_new.rxns,'cycle_accumulation_sec_pdi1p_ero1p_complex'));
enzymedata_new.rxnscoef(x) = enzymedata_new.rxnscoef(x).*n/10;

% change transport
if accstate
    model.ub(strcmp(model.rxns,'YMR297W_degradation_misfolding_c')) = 0;
end

if accstate
state = 'accon';
else
state = 'accauto';
end

TPrxnID = strcat('r_',TP,'_peptide_translation');
model.lb(strcmp(model.rxns,TPrxnID)) = TPrate;
factor_mu_low = 0.1;
factor_mu_high = 0.4;
flux_tmp = zeros(length(model.rxns),1);
while factor_mu_high-factor_mu_low > 0.001
    factor_mu_mid = (factor_mu_low+factor_mu_high)/2;
    mu = factor_mu_mid;
    disp(['Without sf: D = ' num2str(mu) '; TP = ' TP]);
    
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    name = [TP,'_',num2str(TPrate*1E6),'_mis',num2str(misP*10),'_',num2str(mu*100),'_',num2str(n),'_cons',num2str(extraconstrait),state];

    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name,[extraconstrait]);
    %command = sprintf('/home/f/feiranl/tools/soplex-4.0.0/build/bin/soplex -s0 -g5 -t3000 -f1e-17 -o1e-17 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    command = sprintf('/cephyr/users/feiranl/Hebbe/tools/build/bin/soplex -s0 -g5 -t3000 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = [fileName,'.out'];
    [~,solME_status,solME_full] = readSoplexResult(fileName_out,model_tmp);
    
    if strcmp(solME_status,'optimal') && ~isempty(solME_full)
        factor_mu_low = factor_mu_mid;
        flux_tmp = solME_full;
    else
        factor_mu_high = factor_mu_mid;
    end
end
fluxes_simulated = flux_tmp;
res(1) = mu;
res(2) = TPrate;
res(3) = misP;
res(4) = n;
res(5) = extraconstrait;
resid = {'mu','TPrate','misP','n','extraconstrait'};

save(['res',TP,'_',num2str(TPrate*1E6),'_mis',num2str(misP*10),'_',num2str(n),'extracons',num2str(extraconstrait),'.mat'],'res','fluxes_simulated','resid')

end



