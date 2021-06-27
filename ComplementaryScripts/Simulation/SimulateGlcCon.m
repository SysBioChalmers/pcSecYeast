function SimulateGlcCon(i)
% Simulate chemostat by simulate the carbon concentration

initcluster
% mkdir('SimulateGlcCon')
% cd('SimulateGlcCon')
%197/s HXT7 KM 2.5mM pubmed: 11561293 15345416
% 53/s HXT2  pubmed:10191260
% https://www.jbc.org/action/showPdf?pii=S0021-9258%2819%2973030-0
mkdir('SimulateGlcCon')
cd SimulateGlcCon
res_glc = [1E-6 1E-5 1E-4 1E-3 1E-2 0.1 0.2 0.5 1 2 3 4 5 10 20 50 100 150 500];
Km_hxt = [110;1.5;34;9.3;2.5]; % bionumber 110954 110739 PMID 2482015 10336421 https://doi.org/10.1101/2020.06.22.165753 table
hxt = {'YHR094C';'YMR011W';'YDR345C';'YHR092C';'YDR342C'};
kcatmax = [1012 53 479 155 197]*3600; % calculated from Vmax based on the assumption that the kca tfor the hxt2 is 53/s and hxt7 is 200/s
hxtrxnID = {'r_1166_10_complex';'r_1166_17_complex';'r_1166_5_complex';'r_1166_9_complex';'r_1166_3_complex'}; % glucose transporter rxn
% load model and param
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');
% close dil rxn for misfolding
dilrxn = contains(model.rxns,'_dilution_misfolding');
model.ub(dilrxn) = 0;
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

rxnID = 'dilute_dummy'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
f = tot_protein * f_modeled_protein;
f_mito = 0.1;
f_unmodelER = tot_protein * 0.046;
clear tot_protein f_modeled_protein;
f_erm = 0.008;
factor_k = 1;
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id
lowkcat = 3600;
enzymedata.kcat(enzymedata.kcat< lowkcat) = lowkcat;

enzymedata.kcat(contains(enzymedata.enzyme,'r_1166')) = 0;

kcat_hxt = kcatmax.*(res_glc(i)./(Km_hxt + res_glc(i)))';
[~,Idx] = ismember(hxtrxnID,enzymedata.enzyme);
enzymedata.kcat(Idx) = kcat_hxt;

enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);

factor_mu_low = 0.005;
factor_mu_high = 0.38;


while factor_mu_high-factor_mu_low > 0.001
    factor_mu_mid = (factor_mu_low+factor_mu_high)/2;
    mu = factor_mu_mid;
    disp(['Without sf: D = ' num2str(mu) '; factor_glcKcat = ' num2str(res_glc(i))]);
    
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    name = ['GlcCon_',num2str(res_glc(i)),'_',num2str(mu*100)];
    fileName = writeLP(model_tmp,mu,f,f_mito,f_unmodelER,f_erm,osenseStr,rxnID,enzymedata_all,factor_k,name);
    %command = sprintf('/home/f/feiranl/tools/soplex-4.0.0/build/bin/soplex -s0 -g5 -t3000 -f1e-17 -o1e-17 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    command = sprintf('/cephyr/users/feiranl/Hebbe/tools/build/bin/soplex -s0 -g5 -t3000 -f1e-17 -o1e-17 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = [fileName,'.out'];
    [~,solME_status,solME_full] = readSoplexResult(fileName_out,model_tmp);
    
    if strcmp(solME_status,'optimal')
        factor_mu_low = factor_mu_mid;
        flux_tmp = solME_full;
    else
        factor_mu_high = factor_mu_mid;
    end
end
fluxes_simulated_without_sf = flux_tmp;
glc_conc_without_sf(1) = factor_mu_mid;
glc_conc_without_sf(2) = res_glc(i);
save(['res',num2str(res_glc(i)),'.mat'],'fluxes_simulated_without_sf','glc_conc_without_sf','kcat_hxt')
end


