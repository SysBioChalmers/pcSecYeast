function SimulateGlcConWithoutSec(a,b)
% Simulate chemostat by simulate the carbon concentration

initcluster
% mkdir('SimulateGlcCon')
% cd('SimulateGlcCon')
%197/s HXT7 KM 2.5mM pubmed: 11561293 15345416
% 53/s HXT2  pubmed:10191260
% https://www.jbc.org/action/showPdf?pii=S0021-9258%2819%2973030-0

res_glc = [0.052631579 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.65 0.8 1 5 10 20 50 100 150 200 500 1000];
Km_hxt7 = 2.5;
Km_hxt1 = 100;

% load model and param
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeastWithoutSEC.mat')
load('enzymedataDummyER.mat');

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
f_unmodelER = 0.046;
clear tot_protein f_modeled_protein;

factor_k = 1;
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id
for i = a:b
    enzymedata.kcat(contains(enzymedata.enzyme,'r_1166')) = 0;
    kcat_hxt7 = 720000*(res_glc(i)/(Km_hxt7 + res_glc(i)));
    kcat_hxt1 = 720000*(res_glc(i)/(Km_hxt1 + res_glc(i)));
    enzymedata.kcat(contains(enzymedata.enzyme,'r_1166_3')) = kcat_hxt7;
    enzymedata.kcat(contains(enzymedata.enzyme,'r_1166_10')) = kcat_hxt1;
        
        factor_mu_low = 0.005;
        factor_mu_high = 0.38;
        
        while factor_mu_high-factor_mu_low > 0.001
            factor_mu_mid = (factor_mu_low+factor_mu_high)/2;
            mu = factor_mu_mid;
            disp(['Without sf: D = ' num2str(mu) '; factor_glcKcat = ' num2str(res_glc(i))]);
        
            model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
            name = ['GlcCon_withoutSec',num2str(res_glc(i)),'_',num2str(mu*100)];
            fileName = writeLPWithoutSec(model_tmp,mu,f,f_mito,f_unmodelER,osenseStr,rxnID,enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER,factor_k,name);
  
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
        kcat = [kcat_hxt7,kcat_hxt1];
        glc_conc_without_sf(1) = factor_mu_mid;
        glc_conc_without_sf(2) = res_glc(i);
        save(['reswithoutSec',num2str(res_glc(i)),'.mat'],'fluxes_simulated_without_sf','glc_conc_without_sf','kcat')
end
end
