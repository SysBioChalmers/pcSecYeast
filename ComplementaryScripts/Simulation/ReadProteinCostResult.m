function ReadProteinCostResult(a,b)

initcluster
mkdir('SimulateProteinCostResult')
cd('SimulateProteinCost')
[~,~,protein_info] = xlsread('../../ComplementaryData/Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
ERprotein = [ERprotein;{'Amylase';'Insulin';'HBsAg'}];
% load model and param
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')

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

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
factor_k = 1;

allname = cell(0,1);
mulist = [0.05,0.1];
ratios = [1E-6,1E-5,5E-5,1E-4,1E-3];
file = dir('*.out');
filename = {file.name};
res_riboCost = zeros((b-a),1);
res_TP = zeros((b-a),1);
res_mu = zeros((b-a),1);
res_slope = zeros((b-a),1);
res_glc = zeros((b-a),1);
for i = a:b
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(i));
    display([num2str(i) '/' num2str(length(ERprotein))]);

    
    fluxes = zeros(length(model_tmp.rxns),0);
    allfile = filename(contains(filename,ERprotein{i}));
    for j = 1:length(allfile)
        [~,sol_status,sol_full] = readSoplexResult(allfile{j},model_tmp);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
        else
            fluxes = [fluxes zeros(length(model_tmp.rxns),1)];
        end
    end
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    
    ribo_ex = 'Ribosome_complex_dilution';
    
    mu = round(fluxes(ismember(model_tmp.rxns,'r_2111'),:),2);
    munique = unique(mu);
    for k = 1:length(munique)
        fluxes_tmp = fluxes(:,mu == munique(k));
    ribo = fluxes_tmp(ismember(model_tmp.rxns,ribo_ex),:);
    tp = fluxes_tmp(ismember(model_tmp.rxns,[ERprotein{i},' exchange']),:);
    glc = fluxes_tmp(ismember(model_tmp.rxns,'D-glucose exchange'),:);
    res_riboCost(i,1:length(tp),k) = ribo/munique(k);
    res_TP(i,1:length(tp),k) = tp;
    res_mu(i,1:length(tp),k) = munique(k);
    p = polyfit(tp,ribo/munique(k),1);
    res_glc(i,1:length(tp),k) = glc;
    res_slope(i,1,k) = p(1);
    end
   
end
 save(['../ProteinCostResult/',ERprotein{a},'_',num2str(a),'.mat'],'res_riboCost','res_TP','res_mu','res_slope','a','b','ERprotein','res_glc');
cd ../