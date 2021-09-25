function ReadProteinCostResult(a,b)

initcluster
cd('SimulateProteinCost')
[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strrep(ERprotein,'-','');
ERprotein = strcat(ERprotein,'new');
%TP = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'HGCSF'};
%ERprotein = [TP;ERprotein];

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
f_erm = 0.083;
allname = cell(0,1);
mulist = [0.1,0.2];
ratios = [3E-7,6E-7,9E-7,3E-6,6E-6];
file = dir('*.out');
filename = {file.name};
res_riboCost = zeros((b-a),1);
res_TP = zeros((b-a),1);
res_mu = zeros((b-a),1);
res_slope = zeros((b-a),1);
res_glc = zeros((b-a),1);
res_glcCost = zeros((b-a),1);
for i = a:b
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(i),true);
    display([num2str(i) '/' num2str(length(ERprotein))]);

    
    fluxes = zeros(length(model_tmp.rxns),0);
    allfile = filename(contains(filename,strrep(ERprotein{i},'new','')));
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
    
    ribo_ex = 'Mach_Ribosome_complex_dilution';
    
    mu = round(fluxes(ismember(model_tmp.rxns,'r_2111'),:),2);
    munique = unique(mu);

    for k = 1:length(munique)
        fluxes_tmp = fluxes(:,mu == munique(k));
    ribo = fluxes_tmp(ismember(model_tmp.rxns,ribo_ex),:);
    tp = fluxes_tmp(ismember(model_tmp.rxns,[ERprotein{i},' exchange']),:);
    glc = fluxes_tmp(ismember(model_tmp.rxnNames,'D-glucose exchange'),:);
    res_glcCost(i,1:length(tp),k) = glc;
    res_riboCost(i,1:length(tp),k) = ribo/munique(k);
    res_TP(i,1:length(tp),k) = tp;
    res_mu(i,1:length(tp),k) = munique(k);
    p = polyfit(tp,ribo/munique(k),1);
    res_slope_ribo(i,1,k) = p(1);
    p = polyfit(tp,glc,1);
    res_slope_glc(i,1,k) = p(1);
    end
    
    
end
 save([ERprotein{a},'_',num2str(a),'.mat'],'res_riboCost','res_TP','res_mu','res_slope_ribo','a','b','ERprotein','res_glcCost','res_slope_glc');
cd ../
end