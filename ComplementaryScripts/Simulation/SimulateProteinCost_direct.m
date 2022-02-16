% This function would calculate the energetic cost for Secretory proteins,
% is comparable to the SimulateProteinCost.m. This function only takes into
% account the energetic for itselt, while the other one considers the
% indirect cost for the corresponding shares of catalyzing enzymes under
% the proteome limitation condition

[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
%ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strrep(ERprotein,'-','');
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');

% load model and param
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('enzymedataDummyER.mat');
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
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
factor_k = 1;
f_unmodelER = 0.046;
f_erm = 0.083;
mkdir('SimulateProteincost_direct')
cd SimulateProteincost_direct
allname = cell(0,1);
mulist = [0.05,0.1];
ratios = [5E-7,1E-6,5E-6,1E-5,5E-5,1E-4];
cost = [];
for i = 1:length(ERprotein)
    i
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(i),true);
    mu = mulist(1);
    model_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
    for j = 1:length(ratios)
        ratio = ratios(j);
        model_tmp = changeRxnBounds(model_tmp,[ERprotein{i},' exchange'],ratio,'b');
        model_tmp = changeObjective(model_tmp,'r_1714',1);% acetaldehyde production
        sol = optimizeCbModel(model_tmp);
        cost(i,j,1) = sol.f;
    end
    p = polyfit(ratios,abs(cost(i,:,1)),1);
    res_glc_final(i,1) = p(1);
end

save('directCost_res.mat')
cd ../
