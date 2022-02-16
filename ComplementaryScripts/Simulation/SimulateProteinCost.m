function SimulateProteinCost(a,b)
% this function is to calculate the protein cost, ER cost of synthesizing each
% protein

initcluster
[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
%ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strrep(ERprotein,'-','');
ERprotein = strcat(ERprotein,'new');
protein_info(:,2) = strcat(protein_info(:,2),'new');
%TP = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'HGCSF'};
%ERprotein = [TP;ERprotein];

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
mkdir('SimulateProteinCost')
cd SimulateProteinCost
allname = cell(0,1);
mulist = [0.05,0.1];
ratios = [5E-7,1E-6,5E-6,1E-5,2E-5];
for i = a:b
    [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein(i),true);
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


    for k = 1:length(mulist)
        mu = mulist(k);
        model_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
    for j = 1:length(ratios)
        ratio = ratios(j);
        model_tmp = changeRxnBounds(model_tmp,[ERprotein{i},' exchange'],ratio,'b');
    
    % Set optimization
    rxnID = 'dilute_dummy'; %minimize glucose uptake rate
    osenseStr = 'Maximize';

    name = [strrep(ERprotein{i},'new',''),'_',num2str(mu*100),'_ratio',num2str(1/ratio)];
    fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);
    allname = [allname;{fileName}];
    end
    end
end

% write cluster file
writeclusterfileLP(allname,['sub_proteinCost_',num2str(a)])
display([num2str(a) '/' num2str(length(ERprotein))]);
cd ../
end



% [model_tmp,enzymedataTP] = addTargetProtein(model,ERprotein);
% enzymedata_new = enzymedata;
% enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
% enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
% enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
% enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
% 
% allgenes = enzymedata_new.proteins;
% 
% for i = 1:length(enzymedata.enzyme)
%    [~,idx] = ismember(strrep(enzymedata.enzyme(i),'_complex',''),model.rxns);
%    subunitlist = enzymedata.subunit(i,:);
%    subunit_num = cell2mat(cellfun(@any,subunitlist,'UniformOutput',false));
%    subunitlist = subunitlist(1,subunit_num);
%    subunitcoef = enzymedata.subunit_stoichiometry(i,subunit_num);
%    [~,idx2] = ismember(subunitlist,allgenes);
%    matrix(idx2,idx) = subunitcoef./enzymedata.kcat(i);
% end
% 
% 
% for i = 1:length(enzymedataSEC.enzyme)
%    idx = find(endsWith(model.rxns,['_',enzymedataSEC.enzyme{i}]));
%    subunitlist = enzymedataSEC.subunit(i,:);
%    subunit_num = cell2mat(cellfun(@any,subunitlist,'UniformOutput',false));
%    subunitlist = subunitlist(1,subunit_num);
%    subunitcoef = enzymedataSEC.subunit_stoichiometry(i,subunit_num);
%    [~,idx2] = ismember(subunitlist,allgenes);
%    matrix(idx2,idx) = subunitcoef./enzymedataSEC.kcat(i,13); % max one
% end
% 
% % Ribosome and Ribosome assembly
% trans_rxns = model.rxns(endsWith(model.rxns,'_translation'));
% for i = 1:length(trans_rxns)
%     rxn_id = trans_rxns{i};
%     comp_name = strrep(rxn_id,'_peptide_translation','');
%     comp_name = strrep(comp_name,'r_','');
%     idx = find(strcmp(model.rxns,rxn_id));
%     prot_leng = enzymedata_all.proteinLength(ismember(enzymedata_all.proteins,strrep(comp_name,'_','-')));
%     subunitlist = enzymedataMachine.subunit(1,:);
%     subunit_num = cell2mat(cellfun(@any,subunitlist,'UniformOutput',false));
%     subunitlist = subunitlist(1,subunit_num);
%     subunitcoef = enzymedataMachine.subunit_stoichiometry(1,subunit_num);
%     [~,idx2] = ismember(subunitlist,allgenes);
%     matrix2(idx2,idx) = prot_leng.*subunitcoef./enzymedataMachine.kcat(1); % ribosome
% end
%     
