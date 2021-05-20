function [abun_complex_id,abun_protein_id,abun_complex,abun_protein,abun_complex_mg,abun_protein_mg] = abundanceCalculation(model,flux)
% fileName is the lp.out file
% abun_comples 1st column is in mmol 2 column is in mg
% abun_protein 1st column in mmol column 2nd column in mg
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('enzymedata.mat')
load('enzymedataDummyER.mat')
% [sol_obj,sol_status,sol_full] = readSoplexResult(fileName,model);
sol_full = flux;
mu = sol_full(strcmp(model.rxns,'r_2111'),:);

enzymedata_all.enzyme = [enzymedata.enzyme;enzymedataSEC.enzyme;enzymedataMachine.enzyme];
enzymedata_all.enzyme_MW = [enzymedata.enzyme_MW;enzymedataSEC.enzyme_MW;enzymedataMachine.enzyme_MW];
enzymedata_all.proteins = [enzymedata.proteins;enzymedataSEC.proteins;enzymedataMachine.proteins];
enzymedata_all.proteinMWs = [enzymedata.proteinMWs;enzymedataSEC.proteinMWs;enzymedataMachine.proteinMWs];
enzymedata_all.proteinLength = [enzymedata.proteinLength;enzymedataSEC.proteinLength;enzymedataMachine.proteinLength];
enzymedata_all.subunit = [enzymedata.subunit;enzymedataSEC.subunit;enzymedataMachine.subunit];
enzymedata_all.subunit_stoichiometry = [enzymedata.subunit_stoichiometry;enzymedataSEC.subunit_stoichiometry;enzymedataMachine.subunit_stoichiometry];
dil_rxns = model.rxns(contains(model.rxns,'complex_dilution'));
abun_protein_id = unique(enzymedata_all.proteins,'stable');
abun_complex = zeros(length(dil_rxns),length(flux(1,:)));
abun_protein = zeros(length(abun_protein_id),length(flux(1,:)));
abun_protein_mg = zeros(length(abun_protein_id),length(flux(1,:)));
abun_complex_mg = zeros(length(dil_rxns),length(flux(1,:)));
abun_complex_id = strrep(dil_rxns,'_dilution','');


    for i = 1:length(dil_rxns)
        rxn_id = dil_rxns{i};
        comp_name = strrep(rxn_id,'_dilution','');
        idx = strcmp(model.rxns,rxn_id);
        MW = enzymedata_all.enzyme_MW(contains(enzymedata_all.enzyme,comp_name));
        if isempty(MW)
            MW = enzymedata_all.enzyme_MW(contains(enzymedata_all.enzyme,strrep(comp_name,'_complex','')));
        end
        coeff = MW(1);
        abun_complex(i,:) = sol_full(idx,:)./mu; % mmol
        abun_complex_mg(i,:) = sol_full(idx,:)./mu*coeff; % mg
        % find subunit abundance
        idx2 = strcmp(enzymedata_all.enzyme,comp_name) | strcmp(enzymedata_all.enzyme,strrep(comp_name,'_complex',''));
        subunit = enzymedata_all.subunit(idx2,:);
        subunit_stoi = enzymedata_all.subunit_stoichiometry(idx2,:);
        [~,idx3] = ismember(subunit,abun_protein_id);
        idx3 = idx3(idx3~=0);
        abun_protein(idx3,:) = abun_protein(idx3,:) + transpose(abun_complex(i,:)'.*subunit_stoi(1:length(idx3)));
        abun_protein_mg(idx3,:) = abun_protein(idx3,:).*enzymedata_all.proteinMWs(idx3);
    end
    
    
% find dummy protein abun
[~,idx] = ismember({'dilute_dummy','dilute_dummyER'},model.rxns);
abun_protein = [abun_protein;sol_full(idx,:)./mu];
abun_complex = [abun_complex;sol_full(idx,:)./mu];
abun_protein_id = [abun_protein_id;{'dummy';'dummyER'}];
abun_complex_id = [abun_complex_id;{'dummy';'dummyER'}];
abun_complex_mg = [abun_complex_mg;(sol_full(idx,:)./mu).*46000];
abun_protein_mg = [abun_protein_mg;(sol_full(idx,:)./mu).*46000];
end


