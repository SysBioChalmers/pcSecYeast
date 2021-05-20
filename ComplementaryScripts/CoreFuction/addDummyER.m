%% addDummy
function model = addDummyER(model,pseudo_protein_rxn_id,protein_id,protein_info)
%% Add translation for the dummy complex
% Assuming that the dummy complex has the same AA composition as biomass
% protein, then we can just copy the protein pseudoreaction from the GEM
% and change the product to the dummy complex.
rxn_idx = ismember(model.rxns,pseudo_protein_rxn_id);

coeflist_tmp = full(model.S(:,rxn_idx));
metlist = model.mets(coeflist_tmp ~= 0);

metlist(ismember(metlist,protein_id)) = {'dummyER_peptide[c]'};
coeflist = coeflist_tmp(coeflist_tmp ~= 0);
coeflist(~ismember(metlist,'dummyER_peptide[c]')) = 100*coeflist(~ismember(metlist,'dummyER_peptide[c]'));
rxnid = 'translate_dummyER';
model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);

%% Add dilution for the dummy complex
rxnid = 'dilute_dummyER';
metlist = {'dummyER_folding[c]'};
coeflist = -1;
model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);

%% Add Modification
model = addfolding(model,{'dummyER'},protein_info);% add folding process for those peptide xx_peptide --> xx_folding

%% add an extra exchange rxn for dummyER_sp since it doesn't have the info for the SP length
metid = {'dummyER_sp[c]'};
rxnid = 'dummyER_sp exchange';
model = addReaction(model,rxnid,'metaboliteList',metid,'stoichCoeffList',-1,'reversible',false);

end