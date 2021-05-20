%% addDummy 
function model = addDummy(model,pseudo_protein_rxn_id,protein_id)
%% Add translation for the dummy complex
% Assuming that the dummy complex has the same AA composition as biomass
% protein, then we can just copy the protein pseudoreaction*100 from the GEM
% and change the product to the dummy complex.

rxn_idx = ismember(model.rxns,pseudo_protein_rxn_id);

coeflist_tmp = full(model.S(:,rxn_idx));
metlist = model.mets(coeflist_tmp ~= 0);

metlist(ismember(metlist,protein_id)) = {'dummy[c]'};
coeflist = coeflist_tmp(coeflist_tmp ~= 0);
coeflist(~ismember(metlist,'dummy[c]')) = 100*coeflist(~ismember(metlist,'dummy[c]'));

rxnid = 'translate_dummy';
model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);

%% Add dilution for the dummy complex
rxnid = 'dilute_dummy';
metlist = {'dummy[c]'};
coeflist = -1;
model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
