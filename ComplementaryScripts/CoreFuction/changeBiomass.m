%% changeBiomass 
function model = changeBiomass(model,f_modeled_protein,pseudo_biomass_rxn_id,protein_id)


rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);

% change protein coefficient
metidx = ismember(model.mets,protein_id);
model.S(metidx,rxnidx) = -1 * (1-f_modeled_protein);


