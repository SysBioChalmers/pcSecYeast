%% extractModeledprotein 
function f_modeled_protein = extractModeledprotein(model,pseudo_biomass_rxn_id,protein_id)

metidx = ismember(model.mets,protein_id);
rxnidx = ismember(model.rxns,pseudo_biomass_rxn_id);
f_modeled_protein = 1 + full(model.S(metidx,rxnidx));
end