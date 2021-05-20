function [newModel,peptide_name,rxns] = mature(model,peptide,onlyrxns)
rxns = [];
reaction{1}.rxns = sprintf('%s_Mature',peptide);
reaction{1}.rxnNames = sprintf('%s Mature',peptide);
reaction{1}.eq = sprintf('%s[g] => %s_mature[g]',peptide,peptide);
if onlyrxns == 1
    rxns = [rxns;{reaction{1}.rxns}];
else
    model=addYeastReaction(model,reaction{1}.eq,{reaction{1}.rxns},{reaction{1}.rxnNames});
    rxns = [rxns;{reaction{1}.rxns}];
end
 newModel = model;
 peptide_name = sprintf('%s_mature',peptide);

end
