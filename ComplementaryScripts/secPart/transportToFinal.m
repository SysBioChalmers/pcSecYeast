function [newModel,peptide_name,rxns] = transportToFinal(model,peptide,peptide_org,compartment,onlyrxns)
rxns = [];
if strcmp(compartment,'er')==1 || strcmp(compartment,'erm') ==1
    Reaction = transportFromGolgiToER(peptide,peptide_org,compartment);
elseif strcmp(compartment,'vm')==1
    Reaction = transportFromGolgiToVM(peptide,peptide_org);
elseif strcmp(compartment,'v')==1
    Reaction = transportFromGolgiToV(peptide,peptide_org);
elseif strcmp(compartment,'ce')==1
    Reaction = transportFromGolgiToCe(peptide,peptide_org);
elseif  strcmp(compartment,'e')==1
    Reaction = transportFromGolgiToS(peptide,peptide_org);
else
    Reaction = transportFromGolgiToOther(peptide,peptide_org,compartment);
end 
if exist('Reaction','var')
    for i=1:length(Reaction)
        if onlyrxns == 1
            rxns = [rxns;{Reaction{i}.rxns}];
        else
            model=addYeastReaction(model,Reaction{i}.eq,{Reaction{i}.rxns},{Reaction{i}.rxnNames});
            rxns = [rxns;{Reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s',peptide); %here the peptide name doesn't change to the final metaboloites due to the reason to link the nexr NGmisfolding reactions
else
    newModel = model;
    peptide_name = peptide;
end
end