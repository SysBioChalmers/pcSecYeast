function [newModel,peptide_name,rxns] = addOG(model,peptide,peptide_org,Length_total,OG,onlyrxns)
rxns = [];
if OG >0
    Length = OG;
    
    reaction{1}.rxns = sprintf('%s_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_Pmt3p_complex',peptide_org);
  
    reaction{1}.rxnNames = sprintf('%s_OG_EROG_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_Pmt3p_complex ER O-glycosylation',peptide);
     
    reaction{1}.eq = sprintf('%s[er] + %d dolichyl D-mannosyl phosphate[er] => %s_OG_M1[er] + %d dolichyl phosphate[er]',peptide,OG,peptide,OG);
    
    for i=1:1
         if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
         else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
         end
    end
    newModel = model;
    peptide_name = sprintf('%s_OG_M1',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end