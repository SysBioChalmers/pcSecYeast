function [newModel,peptide_name,rxns] = addDSB(model,peptide,Length_total,DSB,onlyrxns)
rxns = [];
if DSB >0
    Length = Length_total/40;
    reaction{1}.rxns = sprintf('%s_DSB_sec_PDI_BIP_NEFS_complex',peptide);
    reaction{2}.rxns = sprintf('%s_DSB_PDI_II_sec_PDI1_ERV2_Ero1p_complex',peptide);
    reaction{3}.rxns = sprintf('%s_DSB_PDI_III_sec_PDI1_ERV2_Ero1p_complex',peptide);


    reaction{1}.rxnNames = sprintf('%s_DSB_PDI_BIP_NEFS_complex Formation of folding complex with Kar2ATP and disulfide bond formation',peptide);
    reaction{2}.rxnNames =sprintf('%s_DSB_PDI_II_PDI1_ERV2_Ero1p_complex Formation of folding complex with Kar2ATP and disulfide bond formation',peptide);
    reaction{3}.rxnNames =sprintf('%s_DSB_PDI_III_PDI1_ERV2_Ero1p_complex Formation of folding complex with Kar2ATP and disulfide bond formation',peptide);


    reaction{1}.eq = sprintf('%s[er] + %.15f ATP[er] + %.15f H2O[er] => %s_Kar2ATPcplx[er] + %.15f ADP[er] + %.15f H+[er] + %.15f phosphate[er]',peptide,Length,Length,peptide,Length,Length,Length);
    reaction{2}.eq=  sprintf('%s_Kar2ATPcplx[er] + %.15f PDI-ox[er] => %s-ds-PDI-ox[er]',peptide,DSB,peptide);
    reaction{3}.eq=  sprintf('%s-ds-PDI-ox[er] => %s_DSB[er] + %.15f PDI[er]',peptide,peptide,DSB);

    for i=1:3
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s_DSB',peptide);
else
    newModel = model;
    peptide_name = peptide;

end
