function reaction = transportFromGolgiToVM(peptide,peptide_org)

%ALP pathway
reaction{1}.rxns = sprintf('%s_ALPtransport_sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_Direct vacuol transit pathway_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p',peptide);


reaction{1}.eq = sprintf('%s[g] + 4 GTP[c] + 4 H2O[c] => %s_folding[vm] + 4 GDP[c] + 4 phosphate[c] + 4 H+[c]',peptide,peptide_org);

end