function reaction = transportFromGolgiToCe(peptide,peptide_org)

%LDSV pathway
reaction{1}.rxns = sprintf('%s_LDSV_sec_Arf1p_Sec3p_Sec5p_Sec6p_Sec8p_Sec10p_Sec15p_Exo70p_Exo84p_Sec4p_Chc1p_Clc1p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_LDSV_Arf1p_Sec3p_Sec5p_Sec6p_Sec8p_Sec10p_Sec15p_Exo70p_Exo84p_Sec4p_Chc1p_Clc1p',peptide);


reaction{1}.eq = sprintf('%s[g] + GTP[c] + H2O[c] => %s_folding[ce] + GDP[c] + phosphate[c] + H+[c]',peptide,peptide_org);

end