function reaction = transportFromGolgiToER(peptide,peptide_org,compartment)

% add 3 common reactions
reaction{1}.rxns = sprintf('%s_GLER_COPI_formation_sec_Arf1p_Gea1p_Gea2p_Rer1p_Erd2p_Cop1p_Sec26p_Sec27p_Sec21p_Ret2p_Sec28p_Ret3p_complex',peptide_org);
reaction{2}.rxns = sprintf('%s_GLER_COPI_uncoating_and_fission_sec_Rer1p_Ret2p_Cop1p_Sec27p_Sec21p_Bet1p_complex',peptide_org);


reaction{1}.rxnNames = sprintf('%s_GLER_COPI_formation_Arf1p_Gea1p_Gea2p_Rer1p_Erd2p_Cop1p_Sec26p_Sec27p_Sec21p_Ret2p_Sec28p_Ret3p',peptide);
reaction{2}.rxnNames = sprintf('%s_GLER_COPI_uncoating_and_fission_Rer1p_Ret2p_Cop1p_Sec27p_Sec21p_Bet1p',peptide);


reaction{1}.eq = sprintf('%s[g] + 2 GTP[c] + 2 H2O[c] => %s_COPI_G1[c] + 2 GDP[c] + 2 H+[c] + 2 phosphate[c]',peptide,peptide);
reaction{2}.eq=  sprintf('%s_COPI_G1[c] => %s_mature[er]',peptide,peptide_org);
  if strcmp(compartment,'erm')== 1
    reaction{3}.rxns = sprintf('%s_GLER3_Final_demand',peptide_org);
    reaction{3}.rxnNames = sprintf('%s_GLER3_final demand',peptide);
    reaction{3}.eq=  sprintf('%s_mature[er] => %s_folding[erm]',peptide_org,peptide_org);
  else
    reaction{3}.rxns = sprintf('%s_GLER3_Final_demand',peptide_org);
    reaction{3}.rxnNames = sprintf('%s_GLER3_final demand',peptide);
    reaction{3}.eq=  sprintf('%s_mature[er] => %s_folding[er]',peptide_org,peptide_org);
    
   end   
end