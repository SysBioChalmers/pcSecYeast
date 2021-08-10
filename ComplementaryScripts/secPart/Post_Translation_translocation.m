function reaction = Post_Translation_translocation(peptide,Length)
%This function makes five reactions that translocates a peptide using post
%translation
%The function returns these five reaction

reaction{1}.rxns = sprintf('%s_Post_translation_PSTA_sec_RAC_complex',peptide);
reaction{2}.rxns = sprintf('%s_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex',peptide);
reaction{3}.rxns = sprintf('%s_Post_translation_PSTA_sec_SEC61SEC63C_complex',peptide);
reaction{4}.rxns = sprintf('%s_Post_translation_PSTA_sec_BIP_NEFS_complex',peptide);

reaction{1}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_RAC_complex',peptide);
reaction{2}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_Ssa1_Ydj1_Snl1_complex',peptide);
reaction{3}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_SEC61SEC63C_complex',peptide);
reaction{4}.rxnNames = sprintf('%s_Post_translation_PSTA_sec_BIP_NEFS_complex',peptide);


reaction{1}.eq = sprintf('%s_peptide[c] => %s_tanslocate_1[c]',peptide,peptide);
reaction{2}.eq= sprintf('%s_tanslocate_1[c] + ATP[c] + H2O[c] => %s_tanslocate_2[c] + ADP[c] + H+[c] + phosphate[c]',peptide,peptide);
reaction{3}.eq= sprintf('%s_tanslocate_2[c] => %s_tanslocate_3[c]',peptide,peptide);
reaction{4}.eq= sprintf('%s_tanslocate_3[c] + %.15f ATP[c] + %.15f H2O[c] => %s[er] + %.15f ADP[c] + %.15f H+[c] + %d phosphate[c]',peptide,Length,Length,peptide,Length,Length,Length);
end


