function reaction = Post_Translation_translocation_tail(peptide)
%This function makes five reactions that translocates a peptide using post
%translation
%The function returns these five reaction

reaction{1}.rxns = sprintf('%s_Post_translation_PSTB_sec_Sgt2_Get4_Get5_complex',peptide);
reaction{2}.rxns = sprintf('%s_Post_translation_PSTB_sec_Get3_complex',peptide);
reaction{3}.rxns = sprintf('%s_Post_translation_PSTB_sec_Get1_Get2_complex',peptide);


reaction{1}.rxnNames = sprintf('%s_Post_translation_PSTB_sec_Sgt2_Get4_Get5_complex',peptide);
reaction{2}.rxnNames =sprintf('%s_Post_translation_PSTB_sec_Get3_complex',peptide);
reaction{3}.rxnNames = sprintf('%s_Post_translation_PSTB_sec_Get1_Get2_complex',peptide);


reaction{1}.eq = sprintf('%s_peptide[c] => %s_tanslocate_1[c]',peptide,peptide);
reaction{2}.eq= sprintf('%s_tanslocate_1[c] + ATP[c] + H2O[c] => %s_tanslocate_2[c] + ADP[c] + H+[c] + phosphate[c]',peptide,peptide);
reaction{3}.eq= sprintf('%s_tanslocate_2[c] => %s[er]',peptide,peptide);
end
