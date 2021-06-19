function reaction = co_Translation_translocation(peptide)
%This function makes five reactions that translocates a peptide hasing SP
%signal.
%The function returns these five reaction
%Ibrahim Elsemman 07-12-2018


reaction{1}.rxns = sprintf('%s_co_translation_TC_sec_SRPC_complex',peptide);
reaction{2}.rxns = sprintf('%s_co_translation_TC_sec_SRC_complex',peptide);
reaction{3}.rxns = sprintf('%s_co_translation_TC_sec_SEC61C_complex',peptide);
reaction{4}.rxns = sprintf('%s_co_translation_TC_sec_SSH1C_complex',peptide);
reaction{5}.rxns = sprintf('%s_co_translation_TC_sec_SPC_complex',peptide);
reaction{6}.rxns = sprintf('%s_export_sp_to_c',peptide);

reaction{1}.rxnNames = sprintf('%s signal peptide recognition',peptide);
reaction{2}.rxnNames = sprintf('%s signal peptide recognition',peptide);
reaction{3}.rxnNames = sprintf('Biding of %s-SRPC-SRC to the translocator SEC61C SRCPC and SRC dissociation',peptide);
reaction{4}.rxnNames = sprintf('Biding of %s-SRPC-SRC to the translocator SEC61C alternative SRCPC and SRC dissociation',peptide);
reaction{5}.rxnNames = sprintf('%s Signal peptidase',peptide);
reaction{6}.rxnNames = sprintf('%s_export sp to cytosol',peptide);

reaction{1}.eq = sprintf('%s_peptide[c] => %s_tanslocate_1[c]',peptide,peptide);
reaction{2}.eq= sprintf('%s_tanslocate_1[c] => %s_tanslocate_2[c]',peptide,peptide);
reaction{3}.eq= sprintf('%s_tanslocate_2[c] + 2 GTP[c] + 2 H2O[c] => %s_tanslocate_3[c] + 2 GDP[c] + 2 phosphate[c] + H+[c]',peptide,peptide);
reaction{4}.eq= sprintf('%s_tanslocate_2[c] + 2 GTP[c] + 2 H2O[c] => %s_tanslocate_3[c] + 2 GDP[c] + 2 phosphate[c] + H+[c]',peptide,peptide);
reaction{5}.eq= sprintf('%s_tanslocate_3[c] + H2O[c] => %s[er] + %s_sp[er]',peptide,peptide,peptide);
reaction{6}.eq = sprintf('%s_sp[er] => %s_sp[c]',peptide,peptide);
