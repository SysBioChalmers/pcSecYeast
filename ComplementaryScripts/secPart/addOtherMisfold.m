function [newModel,peptide_name,rxns] = addOtherMisfold(model,peptide,peptide_org,Length,NG,OG,DSB,GPI,onlyrxns)
rxns = [];
Length = round(Length/40);

if NG == 0 && GPI == 0
    if DSB > 0
        reaction{1}.rxns = sprintf('%s_misfold_er_FL6_sec_Kar2p_complex',peptide_org);
        reaction{2}.rxns = sprintf('%s_ERADegradation2A_sec_Pdi1p_complex',peptide);
        reaction{1}.rxnNames = sprintf('%s_FL5_misfolding_sec_Kar2p_complex',peptide);
        reaction{2}.rxnNames = sprintf('%s_ERADegradation2A_sec_Pdi1p_complex',peptide);
        reaction{1}.eq = sprintf('%s[er] + %d ATP[er] + %d H2O[er] => %s_misf[er] + %d ADP[c] + %d H+[c] + %d phosphate[c]',peptide,Length,Length,peptide,Length,Length,Length);

        %reaction{2}.eq=  sprintf('%s_misf[er] + %d FADH2[er] => %s_misf_G1[er] + %d FAD[er] + %d H+[er]',peptide,2*DSB,peptide,2*DSB,2*DSB);
        reaction{2}.eq=  sprintf('%s_misf[er] + %d glutathione[er] => %s_misf_G1[er] + %d glutathione disulfide[er] + %d H+[er]',peptide,2*DSB,peptide,DSB,2*DSB);

    elseif DSB == 0
        reaction{1}.rxns = sprintf('%s_misfold_er_FL5_sec_Kar2p_complex',peptide_org);
        reaction{2}.rxns = sprintf('%s_ERADegradation2B',peptide);
        reaction{1}.rxnNames = sprintf('%s_FL5_misfolding',peptide);
        reaction{2}.rxnNames = sprintf('%ss_ERADegradation2A',peptide);
        reaction{1}.eq = sprintf('%s[er] + %d ATP[er] + %d H2O[er] => %s_misf[er] + %d ADP[c] + %d H+[c] + %d phosphate[c]',peptide,Length,Length,peptide,Length,Length,Length);

        reaction{2}.eq=  sprintf('%s_misf[er] => %s_misf_G1[er]',peptide,peptide);
    end

 %refolding
    reaction{3}.rxns = sprintf('%s_refolding_er',peptide);
    reaction{3}.rxnNames = sprintf('%s_refolding_er',peptide);
    reaction{3}.eq = sprintf('%s_misf[er] => %s[er]',peptide,peptide);

    reaction{4}.rxns = sprintf('%s_ERADegradation3_sec_Yos9p_Hrd3p_Hrd1p_Usa1p_Der1p_complex',peptide);
    reaction{5}.rxns = sprintf('%s_ERADegradatio4_sec_Sbh1p_Sss1p_Ssh1p_Ubc7p_Cue1p_Ubx2p_Cdc48p_Ufd1p_Npl4p_Hrd3p_Hrd1p_Usa1p_Der1p_complex',peptide);

    reaction{4}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd3p_Hrd1p_Usa1p_Der1p',peptide);
    reaction{5}.rxnNames = sprintf('%s_ERADegradation4',peptide);

    reaction{4}.eq = sprintf('%s_misf_G1[er] => %s_misf_G2[er]',peptide,peptide);
    reaction{5}.eq = sprintf('%s_misf_G2[er] + 4 Ubiquitin_for_Transfer[c] => %s_misf_G3[c] + 4 Ubiquitin[c]',peptide,peptide);
   if OG > 0
        reaction{6}.rxns = sprintf('%s_ERADegradation5_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
        reaction{6}.rxnNames = sprintf('%s_ERADegradation5A_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
        reaction{6}.eq = sprintf('%s_misf_G3[c] => %s_misfolding[c] + %d D-mannose[er]',peptide,peptide_org,OG);
   else
        reaction{6}.rxns = sprintf('%s_ERADegradation5',peptide);
        reaction{6}.rxnNames = sprintf('%s_ERADegradation5B',peptide);
        reaction{6}.eq = sprintf('%s_misf_G3[c] => %s_misfolding[c]',peptide,peptide_org);
   end
    n = length(reaction) + 1;
    reaction{n}.rxns = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.rxnNames = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.eq = sprintf('%s_misfolding[c] => %s_subunit[c]',peptide_org,peptide_org);
    %misfolding dilution
%     reaction{n+1}.rxns = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.eq = sprintf('%s_misfolding[c] => ',peptide_org);
   %misfolding retension cycle
    reaction{n+1}.rxns = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide);
    reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide);
    reaction{n+1}.eq = sprintf('%s_misf[er] +  + %d glutathione[er] + %d oxygen[er] => %d hydrogen peroxide[er] + %d glutathione disulfide[er] + %s_misfolding_acc[er]',peptide,10,5,5,5,peptide);
    %misfolding dilution
    reaction{n+2}.rxns = sprintf('%s_dilution_misfolding_er',peptide_org);
    reaction{n+2}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide);
    reaction{n+2}.eq = sprintf('%s_misfolding_acc[er] => ',peptide);

    for i=1:length(reaction)
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s',peptide);
else
    newModel = model;
    peptide_name = peptide;
end
end
