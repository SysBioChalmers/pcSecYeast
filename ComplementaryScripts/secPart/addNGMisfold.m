function [newModel,peptide_name,rxns] = addNGMisfold(model,peptide,peptide_org,Length,NG,OG,DSB,GPI,Trans,compartment,onlyrxns)

rxns = [];
Length = round(Length/40);
if NG > 0 && GPI == 0
    S = regexp(peptide, '_M8', 'split');
    peptidenew = S{1};
    man = 8*NG + OG;
    nac = 2*NG;
    if DSB > 0
    reaction{1}.rxns = sprintf('%s_misfold_er_FL5_sec_Kar2p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_ERAD1A_sec_Mns1p_Pdi1p_complex',peptidenew);
    reaction{1}.rxnNames = sprintf('%s_FL5_misfold_sec_Kar2p_complex',peptidenew);
    reaction{2}.rxnNames = sprintf('%s_ER demanosylation_II_sec_Kar2p_Mns1p_Pdi1p_complex',peptidenew);
    reaction{1}.eq = sprintf('%s_M9[er] + %d ATP[er] + %d H2O[er] => %s_M9_misf[er] + %d ADP[er] + %d H+[er] + %d phosphate[er]',peptidenew,Length,Length,peptidenew,Length,Length,Length);

    %reaction{2}.eq=  sprintf('%s_M9_misf[er] + %d FADH2[er] + %d H2O[er] => %s_M8_misf[er] + %d D-mannose[er] + %d FAD[er] + %d H+[er]',peptidenew,2*DSB,NG,peptidenew,NG,2*DSB,2*DSB);
    reaction{2}.eq=  sprintf('%s_M9_misf[er] + %d glutathione[er] + %d H2O[er] => %s_M8_misf[er] + %d D-mannose[er] + %d glutathione disulfide[er] + %d H+[er]',peptidenew,2*DSB,NG,peptidenew,NG,DSB,2*DSB);

    elseif DSB == 0
    reaction{1}.rxns = sprintf('%s_misfold_er_FL5_sec_Kar2p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_ERAD1B_sec_Mns1p_complex',peptidenew);
    reaction{1}.rxnNames = sprintf('%s_FL5_misfold_sec_Kar2p_complex',peptidenew);
    reaction{2}.rxnNames = sprintf('%s_ERNG_FL_NG_ ER Glycan trimming II',peptidenew);
    reaction{1}.eq = sprintf('%s_M9[er] + %d ATP[er] + %d H2O[er] => %s_M9_misf[er] + %d ADP[er] + %d H+[er] + %d phosphate[er]',peptidenew,Length,Length,peptidenew,Length,Length,Length);

    reaction{2}.eq = sprintf('%s_M9_misf[er] + %d H2O[er] => %s_M8_misf[er] + %d D-mannose[er]',peptidenew,NG,peptidenew,NG);
    end

    reaction{3}.rxns = sprintf('%s_ERADII_sec_Mnl1p_complex',peptidenew);
    reaction{3}.rxnNames = sprintf('%s_ERNG_FL_NG_Htm1p ER demanosylation',peptidenew);
    reaction{3}.eq = sprintf('%s_M8_misf[er] + %d H2O[er] => %s_M7_misf[er] + %d D-mannose[er]',peptidenew,NG,peptidenew,NG);
    %refolding
    reaction{4}.rxns = sprintf('%s_refolding_er',peptidenew);
    reaction{4}.rxnNames = sprintf('%s_refolding_er',peptidenew);
    reaction{4}.eq = sprintf('%s_M9_misf[er] => %s_M9[er]',peptidenew,peptidenew);

    if strcmp(compartment,'s') || strcmp(compartment,'ce')
        reaction{5}.rxns = sprintf('%s_ERADLIII_sec_Yos9p_Hrd3p_Hrd1p_Usa1p_Der1p_complex',peptidenew);
        reaction{6}.rxns = sprintf('%s_ERADLIV_sec_Usa1p_Der1p_Ubx2p_Hrd1p_Ssh1p_Sbh2p_Sss1p_complex',peptidenew);
        reaction{7}.rxns = sprintf('%s_ERADLV_sec_Sbh1p_Sss1p_Ssh1p_Ubc7p_Cue1p_Ubx2p_Cdc48p_Ufd1p_Npl4p_Hrd3p_Hrd1p_Usa1p_Der1p_complex',peptidenew);
        reaction{8}.rxns = sprintf('%s_ERADLVI_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptidenew);

        reaction{5}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd3p_Hrd1p_Usa1p_Der1p',peptidenew);
        reaction{6}.rxnNames = sprintf('%s_Co retrotranslocation ubiquitin addition to the misfold prtotein_Sbh1p_Sss1p_Ssh1p_Ubc7p_Cue1p_Ubx2p_Cdc48p_Ufd1p_Npl4p_Hrd3p_Hrd1p_Usa1p_Der1p',peptidenew);
        reaction{7}.rxnNames = sprintf('%s_ERADLV',peptidenew);
        reaction{8}.rxnNames = sprintf('%s_ERADLVI_Dsk2p_Rad23p_Png1p_Uba1p',peptidenew);


        reaction{5}.eq = sprintf('%s_M7_misf[er] => %s_M7_misf_G1[er]',peptidenew,peptidenew);
        reaction{6}.eq = sprintf('%s_M7_misf_G1[er] => %s_M7_misf_G2[er]',peptidenew,peptidenew);
        reaction{7}.eq = sprintf('%s_M7_misf_G2[er] + 8 Ubiquitin_for_Transfer[c] => %s_M7_misf_G3[c] + 8 Ubiquitin[c]',peptidenew,peptidenew);
        reaction{8}.eq = sprintf('%s_M7_misf_G3[c] => %s_misfolding[c] + %d D-mannose[er] + %d N-acetyl-alpha-D-glucosamine 1-phosphate[c]',peptidenew,peptidenew,man,nac);
    elseif Trans > 0
        reaction{5}.rxns = sprintf('%s_ERADMIII_sec_Yos9p_Hrd1p_Hrd3p_complex',peptidenew);
        reaction{6}.rxns = sprintf('%s_ERADMIV',peptidenew);
        reaction{7}.rxns = sprintf('%s_ERADMV_sec_Sbh1p_Sss1p_Ssh1p_Ubc7p_Cue1p_Ubx2p_Cdc48p_Ufd1p_Npl4p_Hrd3p_Hrd1p_Usa1p_Der1p_complex',peptidenew);
        reaction{8}.rxns = sprintf('%s_ERADMVI_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptidenew);

        reaction{5}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptidenew);
        reaction{6}.rxnNames = sprintf('%s_Co retrotranslocation ubiquitin addition to the misfold prtotein',peptidenew);
        reaction{7}.rxnNames = sprintf('%s_ERADMV_Sbh1p_Sss1p_Ssh1p_Ubc7p_Cue1p_Ubx2p_Cdc48p_Ufd1p_Npl4p_Hrd3p_Hrd1p_Usa1p_Der1p',peptidenew);
        reaction{8}.rxnNames = sprintf('%s_ERADMVI_Dsk2p_Rad23p_Png1p_Uba1p',peptidenew);


        reaction{5}.eq = sprintf('%s_M7_misf[er] => %s_M7_misf_G1[er]',peptidenew,peptidenew);
        reaction{6}.eq = sprintf('%s_M7_misf_G1[er] => %s_M7_misf_G2[er]',peptidenew,peptidenew);
        reaction{7}.eq = sprintf('%s_M7_misf_G2[er] + 8 Ubiquitin_for_Transfer[c] => %s_M7_misf_G3[c] + 8 Ubiquitin[c]',peptidenew,peptidenew);
        reaction{8}.eq = sprintf('%s_M7_misf_G3[c] => %s_misfolding[c] + %d D-mannose[er] + %d N-acetyl-alpha-D-glucosamine 1-phosphate[c]',peptidenew,peptide_org,man,nac);
    else
        reaction{5}.rxns = sprintf('%s_ERADCIII_sec_Cue1p_Ubc6p_Doa10p_Ubc7p_Cdc48p_Ubx2p_complex',peptidenew);
        reaction{6}.rxns = sprintf('%s_ERADCIV',peptidenew);
        reaction{7}.rxns = sprintf('%s_ERADCV_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptidenew);

        reaction{5}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptidenew);
        reaction{6}.rxnNames = sprintf('%s_ERADCIV',peptidenew);
        reaction{7}.rxnNames = sprintf('%s_ERADCV_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptidenew);


        reaction{5}.eq = sprintf('%s_M7_misf[er] => %s_M7_misf_G1[er]',peptidenew,peptidenew);
        reaction{6}.eq = sprintf('%s_M7_misf_G1[er] + 8 Ubiquitin_for_Transfer[c] => %s_M7_misf_G2[c] + 8 Ubiquitin[c]',peptidenew,peptidenew);
        reaction{7}.eq = sprintf('%s_M7_misf_G2[c] => %s_misfolding[c] + %d D-mannose[er] + %d N-acetyl-alpha-D-glucosamine 1-phosphate[c]',peptidenew,peptide_org,man,nac);
    end

    n = length(reaction) + 1;
    reaction{n}.rxns = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.rxnNames = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.eq = sprintf('%s_misfolding[c] => %s_subunit[c]',peptide_org,peptide_org);
%     %misfolding dilution
%     reaction{n+1}.rxns = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.eq = sprintf('%s_misfolding[c] => ',peptide_org);

    %misfolding retension cycle
    if DSB > 0
    reaction{n+1}.rxns = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptidenew);
    reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptidenew);
    reaction{n+1}.eq = sprintf('%s_M9_misf[er] + %d glutathione[er] + %d oxygen[er] => %d hydrogen peroxide[er] + %d glutathione disulfide[er] + %s_M9_misf2[er]',peptidenew,20*DSB,10*DSB,10*DSB,10*DSB,peptidenew);
    else
        reaction{n+1}.rxns = sprintf('%s_cycle_accumulation',peptidenew);
        reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation',peptidenew);
        reaction{n+1}.eq = sprintf('%s_M9_misf[er] => %s_M9_misf2[er]',peptidenew,peptidenew);
    end
    reaction{n+2}.rxns = sprintf('%s_cycle_accumulation_sec_acc_Kar2p_complex',peptidenew);
    reaction{n+2}.rxnNames = sprintf('%s_cycle_accumulation_sec_acc_Kar2p_complex',peptidenew);
    reaction{n+2}.eq = sprintf('%s_M9_misf2[er] + %d ATP[er] + %d H2O[er] => %s_misfolding_acc[er] + %d ADP[er] + %d H+[er] + %d phosphate[er]',peptidenew,10*Length,10*Length,peptidenew,10*Length,10*Length,10*Length);

    %misfolding dilution
    reaction{n+3}.rxns = sprintf('%s_dilution_misfolding_er',peptide_org);
    reaction{n+3}.rxnNames = sprintf('%s_dilution_misfolding_c',peptidenew);
    reaction{n+3}.eq = sprintf('%s_misfolding_acc[er] => ',peptidenew);

    for i=1:length(reaction)
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = peptide;
else
    newModel = model;
    peptide_name = peptide;

end
