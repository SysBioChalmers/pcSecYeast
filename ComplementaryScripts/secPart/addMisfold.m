function [newModel,peptide_name,rxns] = addMisfold(model,peptide,peptide_org,Length,NG,OG,DSB,GPI,Trans,compartment,onlyrxns)

rxns = [];
Length = round(Length/40);
if NG > 0
    S = regexp(peptide, '_M8', 'split');
    peptide = S{1};
    peptide = [peptide,'_M9'];
    man = 7*NG + OG;
    nac = 2*NG;
elseif NG== 0 && OG > 0
     man = OG;
end

if GPI > 1
    GPI = 1;
end
    
  reaction{1}.rxns = sprintf('%s_misfold_ERAD_sec_Kar2p_complex',peptide_org); % misfold is the reserve word for constrain misfolding ratio
  reaction{1}.rxnNames = sprintf('%s_ERAD_misfold_sec_Kar2p_complex',peptide);
  reaction{1}.eq = sprintf('%s[er] + %d ATP[er] + %d H2O[er] => %s_misf[er] + %d ADP[er] + %d H+[er] + %d phosphate[er]',peptide,Length,Length,peptide,Length,Length,Length);

if DSB > 0
        reaction{2}.rxns = sprintf('%s_ERAD2A_sec_Pdi1p_complex',peptide_org);
        reaction{2}.rxnNames = sprintf('%s_ERAD2A_sec_Pdi1p_complex',peptide);
        reaction{2}.eq=  sprintf('%s_misf[er] + %d glutathione[er] => %s_misf_G1[er] + %d glutathione disulfide[er] + %d H+[er]',peptide,2*DSB,peptide,DSB,2*DSB);
else
       reaction{2}.rxns = sprintf('%s_ERAD2B',peptide_org);
        reaction{2}.rxnNames = sprintf('%s_ERAD2B',peptide);
        reaction{2}.eq = sprintf('%s_misf[er] => %s_misf_G1[er]',peptide,peptide);
end

if NG > 0 
    reaction{3}.rxns = sprintf('%s_ERAD3A_sec_Mns1p_complex',peptide_org);
    reaction{3}.rxnNames = sprintf('%s_ER demanosylation_s_ERAD3A_sec_Mns1p_complex',peptide);
    reaction{3}.eq =  sprintf('%s_misf_G1[er] + %d H2O[er] => %s_misf_G2[er] + %d D-mannose[er]',peptide,NG,peptide,NG);
    reaction{4}.rxns = sprintf('%s_ERAD4A_sec_Mnl1p_Pdi1p_complex',peptide_org);
    reaction{4}.rxnNames = sprintf('%s_ERNG_FL_NG_Htm1p_Pdi1P complex ER demanosylation',peptide);
    reaction{4}.eq = sprintf('%s_misf_G2[er] + %d H2O[er] => %s_misf_G3[er] + %d D-mannose[er]',peptide,NG,peptide,NG);

else
    reaction{3}.rxns = sprintf('%s_ERAD3B',peptide_org);
    reaction{3}.rxnNames = sprintf('%s_ER no real function3',peptide);
    reaction{3}.eq=  sprintf('%s_misf_G1[er] => %s_misf_G2[er]',peptide,peptide);
    reaction{4}.rxns = sprintf('%s_ERAD4B',peptide_org);
    reaction{4}.rxnNames = sprintf('%s_ER no real function4',peptide);
    reaction{4}.eq = sprintf('%s_misf_G2[er] => %s_misf_G3[er]',peptide,peptide);
end
if GPI > 0 
    reaction{5}.rxns = sprintf('%s_ERAD5A',peptide_org);
    reaction{5}.rxnNames = sprintf('%s_GPImisfoldII',peptide);
    reaction{5}.eq = sprintf('%s_misf_G3[er] + %.15f H2O[er] => %s_misf_G4[er] + %.15f 6-O-2-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-inositol-P-ceramide C (C26)[er]',peptide,GPI,peptide,GPI);
else
   reaction{5}.rxns = sprintf('%s_ERAD5B',peptide_org);
    reaction{5}.rxnNames = sprintf('%s_GPImisfoldII',peptide);
    reaction{5}.eq = sprintf('%s_misf_G3[er] => %s_misf_G4[er]',peptide,peptide);
end
   if Trans > 0 && strcmp(compartment,'c') && DSB == 0 && NG == 0 % ERADC
        reaction{6}.rxns = sprintf('%s_ERADC_sec_Cue1p_Ubc6p_Ubc7p_Doa10p_complex',peptide_org);
        reaction{6}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein',peptide);
        reaction{6}.eq = sprintf('%s_misf_G4[er] => %s_misf_G5[er]',peptide,peptide);
                    
        reaction{7}.rxns = sprintf('%s_ERADC_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex',peptide_org);
        reaction{7}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptide);
        reaction{7}.eq = sprintf('%s_misf_G5[er] + 8 Ubiquitin_for_Transfer[c] => %s_misf_G6[c] + 8 Ubiquitin[c]',peptide,peptide);
   elseif Trans > 0 && (strcmp(compartment,'er')||strcmp(compartment,'erm')) %ERADM
            
        reaction{6}.rxns = sprintf('%s_ERADM_sec_Cue1p_Ubc6p_Ubc7p_Hrd1p_Hrd3p_Der1p_complex',peptide_org);
        reaction{6}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptide);
        reaction{6}.eq = sprintf('%s_misf_G4[er] => %s_misf_G5[er]',peptide,peptide);
                    
        reaction{7}.rxns = sprintf('%s_ERADM_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex',peptide_org);
        reaction{7}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptide);
        reaction{7}.eq = sprintf('%s_misf_G5[er] + 8 Ubiquitin_for_Transfer[c] => %s_misf_G6[c] + 8 Ubiquitin[c]',peptide,peptide);
   else % ERADL
             
        reaction{6}.rxns = sprintf('%s_ERADL_sec_Cue1p_Ubc6p_Ubc7p_Yos9p_Hrd1p_Hrd3p_Der1p_Usa1p_complex',peptide_org);
        reaction{6}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd3p_Hrd1p_Usa1p_Der1p',peptide);
        reaction{6}.eq = sprintf('%s_misf_G4[er] => %s_misf_G5[er]',peptide,peptide);
        reaction{7}.rxns = sprintf('%s_ERADL_sec_Sbh1p_Sss1p_Ssh1p_Cdc48p_Ubx2p_Ufd1p_Npl4p_complex',peptide_org);
        reaction{7}.rxnNames = sprintf('%s_Survaliancecplx formation of misfoled protein_Yos9p_Hrd1p_Hrd3p',peptide);
        reaction{7}.eq = sprintf('%s_misf_G5[er] + 8 Ubiquitin_for_Transfer[c] => %s_misf_G6[c] + 8 Ubiquitin[c]',peptide,peptide);
   end
   if NG > 0
        reaction{8}.rxns = sprintf('%s_ERAD7A_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide_org);
        reaction{8}.rxnNames = sprintf('%s_ERAD_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
    
        reaction{8}.eq = sprintf('%s_misf_G6[c] => %s_misfolding[c] + %d D-mannose[er] + %d N-acetyl-alpha-D-glucosamine 1-phosphate[c]',peptide,peptide_org,man,nac);
    elseif NG == 0
        if OG > 0
            reaction{8}.rxns = sprintf('%s_ERAD7B_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide_org);
            reaction{8}.rxnNames = sprintf('%s_ERAD_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
            reaction{8}.eq = sprintf('%s_misf_G6[c] => %s_misfolding[c] + %d D-mannose[er]',peptide,peptide_org,man);

        else
            reaction{8}.rxns = sprintf('%s_ERAD7C_sec_Dsk2p_Rad23p_Uba1p_complex',peptide_org);
            reaction{8}.rxnNames = sprintf('%s_ERAD_sec_Dsk2p_Rad23p_Uba1p_complex',peptide);
            reaction{8}.eq = sprintf('%s_misf_G6[c] => %s_misfolding[c]',peptide,peptide_org);
        end
   end
    reaction{9}.rxns = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{9}.rxnNames = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{9}.eq = sprintf('%s_misfolding[c] => %s_subunit[c]',peptide_org,peptide_org);
%     %misfolding dilution
%     reaction{n+1}.rxns = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.eq = sprintf('%s_misfolding[c] => ',peptide_org);

     %refolding
    reaction{10}.rxns = sprintf('%s_refolding_er',peptide_org);
    reaction{10}.rxnNames = sprintf('%s_refolding_er',peptide);
    reaction{10}.eq = sprintf('%s_misf[er] => %s[er]',peptide,peptide);
        n = length(reaction);

    if DSB > 0
        reaction{n+1}.rxns = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide_org);
        reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide);
        reaction{n+1}.eq = sprintf('%s_misf[er] + %d glutathione[er] + %d oxygen[er] => %d hydrogen peroxide[er] + %d glutathione disulfide[er] + %s_misf2[er]',peptide,20*DSB,10*DSB,10*DSB,10*DSB,peptide);
    else
        reaction{n+1}.rxns = sprintf('%s_cycle_accumulation',peptide_org);
        reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation',peptide);
        reaction{n+1}.eq = sprintf('%s_misf[er] => %s_misf2[er]',peptide,peptide);
    end
    reaction{n+2}.rxns = sprintf('%s_cycle_accumulation_sec_acc_Kar2p_complex',peptide_org);
    reaction{n+2}.rxnNames = sprintf('%s_cycle_accumulation_sec_acc_kar2p_complex',peptide);
    reaction{n+2}.eq = sprintf('%s_misf2[er] + %d ATP[er] + %d H2O[er] => %s_misfolding_acc[er] + %d ADP[er] + %d H+[er] + %d phosphate[er]',peptide,10*Length,10*Length,peptide,10*Length,10*Length,10*Length);

    %misfolding dilution
    reaction{n+3}.rxns = sprintf('%s_dilution_misfolding_er',peptide_org);
    reaction{n+3}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide);
    reaction{n+3}.eq = sprintf('%s_misfolding_acc[er] => ',peptide);
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
    if NG > 0
    S = regexp(peptide_name, '_M9', 'split');
    peptide_name = [S{1},'_M8'];
    end

end
