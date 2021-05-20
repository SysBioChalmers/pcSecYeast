function [newModel,peptide_name,rxns] = addGPIMisfold(model,peptide,peptide_org,Length,NG,OG,DSB,GPI,onlyrxns)
%since the misfolding for NG starts before the trimming of M
rxns = [];
Length = round(Length/40);
if GPI > 0 && NG > 0
    S = regexp(peptide, '_M8', 'split');
    peptide = [S{1},'_M9'];
    man = 8*NG + OG;
    nac = 2*NG;
end

if GPI > 0
    if GPI > 1
        GPI = 1;
    end
    reaction{1}.rxns = sprintf('%s_misfold_er_GPI_I_sec_Kar2p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_refolding_er',peptide);
    reaction{3}.rxns = sprintf('%s_GPImisfoldII',peptide);
    
    
    reaction{1}.rxnNames = sprintf('%s_GPImisfoldI_sec_Kar2p_complex',peptide);
    reaction{2}.rxnNames = sprintf('%s_refolding_er',peptide);
    reaction{3}.rxnNames = sprintf('%s_GPImisfoldII',peptide);
    reaction{1}.eq = sprintf('%s[er] + %d ATP[er] + %d H2O[er] => %s_misf[er] + %d ADP[c] + %d H+[c] + %d phosphate[c]',peptide,Length,Length,peptide,Length,Length,Length);

    reaction{2}.eq = sprintf('%s_misf[er] => %s[er]',peptide,peptide);
    reaction{3}.eq = sprintf('%s_misf[er] + %.15f H2O[er] => %s_misf_G1[er] + %.15f 6-O-2-O-alpha-D-mannosyl-(1-2)-{alpha-D-mannosyl-2-O-((2-aminoethyl)phosphoryl)-(1-2)-alpha-D-mannosyl-(1-6)-2-O-((2-aminoethyl)phosphoryl)-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-inositol-P-ceramide C (C26)[er]',peptide,GPI,peptide,GPI);
    if DSB > 0
        reaction{4}.rxns = sprintf('%s_GPImisfoldIIIA_sec_Pdi1p_complex',peptide);
        reaction{4}.rxnNames = sprintf('%s_GPImisfoldIIIA_Pdi1p_complex',peptide);
        %reaction{4}.eq = sprintf('%s_misf_G1[er] + %d FADH2[er] => %s_misf_G2[er] + %d FAD[er] + %d H+[er]',peptide,2*DSB,peptide,2*DSB,2*DSB);
        reaction{4}.eq=  sprintf('%s_misf_G1[er] + %d glutathione[er] => %s_misf_G2[er] + %d glutathione disulfide[er] + %d H+[er]',peptide,2*DSB,peptide,DSB,2*DSB);

    elseif DSB == 0
        reaction{4}.rxns = sprintf('%s_GPImisfoldIIIB',peptide);
        reaction{4}.rxnNames = sprintf('%s_GPImisfoldIIIB',peptide);
        reaction{4}.eq = sprintf('%s_misf_G1[er] => %s_misf_G2[er]',peptide,peptide);
    end
    if NG > 0
        reaction{5}.rxns = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
        reaction{6}.rxns = sprintf('%s_GPImisfoldV_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
        reaction{5}.rxnNames = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
        reaction{6}.rxnNames = sprintf('%s_GPImisfoldV_sec_Dsk2p_Rad23p_Png1p_Uba1p_complex',peptide);
        
        reaction{5}.eq = sprintf('%s_misf_G2[er] + GTP[er] + H2O[er] => %s_misf_G3[c] + GDP[er] + H+[er] + phosphate[er]',peptide,peptide);
        reaction{6}.eq = sprintf('%s_misf_G3[c] => %s_misfolding[c] + %d D-mannose[er] + %d N-acetyl-alpha-D-glucosamine 1-phosphate[c]',peptide,peptide_org,man,nac);
    elseif NG == 0
        if OG > 0
            reaction{5}.rxns = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
            reaction{5}.rxnNames = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
            reaction{5}.eq = sprintf('%s_misf_G2[er] + GTP[er] + H2O[er] => %s_misfolding[c] + GDP[er] + H+[er] + phosphate[er] + %d D-mannose[er]',peptide,peptide_org,OG);
        else
            reaction{5}.rxns = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
            reaction{5}.rxnNames = sprintf('%s_GPImisfoldIIII_sec_Sec12p_Sar1p_Sec23p_Sec24p_Emp24p_Erp1p_Erp2p_Erv25p_complex',peptide);
            reaction{5}.eq = sprintf('%s_misf_G2[er] + GTP[er] + H2O[er] => %s_misfolding[c] + GDP[er] + H+[er] + phosphate[er]',peptide,peptide_org);
        end
        
    end
    n = length(reaction) + 1;
    reaction{n}.rxns = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.rxnNames = sprintf('%s_degradation_misfolding_c',peptide_org);
    reaction{n}.eq = sprintf('%s_misfolding[c] => %s_subunit[c]',peptide_org,peptide_org);
%     %misfolding dilution
%     reaction{n+1}.rxns = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.rxnNames = sprintf('%s_dilution_misfolding_c',peptide_org);
%     reaction{n+1}.eq = sprintf('%s_misfolding[c] => ',peptide_org);
%     %misfolding retension cycle
    reaction{n+1}.rxns = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide);
    reaction{n+1}.rxnNames = sprintf('%s_cycle_accumulation_sec_pdi1p_ero1p_complex',peptide);
    reaction{n+1}.eq = sprintf('%s_misf[er] + %d glutathione[er] + %d oxygen[er] => %d hydrogen peroxide[er] + %d glutathione disulfide[er] + %s_misfolding_acc[er]',peptide,10,5,5,5,peptide);
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
    peptide_name = peptide;
else
    newModel = model;
    peptide_name = peptide;
end

end
