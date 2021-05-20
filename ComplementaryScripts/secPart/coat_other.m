function [newModel,peptide_name,rxns] = coat_other(model,peptide,GPI,trans,onlyrxns)
rxns = [];
if trans == 0 && GPI == 0
    reaction{1}.rxns = sprintf('%s_COPII_normal_ERGLD_sec_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p_complex',peptide);
    reaction{2}.rxns = sprintf('%s_COPII_common_ERGLA_sec_Sec13p_Sec31p_Sec16p_Sed4p_Sec5p_Sec17p_complex',peptide);
    reaction{3}.rxns = sprintf('%s_COPII_common_ERGLA_sec_Ypt1p_Uso1p_bug1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p_complex',peptide);
    
    reaction{1}.rxnNames =  sprintf('%s_COPII_normal_ERGLD_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p Pre budding complex forming for soluble proteins',peptide);
    reaction{2}.rxnNames =  sprintf('%s_COPII_common_ERGLA_Sec12p_Sar1p_Sec23p_Sec24p_Erv29p COPII formation',peptide);
    reaction{3}.rxnNames = sprintf('%s_COPII_common_ERGLA_Ypt1p_Uso1p_bug1p_Bet3p_Bet5p_Trs20p_Trs23p_Trs31p_Trs33p COPII fusion',peptide);
    
    reaction{1}.eq = sprintf('%s[er] + GTP[er] + H2O[er] => %s_COP_coated[er] + GDP[er] + H+[er] + phosphate[er]',peptide,peptide);
    reaction{2}.eq = sprintf('%s_COP_coated[er] => %s_COP_coated[c]',peptide,peptide);
    reaction{3}.eq = sprintf('%s_COP_coated[c] + GTP[c] + H2O[c] => %s[g] + GDP[c] + H+[c] + phosphate[c]',peptide,peptide);
    
    for i=1:3
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