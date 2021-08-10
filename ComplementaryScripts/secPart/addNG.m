function [newModel,peptide_name,rxns] = addNG(model,peptide,peptide_org,Length_total,NG,onlyrxns)
rxns = [];
if NG >0
    Length = Length_total/40;
    
    
    reaction{1}.rxns = sprintf('%s_ERNG_NG_sec_OSTC_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_ERNG_FLI_NG_sec_Cwh41p_complex',peptide_org);
    reaction{3}.rxns = sprintf('%s_ERNG_FLII_NG_sec_Rot2p_complex',peptide_org);
    reaction{4}.rxns = sprintf('%s_ERNG_FLIII_NG_sec_Rot2p_complex',peptide_org);
    reaction{5}.rxns = sprintf('%s_ERNG_FLIV_NG_sec_Mns1p_complex',peptide_org);

    reaction{1}.rxnNames = sprintf('%s_ERNG_NG_OSTC_complex ER N-glycosylation',peptide);
    reaction{2}.rxnNames = sprintf('%s_ERNG_FLI_NG_Cwh41p_complex ER Glycan trimming I',peptide);
    reaction{3}.rxnNames = sprintf('%s_ERNG_FLII_NG_Rot2p_complex ER Glycan trimming II',peptide);
    reaction{4}.rxnNames = sprintf('%s_ERNG_FLIII_NG_Rot2p_complex ER Glycan trimming III',peptide);
    reaction{5}.rxnNames = sprintf('%s_ERNG_FLIV_NG_Mns1p_complex ER demanosylation I',peptide);
    
    reaction{1}.eq = sprintf('%s[er] + %d Glucose(3)Mannose(9)GlucoseNAc(2)-PP-dolichol[er] => %s_G3M9[er] + %d dolichyl phosphate[er]',peptide,NG,peptide,NG);
    reaction{2} .eq=  sprintf('%s_G3M9[er] + %d H2O[er] => %s_G2M9[er] + %d D-glucose[er]',peptide,NG,peptide,NG);
    reaction{3} .eq=  sprintf('%s_G2M9[er] + %d H2O[er] => %s_G1M9[er] + %d D-glucose[er]',peptide,NG,peptide,NG);
    reaction{4} .eq=  sprintf('%s_G1M9[er] + %d H2O[er] => %s_M9[er] + %d D-glucose[er]',peptide,NG,peptide,NG);
    reaction{5} .eq=  sprintf('%s_M9[er] + %d H2O[er] => %s_M8[er] + %d D-mannose[er]',peptide,NG,peptide,NG);

    for i=1:5
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s_M8',peptide); %here the peptide name doesn't change to the final metaboloites due to the reason to link the nexr NGmisfolding reactions
else
    newModel = model;
    peptide_name = peptide;
    
end