function [newModel,peptide_name,rxns] = golgiProcessing_N(model,peptide,peptide_org,Length_total,NG,onlyrxns)
rxns = [];
if NG >0
    Length = Length_total/40;
    
    
    reaction{1}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex',peptide_org);
    reaction{2}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex',peptide_org);
    reaction{3}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex',peptide_org);
    reaction{4}.rxns = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn1p_Mnn2p_Mnn5p_complex',peptide_org);
    
    reaction{1}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_I_sec_Och1p_complex',peptide);
    reaction{2}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_MPOLI_complex',peptide);
    reaction{3}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_III_sec_MPoLII_complex',peptide);
    reaction{4}.rxnNames = sprintf('%s_GLNG_Golgi_N_linked_glycosylation_II_sec_Mnn1p_Mnn2p_Mnn5p_complex',peptide);
    
    
    reaction{1}.eq = sprintf('%s[g] + %d GDP-alpha-D-mannose[g] => %s_GNG_G1[g] + %d GDP[g]',peptide,NG,peptide,NG);
    reaction{2}.eq=  sprintf('%s_GNG_G1[g] + %d GDP-alpha-D-mannose[g] => %s_GNG_G2[g] + %d GDP[g]',peptide,9*NG,peptide,9*NG);
    reaction{3}.eq=  sprintf('%s_GNG_G2[g] + %d GDP-alpha-D-mannose[g] => %s_GNG_G3[g] + %d GDP[g]',peptide,30*NG,peptide,30*NG);
    reaction{4}.eq=  sprintf('%s_GNG_G3[g] + %d GDP-alpha-D-mannose[g] => %s_GNG_G4[g] + %d GDP[g]',peptide,NG,peptide,NG);

    for i=1:4
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s_GNG_G4',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end