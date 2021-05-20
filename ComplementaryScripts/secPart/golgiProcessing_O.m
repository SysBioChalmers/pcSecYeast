function [newModel,peptide_name,rxns] = golgiProcessing_O(model,peptide,Length_total,NG,onlyrxns)
rxns =[];
if NG >0
    Length = Length_total/40;
    
    reaction{1}.rxns = sprintf('%s_GLOG_Golgi_O_linked_manosylation_I_sec_Kre2p_ktr1p_ktr3p_complex',peptide);
    reaction{2}.rxns = sprintf('%s_GLOG_Golgi_O_linked_manosylation_II_sec_Mnn1p_complex',peptide);
    
    reaction{1}.rxnNames = sprintf('%s_GLOG_Golgi_O_linked_manosylation_I_sec_Kre2p_ktr1p_ktr3p_complex',peptide);
    reaction{2}.rxnNames = sprintf('%s_GLOG_Golgi_O_linked_manosylation_II_sec_Mnn1p_complex',peptide);
    
    
    reaction{1}.eq = sprintf('%s[g] + %d GDP-alpha-D-mannose[g] => %s_GOG_G1[g] + %d GDP[g]',peptide,3*NG,peptide,3*NG);
    reaction{2} .eq=  sprintf('%s_GOG_G1[g] + %d GDP-alpha-D-mannose[g] => %s_GOG_G2[g] + %d GDP[g]',peptide,2*NG,peptide,2*NG);

    for i=1:2
        if onlyrxns == 1
            rxns = [rxns;{reaction{i}.rxns}];
        else
            model=addYeastReaction(model,reaction{i}.eq,{reaction{i}.rxns},{reaction{i}.rxnNames});
            rxns = [rxns;{reaction{i}.rxns}];
        end
    end
    newModel = model;
    peptide_name = sprintf('%s_GOG_G2',peptide);
else
    newModel = model;
    peptide_name = peptide;
    
end