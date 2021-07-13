function reaction = transportFromGolgiToOther(peptide,peptide_org,compartment)
comps = compartment;
%change compartments(only in old version of raven)
% CONValldata = cat(2,model.comps,model.compNames);
% [m, n]      = size(CONValldata);
% for i = 1:m
%     aa = CONValldata(i,1);
%     aa = char(aa);
%     bb = comps;
%     bb = char(bb);
%     if strcmp(bb,aa)
%         comps = CONValldata(i,2);
%     end
% end
%met = [peptide_org,'_folding [',char(comps),']'];
met = [peptide_org,'_folding[',char(comps),']'];
reaction{1}.rxns = sprintf('%s_transportFromGolgiToOthercompartment',peptide_org);

reaction{1}.rxnNames = sprintf('%s_transportFromGolgiToOthercompartment',peptide);

reaction{1}.eq = sprintf('%s[g] => %s',peptide,met);

    
end