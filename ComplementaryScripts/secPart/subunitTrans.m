function   [model,rxns] = subunitTrans(model,peptide,comps)

%this function is to add a subunit trnsport reaction for subunit which is
%for degradation of complex. Note: only the peptide that doesn't exist in
%cytosol will be added

%Feiran Li, 2019-04-12
rxns = [];
if strcmp(comps,'c') ~=1
    rxnID=sprintf('%s_subunit_export_%s',peptide,char(comps));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('export %s subunit from %s',peptide,char(comps));
    r=sprintf('%s_subunit[%s] => %s_subunit[c]',peptide,char(comps),peptide);
    model=addYeastReaction(model,r,rxnID,{rxnName});

end
    
end
                 


