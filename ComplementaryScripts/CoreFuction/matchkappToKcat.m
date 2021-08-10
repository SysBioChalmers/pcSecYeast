function enzymedata = matchkappToKcat(enzymedata,kapp4)

% This function is to map kapp into the enzymedata.kcat
kapp4.max(strcmp(kapp4.rxn,'r_1166')) = [];
kapp4.rxn(strcmp(kapp4.rxn,'r_1166')) = [];
% match reaction & rev
rxnList_kapp4 = extractBefore(kapp4.rxn,7);
rvs_kapp4 = extractAfter(kapp4.rxn,6);
rxnList_enzymedata = extractBefore(enzymedata.enzyme,7);
rvs_enzymedata = extractAfter(enzymedata.enzyme,6);
for i = 1:length(rxnList_kapp4)
    
    if ~isempty(rvs_kapp4{i})
        rxnTmp = find(contains(rxnList_enzymedata,rxnList_kapp4{i}) &contains(rvs_enzymedata,rvs_kapp4{i}));
    else
        rxnTmp = find(contains(rxnList_enzymedata,rxnList_kapp4{i}));
    end
    stoi = enzymedata.subunit_stoichiometry(rxnTmp,1);
    enzymedata.kcat(rxnTmp) = stoi*kapp4.max(i)*3600;
end

end