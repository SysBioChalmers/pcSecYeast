function enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER)

enzymedata_all.proteinLoc = [enzymedata.proteinLoc;enzymedataSEC.proteinLoc;enzymedataMachine.proteinLoc];
enzymedata_all.enzyme = [enzymedata.enzyme;enzymedataSEC.enzyme;enzymedataMachine.enzyme];
enzymedata_all.enzyme_MW = [enzymedata.enzyme_MW;enzymedataSEC.enzyme_MW;enzymedataMachine.enzyme_MW];
enzymedata_all.proteins = [enzymedata.proteins;enzymedataSEC.proteins;enzymedataMachine.proteins];
enzymedata_all.proteinMWs = [enzymedata.proteinMWs;enzymedataSEC.proteinMWs;enzymedataMachine.proteinMWs];
enzymedata_all.proteinLength = [enzymedata.proteinLength;enzymedataSEC.proteinLength;enzymedataMachine.proteinLength];
enzymedata_all.subunit = [enzymedata.subunit;enzymedataSEC.subunit;enzymedataMachine.subunit];
enzymedata_all.subunit_stoichiometry = [enzymedata.subunit_stoichiometry;enzymedataSEC.subunit_stoichiometry;enzymedataMachine.subunit_stoichiometry];
enzymedata_all.kcat = [enzymedata.kcat;enzymedataSEC.kcat;enzymedataMachine.kcat];
enzymedata_all.label = [repmat({'Met'},length(enzymedata.enzyme),1);repmat({'Sec'},length(enzymedataSEC.enzyme),1);repmat({'Machine'},length(enzymedataMachine.enzyme),1)];
enzymedata_all.proteinPST = enzymedata.proteinPST;
enzymedata_all.proteinExtraMW = enzymedata.proteinExtraMW;
enzymedata_all.rxnscoef = [enzymedata.rxnscoef;enzymedataDummyER.rxnscoef];
enzymedata_all.rxns = [enzymedata.rxns;enzymedataDummyER.rxns];
enzymedata_all.comp = [repmat({''},length(enzymedata.enzyme),1);enzymedataSEC.comp];
enzymedata_all.kdeg = enzymedata.kdeg;
[enzymedata_all.proteins,a] = unique(enzymedata_all.proteins,'stable');
enzymedata_all.proteinLength = enzymedata_all.proteinLength(a);
enzymedata_all.proteinMWs = enzymedata_all.proteinMWs(a);
enzymedata_all.proteinLoc = enzymedata_all.proteinLoc(a);
enzymedata_all.proteinPST = enzymedata.proteinPST(a);
enzymedata_all.proteinExtraMW = enzymedata.proteinExtraMW(a);

for j = 1:length(enzymedata_all.enzyme)
    subunit = enzymedata_all.subunit(j,:);
    subunit_stoi = enzymedata_all.subunit_stoichiometry(j,:);
    [~,idx3] = ismember(subunit,enzymedata_all.proteins);
    idx3 = idx3(idx3~=0);
    enzymedata_all.enzymeSubunitMatrix(j,idx3) = subunit_stoi(subunit_stoi~=0);
end