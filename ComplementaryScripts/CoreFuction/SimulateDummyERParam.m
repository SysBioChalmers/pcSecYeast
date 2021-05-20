function [protein_info,enzymedataDummyER] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC)


protein_info = [protein_info;[{'dummyER'},{'dummyER'},1,1,meanprotein_info(3),meanprotein_info(4),meanprotein_info(5),0,meanprotein_info(7),'c',{''},1,{''},{''}]];
enzymedataDummyER.proteins = {'dummyER'};
enzymedataDummyER.proteinLength = 423;
enzymedataDummyER.proteinPST = meanprotein_info([3:5,7]);
enzymedataDummyER.proteinPSTInfo = {'DSB','NG','OG','GPI'};
enzymedataDummyER.proteinExtraMW_specific = enzymedataDummyER.proteinPST.*[0,9265,1080,2009];
enzymedataDummyER.proteinExtraMW = sum(enzymedataDummyER.proteinPST.*[0,9265,1080,2009]);
enzymedataDummyER.proteinMWs = 46000;

enzymedataDummyER.rxnscoef = [];
enzymedataDummyER.rxns = [];
SecComplex = enzymedataSEC.enzyme;
SecComplex_coef = enzymedataSEC.coefref;
for i = 1:length(SecComplex)
     rxns_by_this_complex = find(endsWith(model.rxns,['_',SecComplex{i}])& startsWith(model.rxns,'dummyER'));
     coefref_tmp = SecComplex_coef{i};
    %print kcat constriants in lp file
    if ~isempty(rxns_by_this_complex)
       rxnList = model.rxns(rxns_by_this_complex);

    DSB = enzymedataDummyER.proteinPST(1,strcmp(enzymedataDummyER.proteinPSTInfo,'DSB'));
    NG = enzymedataDummyER.proteinPST(1,strcmp(enzymedataDummyER.proteinPSTInfo,'NG'));
    OG = enzymedataDummyER.proteinPST(1,strcmp(enzymedataDummyER.proteinPSTInfo,'OG'));
    GPI = enzymedataDummyER.proteinPST(1,strcmp(enzymedataDummyER.proteinPSTInfo,'GPI'));
    if contains(coefref_tmp,'proteinLength')
        enzymedataDummyER.rxnscoef = [enzymedataDummyER.rxnscoef;enzymedataDummyER.proteinLength/467]; % bionumber 105224 average length
        enzymedataDummyER.rxns = [enzymedataDummyER.rxns;rxnList];
    elseif strcmp(coefref_tmp,'ProteinMW')
        enzymedataDummyER.rxnscoef = [enzymedataDummyER.rxnscoef;enzymedataDummyER.proteinMWs/54580]; % bionumber 115091 average weight
        enzymedataDummyER.rxns = [enzymedataDummyER.rxns;rxnList];
    else
        enzymedataDummyER.rxnscoef = [enzymedataDummyER.rxnscoef;repmat(eval(coefref_tmp),length(rxnList),1)]; 
        enzymedataDummyER.rxns = [enzymedataDummyER.rxns;rxnList];
    end
    end
end
end
