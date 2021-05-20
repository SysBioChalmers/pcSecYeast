%% updateComplexformation
function model = updateComplexformation(model,protein_info)

[~,~,raw] = xlsread('manual_update.xlsx','stoichiometry');

changed_enzymes = raw(2:end,1);
subunits = raw(2:end,2);
coeffs = raw(2:end,3);

for i = 1:length(changed_enzymes)
    enzymename = changed_enzymes{i};
    enzymeformationrxn = strcat(enzymename,'_formation');
    rxnidx = ismember(model.rxns,enzymeformationrxn);
    
    changed_subunits = subunits{i};
    changed_subunits = split(changed_subunits,'; ');
    [~,geneidx] =ismember(changed_subunits,protein_info(:,2));
    peptide_comp = protein_info(geneidx,10); %peptide compartment
    changed_subunits = cellfun(@(x) strcat(x,'_folding'),changed_subunits,'UniformOutput',false);
    changed_subunits = cellfun(@(x) strrep(x,'-','_'),changed_subunits,'UniformOutput',false);
    
    changed_subunits = strcat(changed_subunits,strcat(repmat({'['},length(changed_subunits),1),peptide_comp,repmat({']'},length(changed_subunits),1)));
    
    changed_coeffs = coeffs{i};
    
    if length(changed_subunits) > 1
        changed_coeffs = split(changed_coeffs,'; ');
        for j = 1:length(changed_subunits)
            subunit_tmp = changed_subunits(j);
            coeff_tmp = -1 * str2double(changed_coeffs{j});
            metidx_tmp = ismember(model.mets,subunit_tmp);
            model.S(metidx_tmp,rxnidx) = coeff_tmp;
        end
    else
        if ischar(changed_coeffs)
            changed_coeffs = str2double(changed_coeffs);
        end
        coeff_tmp = -1 * changed_coeffs;
        metidx_tmp = ismember(model.mets,changed_subunits);
        model.S(metidx_tmp,rxnidx) = coeff_tmp;
    end
    
    % degradation
    enzymedegradationrxn = strcat(enzymename,'_degradation');
    rxnidx = ismember(model.rxns,enzymedegradationrxn);
    
    changed_subunits = subunits{i};
    changed_subunits = split(changed_subunits,'; ');
    [~,geneidx] =ismember(changed_subunits,protein_info(:,2));
    peptide_comp = protein_info(geneidx,10); %peptide compartment
    changed_subunits = cellfun(@(x) strcat(x,'_subunit'),changed_subunits,'UniformOutput',false);
    changed_subunits = cellfun(@(x) strrep(x,'-','_'),changed_subunits,'UniformOutput',false);
    
    changed_subunits = strcat(changed_subunits,strcat(repmat({'['},length(changed_subunits),1),peptide_comp,repmat({']'},length(changed_subunits),1)));
    
    changed_coeffs = coeffs{i};
    
    if length(changed_subunits) > 1
        changed_coeffs = split(changed_coeffs,'; ');
        for j = 1:length(changed_subunits)
            subunit_tmp = changed_subunits(j);
            coeff_tmp = str2double(changed_coeffs{j});
            metidx_tmp = ismember(model.mets,subunit_tmp);
            model.S(metidx_tmp,rxnidx) = coeff_tmp;
        end
    else
        coeff_tmp = changed_coeffs;
        metidx_tmp = ismember(model.mets,changed_subunits);
        model.S(metidx_tmp,rxnidx) = coeff_tmp;
    end
end
