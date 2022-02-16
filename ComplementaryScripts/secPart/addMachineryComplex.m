function model=addMachineryComplex(model,proteins,protein_info,Protein_stoichiometry,ProteinSequence)
%load secretory machnism genes and add trnalstion and degradation for those
%unique genes
%load secretory complexes information and extract subunit information from
%subunit sheet
[~,idx] = unique(proteins(:,1),'stable');
complex_list = proteins(idx,1);
compartment = proteins(idx,24);

for i = 2:numel(complex_list)%1 is the title line
    disp(['Adding complex:' num2str(i) '/' num2str(length(complex_list))]);
    
    % add complex formation
    metrxnname = complex_list(i);
    cmplxid = strcat(cell2mat(metrxnname),'[',compartment{i},']');
    rxnid = strcat(cell2mat(metrxnname),'_formation');
    
    gpr_tmp = proteins((find(strcmp(proteins(:,1),complex_list(i)))),2);
    gpr_tmp_tmp = cellfun(@(x) strrep(x,'-','_'),gpr_tmp,'UniformOutput',false);
    
    [~,geneidx] =ismember(gpr_tmp,protein_info(:,2));
    peptide_comp = protein_info(geneidx,10); %peptide compartment
    subs_id = cellfun(@(x) strcat(x,'_folding'),gpr_tmp_tmp,'UniformOutput',false);
    subs_id = strcat(subs_id,strcat(repmat({'['},length(subs_id),1),peptide_comp,repmat({']'},length(subs_id),1)));
    
    % check if protein stoichiometry is determined
    idx_tmp = ismember(gpr_tmp,Protein_stoichiometry.protein);
    determined = gpr_tmp(idx_tmp);
    determined_subs_id = subs_id(idx_tmp);
    [~,b] = ismember(determined,Protein_stoichiometry.protein);
    determined_stoich = Protein_stoichiometry.stoichiometry(b);
    undetermined_subs_id = subs_id(~idx_tmp);
    undetermined_stoich = ones(length(undetermined_subs_id),1); % assumed to be 1 if not determined
    
    metlist = [determined_subs_id' undetermined_subs_id' cmplxid];
    coeflist = [-1*determined_stoich' -1*undetermined_stoich' 1];
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    
%   % complex degradation complex is not degraded for now in the model (only protein subunit is degraded), but which
%     can be added back in the model later if needed
%     metlist_deg = strrep(metlist,'folding','subunit');
%     coeflist_deg = [1*determined_stoich' 1*undetermined_stoich' -1];
%     rxnid_deg = strcat(cell2mat(metrxnname),'_complex_degradation');
%     model = addReaction(model,rxnid_deg,'metaboliteList',metlist_deg,'stoichCoeffList',coeflist_deg,'reversible',false);
%     
    % complex dilution
    metlist_dil = {cmplxid};
    coeflist_dil = -1;
    rxnid_dil = strcat(cell2mat(metrxnname),'_dilution');
    model = addReaction(model,rxnid_dil,'metaboliteList',metlist_dil,'stoichCoeffList',coeflist_dil,'reversible',false);
end

% adding translation and degradation rxns for all peptides
geneList = unique(proteins(2:end,2)); % the first one is title
model = addfolding(model,geneList,protein_info);% add folding process for those peptide xx_peptide --> xx_folding
model = addTranslationRxns(model,ProteinSequence,geneList);
model = addDegradationRxns(model,ProteinSequence,protein_info,geneList);

end