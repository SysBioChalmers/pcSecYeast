%% addComplexFormationRxns 
function model = addComplexRxns(model,Protein_stoichiometry,protein_info,ProteinSequence)

% find reactions with GPR
idx = ~ismember(model.grRules,'');
metrxnid_list = model.rxns(idx);
gpr_list = model.grRules(idx);

for i = 1:length(metrxnid_list)
    disp(['Adding complex formation:' num2str(i) '/' num2str(length(metrxnid_list))]);
    
    metrxnname = metrxnid_list(i);
    cmplxid = strcat(cell2mat(metrxnname),'_complex[c]');
    rxnid = strcat(cell2mat(metrxnname),'_complex_formation');
    
    gpr = gpr_list(i);
    gpr_tmp = split(gpr,' and ');
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
    
    % degradation
%     metlist_deg = strrep(metlist,'folding','subunit');
%     coeflist_deg = [1*determined_stoich' 1*undetermined_stoich' -1];
%     rxnid_deg = strcat(cell2mat(metrxnname),'_complex_degradation');
%     model = addReaction(model,rxnid_deg,'metaboliteList',metlist_deg,'stoichCoeffList',coeflist_deg,'reversible',false);
%     
    % complex dilution
    metlist_dil = {cmplxid};
    coeflist_dil = -1;
    rxnid_dil = strcat(cell2mat(metrxnname),'_complex_dilution');
    model = addReaction(model,rxnid_dil,'metaboliteList',metlist_dil,'stoichCoeffList',coeflist_dil,'reversible',false);
end

% adding translation and degradation rxns for all peptides
geneList = unique(model.genes);
model = addfolding(model,geneList,protein_info); % add folding process for those peptide xx_peptide --> xx_folding
model = addTranslationRxns(model,ProteinSequence,geneList);
model = addDegradationRxns(model,ProteinSequence,protein_info,geneList);
end