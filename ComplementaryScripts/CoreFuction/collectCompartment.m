%% collectCompartment
% collect a list of proteins in a compartment/compartment(s)
function [genes, reactions] = collectCompartment(model,comps_id)

if ischar(comps_id)
    comps_id = {comps_id};
end
comps_id = cellfun(@(x) strcat('[',x,']'),comps_id,'UniformOutput',false);

metidx = contains(model.mets,comps_id);
rxnidx = any(full(model.S(metidx,:)));
reactions = model.rxns(rxnidx' & ~ismember(model.rules,''));
gpr_list = model.rules(rxnidx');
gpr_list = unique(gpr_list);

gene_list = cell(0,1);
for i = 1:length(gpr_list)
    gpr = gpr_list{i};
    if ~ismember(gpr,'')
        gpr = strrep(gpr,'(( ','');
        gpr = strrep(gpr,' ))','');
        gpr = strrep(gpr,'( ','');
        gpr = strrep(gpr,' )','');
        gpr = strtrim(gpr);
        gpr = strrep(gpr,' | ',' & ');
        gpr = strsplit(gpr,' & ');
        gene_list = [gene_list;gpr'];
    end
end
gene_list = unique(gene_list);
gene_list = strrep(gene_list,'x(','');
gene_list = strrep(gene_list,')','');
gene_idx = cell2mat(cellfun(@(x) str2double(x),gene_list,'UniformOutput',false));
genes = model.genes(gene_idx);
genes = unique(genes);

