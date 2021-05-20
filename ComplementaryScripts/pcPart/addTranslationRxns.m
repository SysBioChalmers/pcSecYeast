%% addTranslationRxns 
function model = addTranslationRxns(model,ProteinSequence,geneList)
% There are several assumptions:
% 1. all proteins are translated in cytoplasm;
% 2. 2ATP 2GTP were for elongation for each aa as energy cost.
% 3. no machineries included, e.g., ribosomes.

%% Import AA IDs
[~,raw,~] = xlsread('aa_id.xlsx','cytoplasm');
aa_list = struct();
aa_list.aa = raw(2:end,1);
aa_list.subs = raw(2:end,3);
aa_list.prod = raw(2:end,5);

[~,raw,~] = xlsread('aa_id.xlsx','energy');
e_list = struct();
e_list.subs = raw(2:end-1,3); % h2o atp gtp
e_list.prod = raw(2:end,5); % H+[c] ADP[c] GDP[c] phosphate[c]

%% Add reactions
for i = 1:length(geneList)
    disp(['Adding translation:' num2str(i) '/' num2str(length(geneList))]);
    geneid = geneList(i);
    seq = ProteinSequence.seq(ismember(ProteinSequence.id,geneid));
    [sum, energy] = countAA(seq,aa_list,e_list);
    protid = cell2mat(geneid);
    protid = strrep(protid,'-','_');
    protid = strcat(protid,'_peptide[c]');
    protid_tmp = strrep(protid,'[c]','');
    
    % translation
    rxnid = strcat('r_',protid_tmp,'_translation');
    metlist = [sum.subs' energy.subs' energy.prod' sum.prod' protid];
    coeflist = [-1*sum.num' -1*energy.subnum' energy.prodnum' sum.num' 1];
    model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    
end