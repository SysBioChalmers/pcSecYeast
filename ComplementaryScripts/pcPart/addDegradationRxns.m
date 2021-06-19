function  model = addDegradationRxns(model,ProteinSequence,protein_info,geneList)
%% Import AA IDs
[~,raw,~] = xlsread('aa_id.xlsx','cytoplasm');
aa_list = struct();
aa_list.aa = raw(2:end,1);
aa_list.subs = raw(2:end,3);
aa_list.prod = raw(2:end,5);
aa_list.prod_deg = raw(2:end,7);

[~,raw,~] = xlsread('aa_id.xlsx','energy');
e_list = struct();
e_list.subs = raw(2:end-1,3); % h2o atp gtp
e_list.prod = raw(2:end,5); % H+[c] ADP[c] GDP[c] phosphate[c]

%degradate Siginal peptide if it exists
for i = 1:length(geneList)
    disp(['Adding degradation rxn:' num2str(i) '/' num2str(length(geneList))]);
    geneid = geneList(i);
    protid = cell2mat(geneid);
    protid = strrep(protid,'-','_');
    protid = strcat(protid,'_peptide[c]');
    protid_tmp = strrep(protid,'_peptide[c]','');
    seq = ProteinSequence.seq(ismember(ProteinSequence.id,geneid));
    
    [~,geneidx] =ismember(geneid,protein_info(:,2));
    % can find this gene && through ER && SP == 1
    if geneidx ~= 0 && cell2mat(protein_info(geneidx,3)) == 1 && cell2mat(protein_info(geneidx,4)) == 1
        % SP degradation
        SP = cell2mat(protein_info(geneidx,14));
        seq = cell2mat(seq);
        SP_seq = seq(1:SP);
        [sum,energy] = countAA_deg(SP_seq,aa_list,e_list,false);
        protid_sp =  strrep(protid,'peptide','sp');
        rxnid = strcat('r_',protid_tmp,'_SP_degradation');
        metlist = [sum.prod' energy.subs_deg' energy.prod_deg' protid_sp];
        coeflist = [1*sum.num' -1*energy.subs_degnum' energy.prod_degnum' -1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);

        % subunit degradation taking off the SP
        seq_withoutSP = seq(SP+1:end);
        [sum,energy] = countAA_deg(seq_withoutSP,aa_list,e_list,true);
        protid =  strrep(protid,'peptide','subunit');
        rxnid = strcat('r_',protid_tmp,'_subunit_degradation');
        metlist = [sum.prod' energy.subs_deg' energy.prod_deg' protid];
        coeflist = [1*sum.num' -1*energy.subs_degnum' energy.prod_degnum' -1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    else
        % degradation
        [sum,energy] = countAA_deg(seq,aa_list,e_list,true);
        protid =  strrep(protid,'peptide','subunit');
        rxnid = strcat('r_',protid_tmp,'_subunit_degradation');
        metlist = [sum.prod' energy.subs_deg' energy.prod_deg' protid];
        coeflist = [1*sum.num' -1*energy.subs_degnum' energy.prod_degnum' -1];
        model = addReaction(model,rxnid,'metaboliteList',metlist,'stoichCoeffList',coeflist,'reversible',false);
    end
end
