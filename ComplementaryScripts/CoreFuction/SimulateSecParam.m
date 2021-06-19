function [enzymedataSEC,modeled_ratio,meanprotein_info,missingsecP_ratio] = SimulateSecParam(model,protein_info,ProteinSequence)
% This function is to get all parameters for secretory pathway
% Input of this function is all protein abundance from paxdb database
% Use this to get the simulation for all parameters in the secretory
% pathway
% missingP_ratio is in percent of mg
% Feiran Li -20190507

% Load the proteome data
[num_abd, raw_abd, ~] = xlsread('protein_abundance.xlsx','paxdb');
[num_mw, raw_mw, ~] = xlsread('protein_abundance.xlsx','uniprot');

% Calculate total mass
pax_protein_list = raw_abd(2:end,2);
pax_abundance = num_abd(:,3);
pax_abundance_mg = num_abd(:,5);
% Load secretory complexes
[~,~,chap]=xlsread('TableS1.xlsx','Secretory');
SecComplex = unique(chap(2:end,1),'stable');
Secgene = chap(2:end,1:2);
[~,Idx] = ismember(SecComplex,chap(:,1));
SecComplex_comp = chap(Idx,24);
SecComplex_coef = chap(Idx,25);
SecComplex_func = chap(Idx,26);

%% 1) get all rxnIDs for all proteins in proteome
rxnList = [];
geneError = [];
for i = 1:length(pax_protein_list)
    [model,allrxns,geneError_tmp] = addModificationforSec(model,pax_protein_list(i),protein_info,1);
    rxnList = [rxnList;allrxns];
    geneError = [geneError;geneError_tmp];
end

% match proteome abundance with rxnList
matchingList = cell(length(rxnList),3);

for i = 1:length(rxnList)
    S = regexp(rxnList{i}, '_', 'split');
    if length(S{2}) == 1
        gene = strcat(S{1},'_',S{2});
    else
        gene = {S{1}};
    end
    gene = strrep(gene,'-','_');
    [~,Idx] = ismember(gene,strrep(pax_protein_list,'-','_'));
    if Idx ~=0
        matchingList(i,:) = [gene,rxnList(i),pax_abundance(Idx)];
    else
        matchingList(i,:) = [gene,rxnList(i),0];
    end
end

% match the coef for the rxnList
enzymedata.proteins = setdiff(pax_protein_list,geneError);
enzymedata = calculateMW(enzymedata,ProteinSequence,protein_info);
enzymedataSEC.enzyme = SecComplex;
enzymedataSEC.comp = SecComplex_comp;
enzymedataSEC.coefref = SecComplex_coef;
a.rxns = rxnList;
enzymedata = SimulateRxnKcatCoef(a,enzymedataSEC,enzymedata);
[~,idx] = ismember(enzymedata.rxns,rxnList);
matchingList(idx(idx~=0),4) = num2cell(enzymedata.rxnscoef(idx~=0)); % coef
clear enzymedata


idx_tmp = contains(model.rxns,'_complex_formation');
s_tmp = model.S(:,idx_tmp);
tf_tmp = s_tmp < 0;
max_subunit = max(sum(tf_tmp));

for i = 1:length(enzymedataSEC.enzyme)
    
    disp(['collect enzyme info' num2str(i) '/' num2str(length(enzymedataSEC.enzyme))]);
    
    enzyme_id = enzymedataSEC.enzyme{i};
    
    % add subunits
    enzfmtrxn_id = strcat(enzyme_id,'_formation');
    idx_tmp = ismember(model.rxns,enzfmtrxn_id);
    s_tmp = model.S(:,idx_tmp);
    subunits_tmp = model.mets(s_tmp < 0);
    na_tmp = repelem({''},max_subunit-length(subunits_tmp));
    subunits_tmp = cellfun(@(x) strrep(x,'_folding',''),subunits_tmp,'UniformOutput',false);
    subunits_tmp = cellfun(@(x) x(1:strfind(x,'[')-1),subunits_tmp,'UniformOutput',false);
    subunits_tmp = cellfun(@(x) strrep(x,'_','-'),subunits_tmp,'UniformOutput',false);
    enzymedataSEC.subunit(i,:) = [subunits_tmp' na_tmp];
    
     % add stoichiometry of subunits
    stoichi_tmp = abs(full(s_tmp(s_tmp < 0)));
    na_tmp = repelem(0,max_subunit-length(stoichi_tmp));
    enzymedataSEC.subunit_stoichiometry(i,:) = [stoichi_tmp' na_tmp];
    
    % get subunit abun
    [~,genes_index_proteome] = ismember(subunits_tmp,pax_protein_list);
    subunit_count = pax_abundance(genes_index_proteome(genes_index_proteome ~= 0))';
    E0(i,1:length(subunits_tmp)) = subunit_count;
    enzymedataSEC.subunit_abun(i,1:length(subunits_tmp)) = subunit_count;
    % get catalyzed all enzyme abun * coef
    All_E = find(endsWith(matchingList(:,2),['_',enzyme_id]));
    E_sum(i,1:length(subunits_tmp)) = sum(cell2mat(matchingList(All_E,3)).*cell2mat(matchingList(All_E,4)));

end

%  ERAD should be only 30%  of total protein
u = 0.4;
E_sum(strcmp(SecComplex_func,'ERAD'),:) = E_sum(strcmp(SecComplex_func,'ERAD'),:)*0.045/(u+0.045);

% sum(V) <= Vsyn = kcat[E]
% (mu + kdeq)*sum([E]) <= kcat[E0]
allsecgene = setdiff(unique(enzymedataSEC.subunit(:)),'');
for i = 1:length(allsecgene)
    idx = ismember(enzymedataSEC.subunit,allsecgene(i));
    if length(find(idx)) > 1
        i
    E_sum(idx) = sum(E_sum(idx));
    end
end


kcat_tmp = (E_sum./E0).* enzymedataSEC.subunit_stoichiometry(:,1:length(E_sum(1,:)));

%enzymedataSEC.kcat = median(kcat_tmp,2,'omitnan')*(u+0.045);
enzymedataSEC.kcat = min(kcat_tmp,[],2)*(u+0.045);
enzymedataSEC.proteins = strrep(setdiff(unique(enzymedataSEC.subunit(:)),''),'_','-'); % get all proteins involved in the sec pathway

% calculate modeled_protein coverage
list = endsWith(model.rxns,'_translation');
list = model.rxns(list);
list = strrep(list,'r_','');
modeled_proteins = strrep(list,'_peptide_translation','');

% calculate modeled protein ratio gram
[~,Idx] = ismember(modeled_proteins,strrep(pax_protein_list,'-','_'));
modeled_ratio = sum(pax_abundance_mg(Idx(Idx~=0)))/sum(pax_abundance_mg);

% proteins not take into account
pax_protein_list = strrep(pax_protein_list,'-','_');
missingprotein = setdiff(pax_protein_list,modeled_proteins);
[~,Idx] = ismember(missingprotein,pax_protein_list);
[~,Idx2] = ismember(missingprotein,strrep(protein_info(:,2),'-','_'));
missingprotein = missingprotein(Idx2~=0);
[~,Idx] = ismember(missingprotein,pax_protein_list);
[~,Idx2] = ismember(missingprotein,strrep(protein_info(:,2),'-','_'));

mean = sum(pax_abundance(Idx).*cell2mat(protein_info(Idx2,3:9)))/sum(pax_abundance(Idx));
mean_length = sum(pax_abundance(Idx).*cell2mat(protein_info(Idx2,12)))/sum(pax_abundance(Idx));
meanprotein_info = mean/mean_length*423;% 4.23 is the unmodeled protein secretory length*100 to save unnesaacy GTP and ATP

%% unmodeled protein secretory ratio
rxnList = [];
for i = 1:length(missingprotein)
    [model,allrxns] = addModificationforSec(model,missingprotein{i},protein_info,1);
    rxnList = [rxnList;allrxns];
end

% find unmodeled secretion protein list
pax_protein_list = strrep(pax_protein_list,'-','_');
for i = 1:length(rxnList)
    S = regexp(rxnList{i}, '_', 'split');
    if length(S{2}) == 1
        gene = strcat(S{1},'_',S{2});
    else
        gene = {S{1}};
    end
    gene = strrep(gene,'-','_');
    [~,Idx] = ismember(gene,pax_protein_list);
    missingsecP(i) = gene;
end

% get missingSecP
missingsecP = unique(missingsecP);
[~,Idx] = ismember(missingsecP,pax_protein_list);
missingsecP_ratio = sum(pax_abundance_mg(Idx(Idx~=0)))/sum(pax_abundance_mg);


end
