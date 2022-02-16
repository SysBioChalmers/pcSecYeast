%% collectkcats
function enzymedata = collectkcats(model)

org_name = 'saccharomyces cerevisiae';
%% Load kcats data from BRENDA database
allkcats = loadkcats;

%% Collect EC number
[~, rawUniprot, ~] = xlsread('ec_number.xlsx','Uniprot');
[~, rawKEGG, ~] = xlsread('ec_number.xlsx','KEGG');

id_Uniprot = rawUniprot(:,1);
ec_Uniprot = rawUniprot(:,2);
id_KEGG = rawKEGG(:,1);
ec_KEGG = rawKEGG(:,2);

id_list = unique([id_Uniprot;id_KEGG]);

ecdata = struct();
ecdata.id = cell(0,1);
ecdata.ec = cell(0,1);

for i = 1:length(id_list)
    id = id_list(i);
    if ismember(id,id_Uniprot)
        ec_U_tmp = ec_Uniprot(ismember(id_Uniprot,id));
        ec_U_tmp = split(ec_U_tmp);
        id_U_tmp = repelem(id,length(ec_U_tmp))';
        ecdata.id = [ecdata.id;id_U_tmp];
        ecdata.ec = [ecdata.ec;ec_U_tmp];
    end
    if ismember(id,id_KEGG)
        ec_K_tmp = ec_KEGG(ismember(id_KEGG,id));
        id_K_tmp = repelem(id,length(ec_K_tmp))';
        ecdata.id = [ecdata.id;id_K_tmp];
        ecdata.ec = [ecdata.ec;ec_K_tmp];
    end
end

%% Assign kcat for each reaction with GPR

% collect enzyme data
enzyme_list = model.mets(contains(model.mets,'_complex') & startsWith(model.mets,'r_'));

idx_tmp = contains(model.rxns,'_complex_formation');
s_tmp = model.S(:,idx_tmp);
tf_tmp = s_tmp < 0;
max_subunit = max(sum(tf_tmp));% max subunit

enzymedata = struct();
enzyme_list = cellfun(@(x) x(1:strfind(x,'[')-1),enzyme_list,'UniformOutput',false);
enzymedata.enzyme = enzyme_list;
enzymedata.substrate = cell(length(enzyme_list),1);
enzymedata.subunit = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_stoichiometry = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_ec = cell(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat = zeros(length(enzyme_list),max_subunit);
enzymedata.subunit_kcat_conf = nan(length(enzyme_list),max_subunit);
enzymedata.kcat = zeros(length(enzyme_list),1);
enzymedata.kcat_conf = zeros(length(enzyme_list),1);

list = endsWith(model.rxns,'_translation');
list = model.rxns(list);
list = strrep(list,'r_','');
enzymedata.proteins = strrep(list,'_peptide_translation',''); % all protein
enzymedata.proteins = strrep(enzymedata.proteins,'_','-'); % all protein

for i = 1:length(enzyme_list)
    
    disp(['Assigning kcats:' num2str(i) '/' num2str(length(enzyme_list))]);
    
    enzyme_id = enzyme_list{i};
    
    % add substrates
    metrxn_id = split(enzyme_id,'_complex');
    metrxn_id = metrxn_id(1);
    idx_tmp = ismember(model.rxns,metrxn_id);
    s_tmp = model.S(:,idx_tmp);
    mets_tmp = model.metNames(s_tmp < 0);
    mets_tmp = cellfun(@(x) x(1:strfind(x,' [')-1),mets_tmp,'UniformOutput',false);
    mets_tmp = strjoin(mets_tmp,'; ');
    enzymedata.substrate(i,1) = {mets_tmp};
    
    % add subunits
    enzfmtrxn_id = strcat(metrxn_id,'_complex_formation');
    idx_tmp = ismember(model.rxns,enzfmtrxn_id);
    s_tmp = model.S(:,idx_tmp);
    subunits_tmp = model.mets(s_tmp < 0);
    na_tmp = repelem({''},max_subunit-length(subunits_tmp));
    subunits_tmp = cellfun(@(x) strrep(x,'_folding',''),subunits_tmp,'UniformOutput',false);
    subunits_tmp = cellfun(@(x) x(1:strfind(x,'[')-1),subunits_tmp,'UniformOutput',false);
    subunits_tmp = cellfun(@(x) strrep(x,'_','-'),subunits_tmp,'UniformOutput',false);
    enzymedata.subunit(i,:) = [subunits_tmp' na_tmp];
    
    % add stoichiometry of subunits
    stoichi_tmp = abs(full(s_tmp(s_tmp < 0)));
    na_tmp = repelem(0,max_subunit-length(stoichi_tmp));
    enzymedata.subunit_stoichiometry(i,:) = [stoichi_tmp' na_tmp];
    
    % add ec numbers of subunits
    for j = 1:length(subunits_tmp)
        subunit_tmp = subunits_tmp(j);
        if ismember(subunit_tmp,ecdata.id)
            ec_tmp = unique(ecdata.ec(ismember(ecdata.id,subunit_tmp)));
            ec_tmp = strjoin(ec_tmp,'; ');
        else
            ec_tmp = 'NA';
        end
        enzymedata.subunit_ec(i,j) = {ec_tmp};
    end
    
    % assign kcats and confidence scores for subunits
    for j = 1:length(subunits_tmp)
        ec_tmp = enzymedata.subunit_ec(i,j);
        substrate_tmp = enzymedata.substrate(i,1);
        [finalkcat_tmp, conf_score_tmp] = searchkcat(ec_tmp,substrate_tmp,org_name,allkcats);
        enzymedata.subunit_kcat(i,j) = finalkcat_tmp;
        enzymedata.subunit_kcat_conf(i,j) = conf_score_tmp;
    end
    
    % assign kcats and confidence scores for enzymes
    conf_tmp = max(enzymedata.subunit_kcat_conf(i,:));
    enzymedata.kcat_conf(i,1) = conf_tmp;
    finalkcat_list = enzymedata.subunit_kcat(i,:).*enzymedata.subunit_stoichiometry(i,:); %consider protein stoichiometry
    kcat_tmp = finalkcat_list(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    enzymedata.kcat(i,1) = median(kcat_tmp); %choose median among subunits
    %enzymedata.kcat(i,1) = min(kcat_tmp); %choose minimum among subunits
end

% no kcat assigned enzyme assumed to be the median of the collected dataset
medianvalue = median(enzymedata.kcat(~isnan(enzymedata.kcat)));
enzymedata.kcat(isnan(enzymedata.kcat)) = medianvalue;

% There are a few cases that multiple isozymes are assigned different kcats
% due to different EC numbers, e.g., EC4.2.1.- and EC4.2.1.36. This should
% be avoided by re-assigning the kcat with higher confidence score for the
% isozyme with lower confidence score?



    