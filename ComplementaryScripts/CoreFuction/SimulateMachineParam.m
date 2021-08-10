function enzymedata = SimulateMachineParam(model)
% this function calculate the MW for machinery complex and also assgin kcat
% for those enzymes

[~,~,chap]=xlsread('TableS1.xlsx','Machinery');
MachineComplex = unique(chap(:,1),'stable');
[~,Idx] = ismember(MachineComplex,chap(:,1));
MachineComplex_comp = chap(Idx,24);
MachineComplex = MachineComplex(2:end,:);
MachineComplex_comp = MachineComplex_comp(2:end,:);

enzymedata = struct();

enzymedata.enzyme = MachineComplex;
enzymedata.comp = MachineComplex_comp;
idx_tmp = contains(model.rxns,'_complex_formation');
s_tmp = model.S(:,idx_tmp);
tf_tmp = s_tmp < 0;
max_subunit = max(sum(tf_tmp));% max subunit

enzymedata.subunit = cell(length(enzymedata.enzyme),max_subunit);

% collect enzyme data
for i = 1:length(enzymedata.enzyme)
    
    disp(['collect enzyme info' num2str(i) '/' num2str(length(enzymedata.enzyme))]);
    
    enzyme_id = enzymedata.enzyme{i};
    
    % add subunits
    enzfmtrxn_id = strcat(enzyme_id,'_formation');
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
    
end

% define parameter
enzymedata.kcat(1,1) = 30*60*60; % /h 107785 bionumber % 372763 ribosome rate withgrowth rate aa/h
enzymedata.kcat(2,1) = 2000*60; % /h  106538 bionumber ribosome/h
enzymedata.kcat(3,1) = 1*60; % /h

enzymedata.proteins = strrep(setdiff(unique(enzymedata.subunit(:)),''),'_','-');

end

