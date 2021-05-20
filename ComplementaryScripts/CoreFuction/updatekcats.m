%% updateComplexformation 
function enzymedata = updatekcats(enzymedata)

[num,txt,~] = xlsread('manual_update.xlsx','kcats');
enzyme_list = txt(2:end,1);
kcat_list = num(:,1) * 3600;
conf_list = num(:,3);
for i = 1:length(enzyme_list)
    enzymename = enzyme_list(i);
    kcat_tmp = kcat_list(i);
    conf_tmp = conf_list(i);
    idx = ismember(enzymedata.enzyme,enzymename);
    enzymedata.kcat(idx) = kcat_tmp;
    enzymedata.kcat_conf(idx) = conf_tmp;
end
