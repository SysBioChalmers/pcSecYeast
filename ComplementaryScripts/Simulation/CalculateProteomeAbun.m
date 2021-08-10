% this function is to calculate wthether the ER protein has a certain
% volume, by calculating the the abundance in different proteome data
%% ratio
%0.053 for ERmachine and 1.707683891812542e-01 for total,
%0.017 for ERmembrane
% 9.138417263261976e-03 of all HXTs ratio from the glucose phase QC3
%
[~,~,proteome] = xlsread('Proteome_collected.xlsx');

study = unique(proteome(1,:),'stable');

[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
load('enzymedataSEC.mat')
SecMachine = enzymedataSEC.proteins;
load('enzymedata.mat')
ERmembrane = enzymedata.proteins(strcmp(enzymedata.proteinLoc,'erm'));
translocator = {'YLR378C','YBR283C','YIL030C','YOL013C'}; % SEC61 SSH1 DOA10 HRD1
%ERmembrane = protein_info(strcmp(protein_info(:,10),'erm'),2); % er membrane proteins
%ERmembrane = intersect(ERmembrane,ERmembranemodel);
% calculateHXTabun
HXTs = {'YHR094C';'YFL011W';'YOL156W';'YIL170W';'YEL069C';'YNL318C';'YDL245C';'YJR158W';'YNR072W';'YMR011W';'YDR345C';'YHR092C';'YHR096C';'YDR343C';'YDR342C';'YJL214W';'YJL219W'};
% HXT10';'HXT11';'HXT12';'HXT13';'HXT14';'HXT15';'HXT16';'HXT17';'HXT2';'HXT3';'HXT4';'HXT5';'HXT6';'HXT7';'HXT8';'HXT9;
HXTAbun = [];
ERproteinAbun = [];
ERmembraneAbun = [];
SecMachineAbun = [];
translocatorAbun = [];
id = [];
id2 = [];
for i = 1:length(study)
    idx = ismember(proteome(1,:),study(i));
    proteome_tmp = proteome(3:end,idx);
    column_id = proteome(2:end,idx);
    condition_id = proteome(2,idx);
    condition_id = condition_id(2:end);
    gene_id = proteome_tmp(:,1);
    proteome_tmp = proteome_tmp(:,2:end);
    proteome_tmp(~cell2mat(cellfun(@ischar,gene_id,'UniformOutput',false)),:) = []; % delete empty line
    gene_id(~cell2mat(cellfun(@ischar,gene_id,'UniformOutput',false))) = []; % delete empty line
    abun = cell2mat(proteome_tmp);
    [~,idx] = ismember(ERprotein,gene_id);
    [~,idx2] = ismember(SecMachine,gene_id);
    [~,idx3] = ismember(ERmembrane,gene_id);
    [~,idx4] = ismember(HXTs,gene_id);
    [~,idx5] = ismember(translocator,gene_id);
    ERproteinAbun = [ERproteinAbun,sum(abun(idx(idx~=0),:))./sum(abun)];
    SecMachineAbun = [SecMachineAbun,sum(abun(idx2(idx2~=0),:))./sum(abun)];
    ERmembraneAbun = [ERmembraneAbun,sum(abun(idx3(idx3~=0),:))./sum(abun)];
    HXTAbun = [HXTAbun,sum(abun(idx4(idx4~=0),:))./sum(abun)];
    translocatorAbun = [translocatorAbun,sum(abun(idx5(idx5~=0),:))./sum(abun)]
    id = [id,repmat(study(i),1,length(proteome_tmp(1,:)))];
    id2 = [id2,condition_id];
end


[num_abd, raw_abd, ~] = xlsread('protein_abundance.xlsx','paxdb');
pax_abundance = num_abd(:,5);
pax_protein_list = raw_abd(2:end,2);
[~,idx] = ismember(ERprotein,pax_protein_list);
[~,idx2] = ismember(SecMachine,pax_protein_list);
[~,idx3] = ismember(ERmembrane,pax_protein_list);
 [~,idx4] = ismember(HXTs,pax_protein_list);
 [~,idx5] = ismember(translocator,pax_protein_list);
ERproteinAbun = [ERproteinAbun,sum(pax_abundance(idx(idx~=0),:))./sum(pax_abundance)];
id = [id,{'paxDb'}];
id2 = [id2,{'paxDb'}];
ERmembraneAbun = [ERmembraneAbun,sum(pax_abundance(idx3(idx3~=0),:))./sum(pax_abundance)];
SecMachineAbun = [SecMachineAbun,sum(pax_abundance(idx2(idx2~=0),:))./sum(pax_abundance)];
HXTAbun = [HXTAbun,sum(pax_abundance(idx4(idx4~=0),:))./sum(pax_abundance)];
    translocatorAbun = [translocatorAbun,sum(pax_abundance(idx5(idx5~=0),:))./sum(pax_abundance)]

boxplot(translocatorAbun,id,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);