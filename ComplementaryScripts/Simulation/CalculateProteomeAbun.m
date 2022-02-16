% this function is to calculate wthether the ER protein has a certain
% volume, by calculating the the abundance in different proteome data
%% ratio percent of the total proteome
%0.053 for SECmachine and 0.1708 for ERprotein,
%0.0174 for ERmembrane
% 9.138417263261976e-03 of all HXTs ratio from the glucose phase QC3
% 0.00271 for ERAD_proteinAbun
% 0.00025 8.22e-5 for translocator from SCE_petri

[~,~,proteome] = xlsread('Proteome_collected.xlsx');

study = unique(proteome(1,:),'stable');

[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
load('enzymedataSEC.mat')
SecMachine = enzymedataSEC.proteins;
load('enzymedata.mat')
ERmembrane = enzymedata.proteins(strcmp(enzymedata.proteinLoc,'erm'));
translocator = {'YIL030C','YOL013C'}; % SEC61 SSH1 DOA10 HRD1
transloc_complex = {'YBR201W';'YDR057W';'YER100W';'YIL030C';'YLR207W';'YML029W';'YMR022W';'YMR264W';'YOL013C'};
%ERmembrane = protein_info(strcmp(protein_info(:,10),'erm'),2); % er membrane proteins
%ERmembrane = intersect(ERmembrane,ERmembranemodel);
% calculateHXTabun
HXTs = {'YHR094C';'YFL011W';'YOL156W';'YIL170W';'YEL069C';'YNL318C';'YDL245C';'YJR158W';'YNR072W';'YMR011W';'YDR345C';'YHR092C';'YHR096C';'YDR343C';'YDR342C';'YJL214W';'YJL219W'};
% HXT10';'HXT11';'HXT12';'HXT13';'HXT14';'HXT15';'HXT16';'HXT17';'HXT2';'HXT3';'HXT4';'HXT5';'HXT6';'HXT7';'HXT8';'HXT9;
ERAD_protein = {'YJL034W','YCL043C','YHR204W','YCL043C','YMR264W','YER100W','YMR022W','YDR057W','YLR207W','YOL013C','YBR201W','YML029W','YMR264W','YER100W','YIL030C','YMR022W','YER087C-B','YDR086C','YBR283C','YML013W','YDL126C','YGR048W','YBR170C','YMR276W','YEL037C','YPL096W','YKL210W','YCL043C','YML130C','YJL034W','YMR264W','YER100W','YMR022W','YLR207W','YOL013C','YBR201W','YMR276W','YEL037C','YKL210W'};
% according to the proteins annotaion in [~,~,proteins]=xlsread('TableS1.xlsx','Secretory');
HXTAbun = [];
ERproteinAbun = [];
ERmembraneAbun = [];
SecMachineAbun = [];
translocatorAbun = [];
ERAD_proteinAbun = [];
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
    [~,idx6] = ismember(ERAD_protein,gene_id);
    ERproteinAbun = [ERproteinAbun,sum(abun(idx(idx~=0),:),1)./sum(abun,1)];
    SecMachineAbun = [SecMachineAbun,sum(abun(idx2(idx2~=0),:),1)./sum(abun,1)];
    ERmembraneAbun = [ERmembraneAbun,sum(abun(idx3(idx3~=0),:),1)./sum(abun,1)];
    HXTAbun = [HXTAbun,sum(abun(idx4(idx4~=0),:),1)./sum(abun,1)];
    translocatorAbun = [translocatorAbun,sum(abun(idx5(idx5~=0),:),1)./sum(abun,1)];
    ERAD_proteinAbun = [ERAD_proteinAbun,sum(abun(idx6(idx6~=0),:),1)./sum(abun,1)];
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
     [~,idx6] = ismember(ERAD_protein,pax_protein_list);

ERproteinAbun = [ERproteinAbun,sum(pax_abundance(idx(idx~=0),:),1)./sum(pax_abundance,1)];
id = [id,{'paxDb'}];
id2 = [id2,{'paxDb'}];
ERmembraneAbun = [ERmembraneAbun,sum(pax_abundance(idx3(idx3~=0),:),1)./sum(pax_abundance,1)];
SecMachineAbun = [SecMachineAbun,sum(pax_abundance(idx2(idx2~=0),:),1)./sum(pax_abundance,1)];
HXTAbun = [HXTAbun,sum(pax_abundance(idx4(idx4~=0),:),1)./sum(pax_abundance,1)];
translocatorAbun = [translocatorAbun,sum(pax_abundance(idx5(idx5~=0),:),1)./sum(pax_abundance,1)];
ERAD_proteinAbun = [ERAD_proteinAbun,sum(pax_abundance(idx6(idx6~=0),:),1)./sum(pax_abundance,1)];

boxplot(ERAD_proteinAbun,id,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[56,108,176]/255);