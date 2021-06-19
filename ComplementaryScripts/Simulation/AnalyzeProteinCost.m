%% Plot for Figure 2C 2D 2E 2F

% analyze ProteinCost
[~,~,protein_info] = xlsread('Protein_Information.xlsx');
protein_info = protein_info(2:end,:);
%ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
ERprotein = protein_info(strcmp(protein_info(:,10),'e')|strcmp(protein_info(:,10),'ce'),2); % go through ER
ERprotein = strrep(ERprotein,'-','');
ERprotein = strcat(ERprotein,'new');
TP = {'Amylase';'Insulin';'HBsAg';'Adalimumabcomplex';'Etanerceptcomplex';'Filgrastim';'Hemoglobincomplex';'Nanobody';'Pex';'Rancomplex';'TSHaCGA';'TSHb'};
ERprotein = [TP;ERprotein];

for i = 1:20
load([ERprotein{20*(i-1)+1},'_',num2str(20*(i-1)+1),'.mat'])
res_slope_final(20*(i-1)+1:20*(i),1,:) = res_slope_ribo(20*(i-1)+1:20*(i),1,:);
res_glc_final(20*(i-1)+1:20*(i),1,:) = res_slope_glc(20*(i-1)+1:20*(i),1,:);
end
figure
hold on
h = histogram(abs(1./res_glc_final(:,1)));
h.FaceAlpha = 0.3;
h.FaceColor = [202,0,32]/255;
h.EdgeColor = [202,0,32]/255;
[~,idx] = ismember(TP,ERprotein);
barnum = discretize(abs(1./res_glc_final(idx,1)),h.BinEdges);
x = abs(1./res_glc_final(idx,1)); %TP glc cost
y = h.Values(barnum)+10;
h2 = stem(x,y,'Color','k','MarkerSize',2);
[~,y1] = sort(abs(1./res_glc_final(idx,1)),'ascend');
TPlabelorder = TP(y1);
box on
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Frequency','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Glucose cost','FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
ylim([0,120])
for i = 1:length(TP)
    text(x(i),y(i)+5,TP{i},'FontSize',6,'FontName','Helvetica','Color','k')
end

[num_abd,raw_abd,~] = xlsread('Proteome_collected.xlsx','mRNA_Qi2');
for i = 1:length(num_abd(1,:))
pax_abundance = num_abd(:,i);
pax_protein_list = raw_abd(2:end,1);
[~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
m = [];
m(idx~=0,1) = pax_abundance(idx(idx~=0));
m = m(1:400);
idx = find(m ~= 0 & ~isnan(m)&res_glc_final(:,1) ~= 0);
m = m(idx);
n = res_glc_final(idx,1);
[RHO4,PVAL4] = corr(log10(abs(1./n)),log10(m),'Type','Pearson');
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = [];
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400);
% idx = find(m ~= 0&~isnan(m)&res_slope_final(:,1) ~= 0)
% m = m(idx);
% n = res_slope_final(idx,1);
% [RHO4,PVAL4] = corr(log10(abs(1./n)),log10(m),'Type','Pearson');
x(i,1) = RHO4;
x(i,2) = PVAL4;
x(i,3) = length(m);
end

group = [1,1,2,2,3,3,1,1,2,2,3,3];
figure
hold on
h = boxplot(x(:,1),group,'Symbol','o','OutlierSize',3,'Widths',0.7,'Colors',[197,27,138]/255);
ylim([-0.4,-0.2]);
[r,p] = ttest2(x(group == 1,1),x(group == 2,1));
text(1.75,-0.23,['p = ' num2str(round(p,3))],'FontSize',6,'FontName','Helvetica','Color','k')
ylabel('Pearson Correlation','FontSize',7,'FontName','Helvetica','Color','k');
xlabel(['n = ',num2str(round(median(x(:,3)))),'/',num2str(length(ERprotein))],'FontSize',7,'FontName','Helvetica','Color','k');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
xticklabels({'AAC','MH34','B184'})
x1 = [1:0.001:2.5];
y1 = repmat(-0.25,length(x1),1);
plot(x1,y1,'-','LineWidth',0.75,'Color','k')
y1 = [-0.26:0.001:-0.25];
x1 = repmat(1,length(y1),1);
plot(x1,y1,'-','LineWidth',0.75,'Color','k')
x1 = repmat(2.5,length(y1),1);
plot(x1,y1,'-','LineWidth',0.75,'Color','k')
set(gca,'FontSize',6,'FontName','Helvetica');


figure
i = 3;
pax_abundance = num_abd(:,i);
pax_protein_list = raw_abd(2:end,1);
[~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
m = [];
m(idx~=0,1) = pax_abundance(idx(idx~=0));
m = m(1:400);
idx = find(m ~= 0 & ~isnan(m)&res_glc_final(:,1) ~= 0);
m = m(idx);
n = res_glc_final(idx,1);
h = scatter(log10(abs(1./n)),log10(m),20,'o','filled','LineWidth',0.75,'MarkerEdgeColor',[178, 24, 43]/255,'MarkerFaceColor',[178, 24, 43]/255,'MarkerFaceAlpha',0.3);
legend off;
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel(['mRNA abundace in log10 scale',char(13,10)','(fmol/mg DCW)'],'FontSize',7,'FontName','Helvetica','Color','k');
xlabel('Glucose cost in log10 scale','FontSize',7,'FontName','Helvetica','Color','k');
box on
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);

% [num_abd, raw_abd, ~] = xlsread('protein_abundance.xlsx','paxdb');
% pax_abundance = num_abd(:,5);
% pax_protein_list = raw_abd(3:end,14);
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0 & ~isnan(m)&res_glc_final(:,1) ~= 0)
% m = m(idx)
% n = res_glc_final(idx,1)
% [RHO4,PVAL4] = corr(log10(abs(n)),log10(m),'Type','Pearson')
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0&~isnan(m)&res_slope_final(:,1) ~= 0)
% m = m(idx)
% n = res_slope_final(idx,1)
% [RHO4,PVAL4] = corr(log10(m),log10(abs(n)),'Type','Pearson')
% 
% 
% [num_abd, raw_abd, ~]=xlsread('Proteome_collected.xlsx','mRNA_petri');
% pax_abundance = num_abd(:,5);
% pax_protein_list = raw_abd(2:end,1);
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0 & ~isnan(m)&res_glc_final(:,1) ~= 0)
% m = m(idx)
% n = res_glc_final(idx,1)
% [RHO4,PVAL4] = corr(log10(abs(n)),log10(m),'Type','Pearson')
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0&~isnan(m)&res_slope_final(:,1) ~= 0)
% m = m(idx)
% n = res_slope_final(idx,1)
% [RHO4,PVAL4] = corr(log10(m),log10(abs(n)),'Type','Pearson')

% [num_abd,raw_abd,~] = xlsread('Proteome_collected.xlsx','mRNA_Jianye2021')
% for i = 1:length(num_abd(1,:))
% pax_abundance = num_abd(:,i);
% pax_protein_list = raw_abd(2:end,1);
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0 & ~isnan(m)&res_glc_final(:,1) ~= 0)
% m = m(idx)
% n = res_glc_final(idx,1)
% [RHO4,PVAL4] = corr(log10(abs(1./n)),log10(m),'Type','Pearson')
% [~,idx] = ismember(strrep(ERprotein,'new',''),pax_protein_list);
% m = []
% m(idx~=0,1) = pax_abundance(idx(idx~=0));
% m = m(1:400)
% idx = find(m ~= 0&~isnan(m)&res_slope_final(:,1) ~= 0)
% m = m(idx)
% n = res_slope_final(idx,1)
% [RHO4,PVAL4] = corr(log10(abs(1./n)),log10(m),'Type','Pearson')
% x(i,1) = RHO4
% x(i,2) = PVAL4
% x(i,3) = length(m)
% end