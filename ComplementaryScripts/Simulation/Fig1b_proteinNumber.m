% Figure 1b
load('pcSecYeast.mat')
[~,~,proteins_info]=xlsread('TableS1.xlsx','Secretory');
proteins_info = proteins_info(2:end,[1,2,26]);
processes = unique(proteins_info(:,3));
for i = 1:length(processes)
    res(i,1) = sum(strcmp(proteins_info(:,3),processes(i,1))); % machinery proteins
    complex = proteins_info(strcmp(proteins_info(:,3),processes(i,1)),1);
    proteins = model.rxns(endsWith(model.rxns,complex));
    proteins = extractBefore(proteins,'_');
    proteins = setdiff(proteins,{'dummyER','dummy'});
    res(i,2) = length(unique(proteins));
end
% sort the number
[a,b] = sort(sum(res,2));
processes = processes(b);
res = res(b,:);
h = bar(res,'stacked','LineWidth',0.7,'BarWidth',0.5,'FaceAlpha',0.3);
xticklabels(processes)
set(gca,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
ylabel('Protein number','FontSize',7,'FontName','Helvetica','Color','k');


