% This to for Ploting amylase validation result

[num, txt, ~] = xlsread('Amylase.xlsx','ForFigure');
label_gene = txt(3:7,1);
Amylase_fold = num(1:5,:);
Amylase_mean_fold = mean(Amylase_fold,2);
Amylase_fold_std = std(Amylase_fold,1,2);
control = num([1,10],:);

figure('Name','Amylase');
base(:,1) = [1:1:5];
scatterdev = [base-0.3;base;base+0.3];

hold on;
b = bar(Amylase_mean_fold(1:5),0.7,'LineWidth',0.5,'FaceColor',[253,208,162]/255);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
errorbar(Amylase_mean_fold(1:5),Amylase_fold_std(1:5),'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,Amylase_fold(:),8,'k');
s.LineWidth = 0.5;
xticks([1:1:5]);
xlim([0.25 5.75]);
xticklabels(label_gene(1:5));
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 3.5]);
for i = 1:length(Amylase_fold(:,1))
    if ~strcmp(label_gene(i),'Control')
    if ismember(label_gene(i),{'PCM1';'PEP12';'VPS1';'OCH1';'GNA1';'CRS1';'USO1';'CYS4'})
        [~,p1] = ttest2(control(1,:),Amylase_fold(i,:),'Vartype','unequal');
    elseif ismember(label(i),{'IRE1';'SWA2';'ERV2';'MNS1';'ERO1';'SEC65'})
         [~,p1] = ttest2(control(2,:),Amylase_fold(i,:),'Vartype','unequal');
    end
    if p1 < 0.05 && p1 >= 0.01
        p1lable = '*';
        p1sz = 14;
    elseif p1 < 0.01 && p1 >= 0.001
        p1lable = '**';
        p1sz = 14;
    elseif p1 < 0.001
        p1lable = '***';
        p1sz = 14;
    else
        p1lable = 'n.s.';
        p1sz = 7;
    end
    text(i,Amylase_fold(i)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end



xlabel('Overexpression genes','FontSize',7,'FontName','Helvetica');
ylabel('Amylase yield fold','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 100 150 150]);
set(gca,'position',[0.1 0.2 0.75 0.5]);


label_gene = txt(8:18,1);
Amylase_fold = num(6:16,:);
Amylase_mean_fold = mean(Amylase_fold,2);
Amylase_fold_std = std(Amylase_fold,1,2);
control = num([1,6],:);

clear base scatterdev
figure('Name','Amylase');
base(:,1) = [1:1:11];
scatterdev = [base-0.3;base;base+0.3];

hold on;
b = bar(Amylase_mean_fold,0.7,'LineWidth',0.5,'FaceColor',[216,0,0]/255);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
errorbar(Amylase_mean_fold,Amylase_fold_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,Amylase_fold(:),8,'k');
s.LineWidth = 0.5;
xticks([1:1:11]);
xlim([0.25 11.75]);
xticklabels(label_gene);
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 3.5]);
for i = 1:length(Amylase_fold(:,1))
    if ~strcmp(label_gene(i),'Control')
    if ismember(label_gene(i),{'PCM1';'PEP12';'VPS1';'OCH1';'GNA1';'CRS1';'USO1';'CYS4'})
        [~,p1] = ttest2(control(1,:),Amylase_fold(i,:),'Vartype','unequal');
    elseif ismember(label_gene(i),{'IRE1';'SWA2';'ERV2';'MNS1';'ERO1';'SEC65'})
         [~,p1] = ttest2(control(2,:),Amylase_fold(i,:),'Vartype','unequal');
    end
    if p1 < 0.05 && p1 >= 0.01
        p1lable = '*';
        p1sz = 14;
    elseif p1 < 0.01 && p1 >= 0.001
        p1lable = '**';
        p1sz = 14;
    elseif p1 < 0.001
        p1lable = '***';
        p1sz = 14;
    else
        p1lable = 'n.s.';
        p1sz = 7;
    end
    text(i,Amylase_fold(i)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end



xlabel('Overexpression genes','FontSize',7,'FontName','Helvetica');
ylabel('Amylase yield fold','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 100 300 150]);
set(gca,'position',[0.1 0.2 0.75 0.5]);
box on

% 
% [num, txt, ~] = xlsread('Amylase.xlsx','Result');
% label_gene = txt(3:18,1);
% label_gene = strrep(label_gene,'1D-AAC+pSPGM1','Control');
% 
% OD = num(1:end,1:3);
% Amylase = num(1:end,4:6);
% Amylase_titer = num(1:end,7:9);
% 
% OD_mean = mean(OD,2);
% Amylase_mean = mean(Amylase,2);
% Amylase_titer_mean = mean(Amylase_titer,2);
% 
% OD_std = std(OD,1,2);
% Amylase_titer_std = std(Amylase_titer,1,2);
% 
% Amylase_fold(1:9,:) = Amylase(1:9,:)./Amylase_mean(1);
% Amylase_fold(10:16,:) = Amylase(10:16,:)./Amylase_mean(10);
% Amylase_mean_fold(1:9) = Amylase_mean(1:9)./Amylase_mean(1);
% Amylase_mean_fold(10:16) = Amylase_mean(10:16)./Amylase_mean(10);
% Amylase_fold_std = std(Amylase_fold,1,2);
% 
% 
% Amylase_titer_fold(1:9,:) = Amylase_titer(1:9,:)./Amylase_titer_mean(1);
% Amylase_titer_fold(10:16,:) = Amylase_titer(10:16,:)./Amylase_titer_mean(10);
% Amylase_titer_mean_fold(1:9) = Amylase_titer_mean(1:9)./Amylase_titer_mean(1);
% Amylase_titer_mean_fold(10:16) = Amylase_titer_mean(10:16)./Amylase_titer_mean(10);
% Amylase_titer_fold_std = std(Amylase_titer_fold,1,2);
% 
% % figure('Name','main');
% % hold on;
% % b = bar(OD_mean,0.7,'LineWidth',0.5,'FaceColor',[253,208,162]/255);
% % b.FaceColor = 'flat';
% % b.CData(1,:) = [255,255,255]/255;
% % b.CData(10,:) = [255,255,255]/255;
% % errorbar(OD_mean,OD_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
% % s = scatter(scatterdev,OD(:),8,'k');
% % s.LineWidth = 0.5;
% % xticks([1:1:16]);
% % xlim([0.25 16.75]);
% % xticklabels(label_gene);
% % set(gca,'FontSize',6,'FontName','Helvetica');
% % ylim([0 7.5]);
% % for i = 1:8
% % [~,p1] = ttest2(OD(1,:),OD(i+1,:),'Vartype','unequal');
% % if p1 < 0.05 && p1 >= 0.01
% %     p1lable = '*';
% %     p1sz = 14;
% % elseif p1 < 0.01 && p1 >= 0.001
% %     p1lable = '**';
% %     p1sz = 14;
% % elseif p1 < 0.001
% %     p1lable = '***';
% %     p1sz = 14;
% % else
% %     p1lable = 'n.s.';
% %     p1sz = 7;
% % end
% % text(i+1,OD(i+1)+0.8,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% % end
% % 
% % for i = 10:15
% % [~,p1] = ttest2(OD(10,:),OD(i+1,:),'Vartype','unequal');
% % if p1 < 0.05 && p1 >= 0.01
% %     p1lable = '*';
% %     p1sz = 14;
% % elseif p1 < 0.01 && p1 >= 0.001
% %     p1lable = '**';
% %     p1sz = 14;
% % elseif p1 < 0.001
% %     p1lable = '***';
% %     p1sz = 14;
% % else
% %     p1lable = 'n.s.';
% %     p1sz = 7;
% % end
% % text(i+1,OD(i+1)+0.8,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% % end
% % xlabel('Overexpression genes','FontSize',7,'FontName','Helvetica');
% % ylabel('OD 600','FontSize',7,'FontName','Helvetica');
% % set(gcf,'position',[200 100 300 150]);
% % set(gca,'position',[0.1 0.2 0.75 0.5]);
% 
% figure('Name','Amylase');
% 
% hold on;
% b = bar(Amylase_mean_fold,0.7,'LineWidth',0.5,'FaceColor',[253,208,162]/255);
% b.FaceColor = 'flat';
% b.CData(1,:) = [255,255,255]/255;
% b.CData(10,:) = [255,255,255]/255;
% errorbar(Amylase_mean_fold,Amylase_fold_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
% s = scatter(scatterdev,Amylase_fold(:),8,'k');
% s.LineWidth = 0.5;
% xticks([1:1:16]);
% xlim([0.25 16.75]);
% xticklabels(label_gene);
% set(gca,'FontSize',6,'FontName','Helvetica');
% ylim([0 3.5]);
% for i = 1:8
% [~,p1] = ttest2(Amylase_fold(1,:),Amylase_fold(i+1,:),'Vartype','unequal');
% if p1 < 0.05 && p1 >= 0.01
%     p1lable = '*';
%     p1sz = 14;
% elseif p1 < 0.01 && p1 >= 0.001
%     p1lable = '**';
%     p1sz = 14;
% elseif p1 < 0.001
%     p1lable = '***';
%     p1sz = 14;
% else
%     p1lable = 'n.s.';
%     p1sz = 7;
% end
% text(i+1,Amylase_fold(i+1)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% end
% 
% for i = 10:15
% [~,p1] = ttest2(Amylase_fold(10,:),Amylase_fold(i+1,:),'Vartype','unequal');
% if p1 < 0.05 && p1 >= 0.01
%     p1lable = '*';
%     p1sz = 14;
% elseif p1 < 0.01 && p1 >= 0.001
%     p1lable = '**';
%     p1sz = 14;
% elseif p1 < 0.001
%     p1lable = '***';
%     p1sz = 14;
% else
%     p1lable = 'n.s.';
%     p1sz = 7;
% end
% text(i+1,Amylase_fold(i+1)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% end
% 
% 
% xlabel('Overexpression genes','FontSize',7,'FontName','Helvetica');
% ylabel('Amylase yield fold','FontSize',7,'FontName','Helvetica');
% set(gcf,'position',[200 100 400 150]);
% set(gca,'position',[0.1 0.2 0.75 0.5]);

% figure('Name','amylasetiter');
% hold on;
% b = bar(Amylase_titer_mean_fold,0.7,'LineWidth',0.5,'FaceColor',[253,208,162]/255);
% b.FaceColor = 'flat';
% b.CData(1,:) = [255,255,255]/255;
% b.CData(10,:) = [255,255,255]/255;
% errorbar(Amylase_titer_mean_fold,Amylase_titer_fold_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
% s = scatter(scatterdev,Amylase_titer_fold(:),8,'k');
% s.LineWidth = 0.5;
% xticks([1:1:16]);
% xlim([0.25 16.75]);
% xticklabels(label_gene);
% set(gca,'FontSize',6,'FontName','Helvetica');
% ylim([0 3]);
% for i = 1:8
% [~,p1] = ttest2(Amylase_titer_fold(1,:),Amylase_titer_fold(i+1,:),'Vartype','unequal');
% if p1 < 0.05 && p1 >= 0.01
%     p1lable = '*';
%     p1sz = 14;
% elseif p1 < 0.01 && p1 >= 0.001
%     p1lable = '**';
%     p1sz = 14;
% elseif p1 < 0.001
%     p1lable = '***';
%     p1sz = 14;
% else
%     p1lable = 'n.s.';
%     p1sz = 7;
% end
% text(i+1,Amylase_titer_fold(i+1)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% end
% 
% for i = 10:15
% [~,p1] = ttest2(Amylase_titer_fold(10,:),Amylase_titer_fold(i+1,:),'Vartype','unequal');
% if p1 < 0.05 && p1 >= 0.01
%     p1lable = '*';
%     p1sz = 14;
% elseif p1 < 0.01 && p1 >= 0.001
%     p1lable = '**';
%     p1sz = 14;
% elseif p1 < 0.001
%     p1lable = '***';
%     p1sz = 14;
% else
%     p1lable = 'n.s.';
%     p1sz = 7;
% end
% text(i+1,Amylase_titer_fold(i+1)+0.5,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% end
% xlabel('Overexpression genes','FontSize',7,'FontName','Helvetica');
% ylabel('Amylase titer fold','FontSize',7,'FontName','Helvetica');
% set(gcf,'position',[200 100 300 150]);
% set(gca,'position',[0.1 0.2 0.75 0.5]);