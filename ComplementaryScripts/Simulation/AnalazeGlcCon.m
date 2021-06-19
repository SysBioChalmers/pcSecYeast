%% plot figure 2A and Figure 2B

%% result for GlcCon
cd SimulateGlcCon
res = dir('res*.mat');
res = {res.name};

for i = 1:length(res)
    load(res{i});
    fluxes(:,i) = fluxes_simulated_without_sf;
end
    
load('pcSecYeast.mat')

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
[a,b] = sort(mu,'ascend');
fluxes = fluxes(:,b);
res = res(b);
res = strrep(res,'res','');
res = strrep(res,'.mat','');
mu = fluxes(strcmp(model.rxns,'r_2111'),:);

glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = -1*fluxes(strcmp(model.rxns,'r_1992'),:);
co2 = fluxes(strcmp(model.rxns,'r_1672'),:);

% Yeast exp data (PMID: 9603825)
fluxes_exp_yeast = [0.1  0.15  0.2  0.25  0.28  0.3   0.32  0.35  0.36  0.38 ; % mu
                    1.1  1.67  2.15 2.83  3.24  3.7   5.44  8.09  8.33  10.23; % glucose
                    0    0     0    0     0     0.51  4.42  6.91  6.71  14.91; % ethanol
                    2.73 2.5   5.07 6.8   8.3   8.8   6.83  6.6   7.1   4.19]; % o2

figure('Name','2');
hold on;
box on;

plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(2,:),'o','LineWidth',0.75,'Color',[55,126,184]/255,'MarkerSize',5);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(3,:),'o','LineWidth',0.75,'Color',[255,127,0]/255,'MarkerSize',5);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(4,:),'o','LineWidth',0.75,'Color',[77,175,74]/255,'MarkerSize',5);
plot(mu,glc,'-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu,etoh,'-','LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu,o2,'-','LineWidth',0.75,'Color',[77,175,74]/255);

xlim([0 0.4]);

set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
legend({'Glucose uptake',...
        'Ethanol production',...
        'O2 uptake'},'FontSize',6,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);


%% Plot the HXT change
hxt = {'HXT1';'HXT2';'HXT3';'HXT4';'HXT7'};
hxtrxnID = {'r_1166_10';'r_1166_17';'r_1166_5';'r_1166_9';'r_1166_3'}; % glucose transporter rxn
[~,idx] = ismember(hxtrxnID,model.rxns);
hxt_flux = fluxes(idx,:);
figure
hold on
box on
plot(mu,hxt_flux(5,:),'-','LineWidth',0.75,'Color',[55,126,184]/255)
plot(mu,hxt_flux(1,:),'-','LineWidth',0.75,'Color',[189,0,81]/255)

xticks([0:0.1:0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
ylabel(['Glucose transport fluxes',char(13,10)','(mmol/gDW/h)'],'FontSize',7,'FontName','Helvetica','Color','k');
legend({'HXT7','HXT1'},'FontSize',6,'FontName','Helvetica','location','nw')
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);