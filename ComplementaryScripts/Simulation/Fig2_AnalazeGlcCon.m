%% plot figure 2A and Figure 2B

%% result for GlcCon
cd SimulateGlcCon_res/
res = dir('res*.mat');
res = {res.name};

for i = 1:length(res)
    load(res{i});
    fluxes(:,i) = fluxes_simulated_without_sf;
end
    
load('pcSecYeast.mat')

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
res = strrep(res,'res','');
res = strrep(res,'.mat','');
res = cellfun(@str2num,res);

[a,b] = sort(res,'ascend');
fluxes = fluxes(:,b);
res = res(b);

fluxes = fluxes(:,res <= 100);
res = res(res <= 100);

mu = fluxes(strcmp(model.rxns,'r_2111'),:);

glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = -1*fluxes(strcmp(model.rxns,'r_1992'),:);
co2 = fluxes(strcmp(model.rxns,'r_1672'),:);

%% Figure 2a 
figure('Name','1');
hold on;
hold on;
yyaxis left;
plot(log10(res),glc,'.-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(log10(res),etoh,'.-','LineWidth',0.75,'Color',[255,127,0]/255);
ylabel('Flux [mmol/gCDW/h]','FontSize',7,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'ycolor','k');
ylim([0 40]);

yyaxis right;
plot(log10(res),mu,'.-','LineWidth',0.75,'Color',[186,186,186]/255);
set(gca,'ycolor','k');
ylim([0 0.45]);
ylabel('growth rate [1/h]','FontSize',7,'FontName','Helvetica');

legend({'qglucose',...
        'qethanol','Growth rate'},'FontSize',6,'FontName','Helvetica','location','nw');

xlabel('Extracelluar glucose in log10 scale [mM]','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[0 200 150 100]);
set(gca,'position',[0.2 0.2 0.6 0.6]);


%% Yeast exp data (PMID: 9603825)
fluxes_exp_yeast = [0.1  0.15  0.2  0.25  0.28  0.3   0.32  0.35  0.36  0.38 ; % mu
                    1.1  1.67  2.15 2.83  3.24  3.7   5.44  8.09  8.33  10.23; % glucose
                    0    0     0    0     0     0.51  4.42  6.91  6.71  14.91; % ethanol
                    2.73 2.5   5.07 6.8   8.3   8.8   6.83  6.6   7.1   4.19]; % o2

figure('Name','2');
hold on;

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
set(gcf,'position',[0 200 150 100]);
set(gca,'position',[0.2 0.2 0.6 0.6]);


%% Plot the HXT change Figure 2B
hxt = {'HXT1';'HXT2';'HXT3';'HXT4';'HXT7'};
hxtrxnID = {'r_1166_10';'r_1166_17';'r_1166_5';'r_1166_9';'r_1166_3'}; % glucose transporter rxn

[~,idx] = ismember(hxtrxnID,model.rxns);
hxt_flux = fluxes(idx,:);
hxt = hxt(any(hxt_flux,2)); % remove hxt without flux
hxt_flux = hxt_flux(any(hxt_flux,2),:);
figure

figure
hold on
color = [202,0,32
244,165,130
146,197,222
5,113,176]./255;
for i = 1:length(hxt)
plot(log10(res),hxt_flux(i,:),'.-','LineWidth',0.75,'Color',color(i,:));
end
box on

set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Extracelluar glucose in log10 scale [mM]','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Glucose uptake rate (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
legend({'Hxt1','Hxt3','Hxt7'},'FontSize',6,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
%% plot the MM figure
hxt = {'Hxt1';'Hxt2';'Hxt3';'Hxt4';'Hxt7'};
hxtrxnID = {'r_1166_10';'r_1166_17';'r_1166_5';'r_1166_9';'r_1166_3'}; % glucose transporter rxn
Km_hxt = [110;1.5;34;9.3;2.5]; % bionumber 110954 110739 PMID 2482015 10336421 https://doi.org/10.1101/2020.06.22.165753 table
hxt = {'YHR094C';'YMR011W';'YDR345C';'YHR092C';'YDR342C'};
kcatmax = [1012 53 479 155 197];
glccost = [3002.6 3145.6
2829.9 2963.2
3082.7 3232.0
3007.0 3149.1
3072.0 3220.5];
[~,idx] = ismember(hxtrxnID,model.rxns);
hxt_flux = fluxes(idx,:);
hxt = hxt(any(hxt_flux,2)); % remove hxt without flux
Km_hxt = Km_hxt(any(hxt_flux,2)); % remove hxt without flux
kcatmax = kcatmax(any(hxt_flux,2)); % remove hxt without flux
hxt_flux = hxt_flux(any(hxt_flux,2),:);
hxt_flux = sum(hxt_flux,1);

for i = 1:length(hxt)
glcCost(i,:) = hxt_flux./(kcatmax(i).*(res./(Km_hxt(i) + res)))*glccost(i);
end

figure 
hold on
color = jet(10); % set colorcode
for i = 1:3
m = plot(log10(res),log10(glcCost(i,:)),'.-','LineWidth',0.75,'Color',color(2*i,:));
end
set(gca,'FontSize',6,'FontName','Helvetica');

xlabel('Residue glucose concentration [mM] log10','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Glucose cost for sustain certain uptake','FontSize',7,'FontName','Helvetica','Color','k');
box on
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.4 0.6 0.6]);
