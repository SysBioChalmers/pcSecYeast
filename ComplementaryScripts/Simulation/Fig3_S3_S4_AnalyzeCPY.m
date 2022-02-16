%% Figure 3A 3B
%% CPY
cd SimulateCPY_res
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');

model = setMedia(model,1);% minimal media (Delft media)
[model,enzymedataTP] = addTargetProtein(model,{'YMR297W'});

file = dir('*_5accauto.mat');
filename = {file.name};
idx = contains(filename,'_2_')|contains(filename,'_50_')|contains(filename,'_0_');
filename = filename(~idx);
fluxes = zeros(length(model.rxns),0);

for j = 1:length(filename)
    load(filename{j})
    if ~isempty(fluxes_simulated)
    fluxes(:,j) = fluxes_simulated;
    result(j,:) = res;
    end
end
result = result(any(fluxes,1),:);
fluxes = fluxes(:,any(fluxes,1));
result(:,6) = fluxes(strcmp(model.rxns,'YMR297W_dilution_misfolding_er'),:);
result(:,7) = fluxes(strcmp(model.rxns,'YMR297W_degradation_misfolding_c'),:);
data_abun = unique(result(:,2));
data_misfold = unique(result(:,3));
for i = 1:length(data_abun)
    idx_tmp = result(:,2) == data_abun(i);
    tmp = result(idx_tmp,:);
    [~,Idx] = ismember(data_misfold,tmp(:,3));
    data_mu(i,Idx~=0) = tmp(Idx(Idx~=0),1);
    data_acc(i,Idx~=0) = tmp(Idx(Idx~=0),6);
    data_deg(i,Idx~=0) = tmp(Idx(Idx~=0),7);
end

% 3d plot
TPrate = linspace(3.875E-7,1.9375e-04,3);

figure('Name','1');
surf(data_misfold,data_abun,data_mu,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xticks([0:0.2:1]);
yticks(TPrate)
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Misfold/total CPY ratio','FontSize',7,'FontName','Helvetica');
ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
zlabel(['Maximal growth rate',char(13,10)','(1/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

% misfoled CPY accumulation
figure('Name','2');
surf(data_misfold,data_abun,data_acc,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xticks([0:0.2:1]);
yticks(TPrate)
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Misfold/total CPY ratio','FontSize',7,'FontName','Helvetica');
ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
zlabel(['Accumulation misfolded ',char(13,10)', 'CPY rate[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);


% misfolded CPY degradation rate
figure('Name','3');
surf(data_misfold,data_abun,data_deg,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xticks([0:0.2:1]);
yticks(TPrate)
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Misfold/total CPY ratio','FontSize',7,'FontName','Helvetica');
ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
zlabel(['Degraded misfolded ',char(13,10)', 'CPY rate[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);
cd ../


%% growth reduction analysis CPY test for the Figure 4c
keys = {'2','50'}; % 1x and 25x
for k = 1:2
    key = keys{k};
file = dir(['*_',key,'_*.mat']);
filename = {file.name};
fluxes = zeros(length(model.rxns),0);
result= [];
for j = 1:length(filename)
    load(filename{j})
    if ~isempty(fluxes_simulated)
    fluxes(:,j) = fluxes_simulated;
    result(j,:) = res;
    end
end
% find the order as no misfold misfold 50% misfold 100% retention 10
% retention 20
order(1) = find(result(:,3) == 0);
order(2) = find(result(:,3) == 0.45);
order(3) = find(result(:,3) == 1 & result(:,4) == 5);
order(4) = find(result(:,3) == 1 & result(:,4) == 10);
order(5) = find(result(:,3) == 1 & result(:,4) == 20);
tmp = result(order,1);
fluxes = fluxes(:,order);
tmp = (0.3924-tmp)/0.3924;
h = bar(tmp,'FaceColor','flat','LineWidth',0.5,'FaceAlpha',0.3,'BarWidth',0.5);
xticklabels({'50%deg','100%deg','100%acc5','100%acc10','100%acc20'})
ylim([0 0.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Maximal growth reduction %','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('CPY expression with 0.1% of total proteome','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

% GSH GSSG supplementary figure
rxnName = 'transport GSSG from c to er';
[~,idx] = ismember(rxnName,model.rxnNames);
Trans_flux = fluxes(idx,:);
h = bar(sort(abs(Trans_flux)),'LineWidth',0.75,'FaceColor',[55,126,184]/255);
xticklabels({'50%deg','100%deg','100%acc5','100%acc10','100%acc20'})
ylim([0 6e-3])
set(gca,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
ylabel('Flux of transport GSSG(ER to Cytosol) [mmol/gDW/h]','FontSize',7,'FontName','Helvetica','Color','k');


% figure for Kar2P and Pdi1
protein_test = {'YJL034W','YCL043C'};
id = {'Kar2','Pdi1'};
% abun
[abun_complex_id,abun_protein_id,abun_complex,abun_protein,~,~] = abundanceCalculation(model,fluxes);

[~,idx] = ismember(protein_test,abun_protein_id);
abun_tmp = abun_protein(idx,:);
for i = 1:2
    figure
    h = bar(sort(abs(abun_tmp(i,:))),'LineWidth',0.75,'FaceColor',[55,126,184]/255);
    xticklabels({'50%deg','100%deg','100%acc5','100%acc10','100%acc20'})
    set(gca,'FontSize',6,'FontName','Helvetica');
    set(gcf,'position',[0 200 150 150]);
    set(gca,'position',[0.2 0.2 0.6 0.6]);
    ylabel([id{i},'Abundance [mmol/gDW]'],'FontSize',7,'FontName','Helvetica','Color','k');
    set(gcf,'position',[0 300 150 130]);
    set(gca,'position',[0.25 0.25 0.55 0.65]);
    ylim([0 8e-5])
end
clear result data_mu data_misfold data_acc data_abun data_deg fluxe
end
