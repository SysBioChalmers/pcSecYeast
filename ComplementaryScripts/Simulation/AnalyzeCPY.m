%% Figure 3A 3B
%% CPY
cd CPYres
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');

model = setMedia(model,1);% minimal media (Delft media)
[model,enzymedataTP] = addTargetProtein(model,{'YMR297W'});

file = dir('*300on.mat');
filename = {file.name};
fluxes = zeros(length(model.rxns),0);

for j = 1:length(filename)
    load(filename{j})
    fluxes(:,j) = fluxes_simulated;
    result(j,:) = res;
end
result(:,7) = fluxes(strcmp(model.rxns,'YMR297W_dilution_misfolding_er'),:);
data_abun = unique(result(:,2));
data_misfold = unique(result(:,3));
for i = 1:length(data_abun)
    idx_tmp = result(:,2) == data_abun(i);
    tmp = result(idx_tmp,:);
    [~,Idx] = ismember(data_misfold,tmp(:,3));
    data_mu(i,Idx~=0) = tmp(Idx(Idx~=0),1);
    data_acc(i,Idx~=0) = tmp(Idx(Idx~=0),7);
end

%% 3d plot
figure('Name','1');
surf(data_misfold,data_abun,data_mu,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xticks([0:0.2:1]);
TPrate = linspace(3.08E-7,1.925E-4,3);
yticks(TPrate)
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Misfold/total CPY ratio','FontSize',7,'FontName','Helvetica');
ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
zlabel(['Maximal growth rate',char(13,10)','(1/h)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

%% misfoled CPY accumulation
figure('Name','2');
surf(data_misfold,data_abun,data_acc,'EdgeColor','none');
hold on;
grid on;
az = -145;
el = 25;
view(az, el);
xticks([0:0.2:1]);
TPrate = linspace(3.08E-7,1.925E-4,3);
yticks(TPrate)
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Misfold/total CPY ratio','FontSize',7,'FontName','Helvetica');
ylabel(['CPY expression rate',char(13,10)','[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
zlabel(['Accumulation misfolded ',char(13,10)', 'CPY rate[mmol/gDW/h]'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[0 300 150 130]);
set(gca,'position',[0.25 0.25 0.55 0.65]);

cd ../
