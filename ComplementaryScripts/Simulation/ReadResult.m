%% crabtree

load('pcSecYeast.mat');
cd crabtree_res/;
file = dir('*.out');
filename = {file.name};
fluxes = zeros(length(model.rxns),0);
for j = 1:length(filename)
    [~,sol_status,sol_full] = readSoplexResult(filename{j},model);
    if strcmp(sol_status,'optimal')
        fluxes = [fluxes sol_full];
    else
        fluxes = [fluxes zeros(length(model.rxns),1)];
    end
end
save('fluxes_Crabtree.mat','fluxes');

cd ../

%% TP1
TP = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'HGCSF'};
load('pcSecYeast.mat');
cd SimulateTP_SDAA_res/;
file = dir('*.out');
filename = {file.name};
for i = 1:length(TP)
    display([num2str(i) '/' num2str(length(TP))]);
    load(['model',TP{i},'.mat']);
    
    fluxes = zeros(length(model.rxns),0);
    allfile = filename(startsWith(filename,['Simulation_dilution',TP{i},'_']));
    for j = 1:length(allfile)
        [~,sol_status,sol_full] = readSoplexResult(allfile{j},model);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
        else
            fluxes = [fluxes zeros(length(model.rxns),1)];
        end
    end
    save(['fluxes_',TP{i},'.mat'],'fluxes');
end
% ref result
file = dir('*.out');
filename = {file.name};
load('pcSecYeast.mat');
fluxes = zeros(length(model.rxns),0);
allfile = filename(startsWith(filename,'Simulation_dilutionref_'));
for j = 1:length(allfile)
    [~,sol_status,sol_full] = readSoplexResult(allfile{j},model);
    if strcmp(sol_status,'optimal')
        fluxes = [fluxes sol_full];
    else
        fluxes = [fluxes zeros(length(model.rxns),1)];
    end
end
save('fluxes_ref.mat','fluxes');

cd ../


%% ProteinCost
[~,~,protein_info] = xlsread('../../ComplementaryData/protein_information.xlsx');
protein_info = protein_info(2:end,:);
ERprotein = protein_info(cell2mat(protein_info(:,3)) == 1,2); % go through ER
targetProteins = [ERprotein(1:3);{'Amylase';'Insulin'}]
load('pcSecYeast.mat');
cd SimulateProteinCostResult/;
file = dir('*.out');
filename = {file.name};
for i = 1:length(targetProteins)
    display([num2str(i) '/' num2str(length(targetProteins))]);
   [model_tmp,enzymedataTP] = addTargetProtein(model,targetProteins(i));
    
    fluxes = zeros(length(model_tmp.rxns),0);
    allfile = filename(contains(filename,targetProteins{i}));
    for j = 1:length(allfile)
        [~,sol_status,sol_full] = readSoplexResult(allfile{j},model_tmp);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
        else
            fluxes = [fluxes zeros(length(model_tmp.rxns),1)];
        end
    end
    save(['../FluxesProteinCost/fluxes_',targetProteins{i},'.mat'],'fluxes');
end
fluxes = zeros(length(model.rxns),0);
    allfile = filename(contains(filename,'ref'));
    for j = 1:length(allfile)
        [~,sol_status,sol_full] = readSoplexResult(allfile{j},model);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
        else
            fluxes = [fluxes zeros(length(model.rxns),1)];
        end
    end
    save(['../FluxesProteinCost/fluxes_ref.mat'],'fluxes');
cd ../


%% Crabtree effect

load('pcSecYeast.mat');
mkdir('FluxesCrabtree')

cd crabtreeresult/;
file = dir('*.out');
filename = {file.name};
    fluxes = zeros(length(model.rxns),0);

    for j = 1:length(filename)
        [~,sol_status,sol_full] = readSoplexResult(filename{j},model);
        if strcmp(sol_status,'optimal')
            fluxes = [fluxes sol_full];
        else
            fluxes = [fluxes zeros(length(model.rxns),1)];
        end
    end
    save(['../FluxesCrabtree/fluxes_Crabtree.mat'],'fluxes');

cd ../

%% FakeTP
k = 1:50:1112;
for i = 1:length(k)-1
    load(['res_maxTP',num2str(k(i)),'.mat'])
    res(k(i):k(i)+49,1) = maxTP(k(i):k(i)+49);
end

load(['res_maxTP1051.mat'])
res(k(i):1112,1) = maxTP(k(i):1112);
load('fakeProteinInfo.mat')
fakeProteinInfo(:,15) = num2cell(res);
fakeProteinInfo(find(res == 0),:) = [];
t = cell2table(fakeProteinInfo);
writetable(t,'res_FakeProtein.txt','Delimiter','\t','QuoteStrings',false,'WriteRowNames',true)
    