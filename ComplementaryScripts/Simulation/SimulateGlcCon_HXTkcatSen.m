% sensitivityanalysisHXT for kcat

% test if HXT1 and HXT7 have the same kcat, then which one would be used

load('pcSecYeast.mat');
load('enzymedata.mat');
load('enzymedataSEC.mat');
load('enzymedataDummyER.mat');
load('enzymedataMachine.mat')

%% Set model
% set medium
model = setMedia(model,1);% minimal media (Delft media)
% set carbon source
model = changeRxnBounds(model,'r_1714',-1000,'l');% glucose
% set oxygen
model = changeRxnBounds(model,'r_1992',-1000,'l');
% block reactions
model = blockRxns(model);
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production

%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_unmodelER = tot_protein * 0.046;
clear tot_protein f_modeled_protein;
factor_k = 1;

hxt1rxnID = {'r_1166_10_complex'};
kcathxt7 = 200*3600; %/h
enzymedata.kcat(strcmp(enzymedata.enzyme,hxt1rxnID)) = kcathxt7;

mumax = 0.35;
enzymedata_all = CombineEnzymedata(enzymedata,enzymedataSEC,enzymedataMachine,enzymedataDummyER);

mu = mumax;
model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
disp(['mu = ' num2str(mu)]);
fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_all,factor_k,num2str(mu*100));

writeclusterfileLP(fileName,'sub_1')

% read the result
[sol_obj,sol_status,sol_full] = readSoplexResult([fileName,'.out'],model);

hxt = {'HXT1';'HXT2';'HXT3';'HXT4';'HXT7'};
hxtrxnID = {'r_1166_10';'r_1166_17';'r_1166_5';'r_1166_9';'r_1166_3'}; % glucose transporter rxn

[~,idx] = ismember(hxtrxnID,model.rxns);
hxt_flux = sol_full(idx,:);

figure

% heatmap figure
minclr = [255,255,255]/255;
maxclr = [202,0,32]/255;
tmp1 = linspace(minclr(1),maxclr(1),180)';
tmp2 = linspace(minclr(2),maxclr(2),180)';
tmp3 = linspace(minclr(3),maxclr(3),180)';
clrmap = [tmp1 tmp2 tmp3];

figure
label_tmp = cellfun(@num2str,num2cell(mu),'UniformOutput',false);
label_tmp = [1:1:length(hxt_flux)];
h_tmp = heatmap(hxt,'HXTs',hxt_flux','Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','none');
h_tmp.XLabel = 'Growth rate [1/h]';
set(h_tmp,'FontSize',6,'FontName','Helvetica');
h_tmp.FontColor = 'k';
set(gcf,'position',[0 200 150 20]);
set(gca,'position',[0.2 0.2 0.6 0.6]);


