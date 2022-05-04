% normal GEM simulation of amylase

% data from experiment Mingtao NC paper
glc = 10.83;% mmol/DCW/h data 1.9496 g/gDW/h
eth = 13.24; %  mmol/DCW/h data 0.61 g/gDW/h
gly = 0.173;% 0.016 g/gDW/h
ace = 0.796; % 0.047 g/gDW/h
mu = 0.31; % h-1

% load model
targetProteins = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'HGCSF'};
figure
hold on
for j = 1:length(targetProteins)
    j
load(['model',targetProteins{j},'.mat'])

% constrain the glc eth gly ace in the model; decrease the mu from 0.35 to
% 0.31
model = setMedia(model,1);% minimal media (Delft media)

model_tmp = changeRxnBounds(model,'r_1714',-glc,'l');
model_tmp = changeRxnBounds(model_tmp,'r_1761',eth,'l');
model_tmp = changeRxnBounds(model_tmp,'r_1634',ace,'l');
model_tmp = changeRxnBounds(model_tmp,'r_1808',gly,'l');
model_tmp = changeObjective(model_tmp,'r_2111',1);% a production
sol = optimizeCbModel(model_tmp);
mu_max = sol.f;
samplemu = linspace(0,mu_max,20);
model = changeRxnBounds(model,'r_1714',-glc,'l');
model = changeRxnBounds(model,'r_1761',eth,'l');
model = changeRxnBounds(model,'r_1634',ace,'l');
model = changeRxnBounds(model,'r_1808',gly,'l');
for i = 1:length(samplemu)
    model_tmp = changeRxnBounds(model,'r_2111',samplemu(i),'l');
    model_tmp = changeObjective(model_tmp,[targetProteins{j},' exchange'],1);% a production
    sol = optimizeCbModel(model_tmp);
    TP_res(j,i) = sol.f;
end
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
plot(samplemu,TP_res(j,:),'-','LineWidth',0.75,'Color',color(j,:))
end
 xlim([0 0.4]);
 set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);
box on
legend(targetProteins,'Fontsize',6)
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Protein production rate [mmol/gDW/h]','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate [1/h]','FontSize',7,'FontName','Helvetica');

 