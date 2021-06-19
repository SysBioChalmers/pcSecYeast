%% Targets identification

targetProteins = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'Rancomplex';'HBsAg';'HGCSF'};

cd FluxesTP1/;
for i = 1:length(targetProteins)
        display([num2str(i) '/' num2str(length(targetProteins))]);
  
    load(['../SimulateTPres0618/model',targetProteins{i},'.mat']);
     load(['fluxes_',targetProteins{i},'.mat']);
    tp_ex = [targetProteins{i},' exchange'];
    
    mu = round(fluxes(ismember(model.rxns,'r_2111'),:),2);
    [a,b] = sort(mu,'ascend');
    fluxes = fluxes(:,b);
    %fluxes = fluxes(:,mu > 0.3);
    
    nonzeroidx = any(fluxes);
    fluxes = fluxes(:,nonzeroidx);
    [abun_complex_id,abun_protein_id,abun_complex,abun_protein,~,~] = abundanceCalculation(model,fluxes);
    q_tp = fluxes(ismember(model.rxns,tp_ex),:);
    idx_max = find(q_tp == max(q_tp));
    idx_min = find(q_tp == min(q_tp));
    
    res.genelist = cell(0,1);
    res.r = zeros(0,1);
    res.p = zeros(0,1);
    res.slope = zeros(0,1);
    res.i = zeros(0,1);
    res.ratio = zeros(0,1);
    res.pro_abun = zeros(0,1);
    for m = 1:length(abun_protein_id)
        if ~any(abun_protein(m,:) == 0)
            geneid = abun_protein_id(m);
            res.genelist = [res.genelist;geneid];
            [RHOtmp,PVALtmp] = corr(abun_protein(m,:)',q_tp','Type','Spearman');
            res.r = [res.r;RHOtmp];
            res.p = [res.p;PVALtmp];
            p = polyfit(abun_protein(m,:),q_tp,1);
            res.slope = [res.slope;p(1)];
            res.i = [res.i;m];
            x = abun_protein(m,:);
            res.pro_abun = [res.pro_abun;x(idx_max + 1)];
            res.ratio = [res.ratio;x(idx_max+1)/x(idx_min -1)];
        end
    end
    
    % remove ribosome genes
    load('enzymedataMachine.mat')
    ribosomegene = unique(enzymedataMachine.subunit(1:2,:));
    ribosomegene = setdiff(ribosomegene,'');
    [~,idx] = ismember(ribosomegene,res.genelist);
    res.genelist(idx) = [];
    res.r(idx) = [];
    res.p(idx) = [];
    res.i(idx) = [];
    res.ratio(idx) = [];
    res.pro_abun(idx) = [];
    res.slope(idx) = [];
    % annotation and protein abundace ref from paxdb
    [num_abd, raw_abd, ~] = xlsread('protein_abundance.xlsx','paxdb');
    res.ref_abun = zeros(length(res.genelist),1);
    [~,idx] = ismember(res.genelist,raw_abd(:,2));
    res.ref_abun(idx~=0) = num_abd(idx(idx~=0),3);
    
    [~, ~, raw_anno] = xlsread('Protein_annotation_uniprot.xlsx','uniprot');
    res.annot = cell(length(res.genelist),8);
    [~,idx] = ismember(res.genelist,raw_anno(:,1));
    res.annot(idx~=0,:) = raw_anno(idx(idx~=0),2:end);
    
    % annotation of homologous gene
    [~, ~, raw_homo] = xlsread('Yeast_homologue_genes.xlsx');
    res.homo = cell(length(res.genelist),1);
    [~,idx] = ismember(res.genelist,raw_homo(:,2));
    res.homo(idx~=0,:) = raw_homo(idx(idx~=0),5);
    
    % annotation of protein wther it is metabolic
    res.met = repmat({'metabolic'},length(res.genelist),1);
    load('enzymedataSEC.mat')
    secprotein = unique(enzymedataSEC.subunit);
    [~,idx] = ismember(res.genelist,secprotein);
    res.met(idx ~= 0,1) = {'Secretory'};
    
    % annotation of complex or not
    [~, ~, raw_complex] = xlsread('Protein_annotation_uniprot.xlsx','complex_portal');
    idx = contains(raw_complex(:,20),'|');
    raw_complex = raw_complex(idx,:);
    complex_protein = [];
    res.complex = cell(length(res.genelist),1);
    for j = 1:length(raw_complex)
        complex_protein = split(raw_complex(j,20),'|');
        [~,idx] = ismember(res.genelist,complex_protein);
        res.complex(idx~=0) = raw_complex(j,2);
    end
    
    % difference between protein abundance/paxdb abundance
    res.difference = res.pro_abun./res.ref_abun;
    
    res.priority = zeros(length(res.genelist),1);

    % priority
    idx = find(res.r > 0.85 & res.ratio > 1);
    res.priority(idx) = 4;
    idx = find(res.r > 0.85 & res.ratio > 2);
    res.priority(idx) = 3;
    idx = find(res.r > 0.85 & res.ratio > 2 & res.difference > 1E-10);
    res.priority(idx) = 2;
    idx = find(res.r > 0.85 & res.ratio > 2 & res.difference > 1E-10 & cellfun(@isempty,res.homo)& cellfun(@isempty,res.complex));
    res.priority(idx) = 1;
    res_list = cell('');
    t = table(res.genelist,res.priority,res.r,res.p,res.i,res.ratio,res.slope,res.pro_abun,res.ref_abun,res.difference,res.homo,res.met,res.complex,res.annot(:,1),res.annot(:,2),res.annot(:,3),res.annot(:,4),res.annot(:,5),res.annot(:,6),res.annot(:,7),res.annot(:,8),'VariableNames',{'protein' ,'priority','correlation', 'p_value','index_in_Model','ratio','slope','prot_abun(max_production)','ref_abun(paxdb)','abun_fold','homolougous_gene','metabolic/secretoy','complex/not',raw_anno{1,2:9}});
    res_list(1:length(res.genelist(res.priority == 1 |res.priority == 2)),i) = res.genelist(res.priority == 1 |res.priority == 2);
    res_listAnnotation(1:length(res.genelist(res.priority == 1 |res.priority == 2)),i) = res.annot(res.priority == 1 |res.priority == 2,3);
    res_listFunc(1:length(res.genelist(res.priority == 1 |res.priority == 2)),i) = res.met(res.priority == 1 |res.priority == 2);
    writetable(t,[targetProteins{i},'Target.txt'],'Delimiter','\t','QuoteStrings',false,'WriteRowNames',true)
end
save('res_geneList.mat','res_list','targetProteins');
cd ../

%% All targets
% res_list is to serve as input for DiVenn
res_allgene = res_list(:);
res_allgene(cellfun(@isempty,res_allgene)) = [];
m = tabulate(res_allgene);
m(:,4) = repmat({'metabolic'},length(m(:,1)),1);
load('enzymedataSEC.mat')
secprotein = unique(enzymedataSEC.subunit);
[~,idx] = ismember(m(:,1),secprotein);
m(idx ~= 0,4) = {'Secretory'};
[~,~,proteins]=xlsread('TableS1.xlsx','Secretory');
[~,idx] = ismember(m(:,1),proteins(:,2));
m(idx ~= 0,4) = proteins(idx(idx~=0),26);
cata = unique(m(:,4));
for i = 1:length(targetProteins)
    idx = find(cell2mat(m(:,2)) == i);
    tmp = m(idx,:);
    tmp_cata = tabulate(tmp(:,4));
    [~,idx] = ismember(tmp_cata(:,1),cata);
    cata_res(i,idx) = cell2mat(tmp_cata(:,2));
end
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
colorbar(color)
bar(cata_res,'stacked','LineWidth',0.5,'BarWidth',0.5,color);
legend(cata,'Fontsize',6)
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Overexpression target number','FontSize',7,'FontName','Helvetica');
xlabel('recombinant protein number','FontSize',7,'FontName','Helvetica');
set(gca,'position',[0.2 0.2 0.6 0.6]);
set(gcf,'position',[0 200 150 150]);
box on;

%% Protein secretion
targetProteins = {'Amylase';'Insulin';'HSA';'Hemoglobincomplex';'Humantransferin';'BGL';'PHO';'Rancomplex';'HBsAg';'HGCSF'};

cd FluxesTP1/;
color = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]./255;
figure
hold on
for i = 1:length(targetProteins)
display([num2str(i) '/' num2str(length(targetProteins))]);
load(['../SimulateTPres0618/model',targetProteins{i},'.mat']);
load(['fluxes_',targetProteins{i},'.mat']);
nonzeroidx = any(fluxes);
fluxes = fluxes(:,nonzeroidx);
tp_ex = [targetProteins{i},' exchange'];
mu = round(fluxes(ismember(model.rxns,'r_2111'),:),2);
[mu,b] = sort(mu,'ascend');
fluxes = fluxes(:,b);

q_tp = fluxes(ismember(model.rxns,tp_ex),:);
plot(mu,q_tp,'-','LineWidth',0.75,'Color',color(i,:))
end

legend(targetProteins,'Fontsize',6)
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Protein production rate [mmol/gDW/h]','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate [1/h]','FontSize',7,'FontName','Helvetica');
