% Figure 2 This function is to generate the component numbers for each
% process
load('pcSecYeast.mat')

index = {'DSB_sec_PDI_BIP_NEFS_complex','_ERNG_NG_sec_OSTC_complex','_OG_EROG_sec_Pmt2p_Pmt5p_Pmt1p_Pmt6p_Pmt4p_Pmt3p_complex','_GPIRI_sec_GPIR_complex',...
    '_LDSV_sec_Arf1p_Sec3p_Sec5p_Sec6p_Sec8p_Sec10p_Sec15p_Exo70p_Exo84p_Sec4p_Chc1p_Clc1p_complex','_HDSVI_sec_Arf1p_Pep12p_Swa2p_Chc1p_Clc1p_Apl4p_Apl2p_Apm1p_Aps1p_complex',...
    '_CPYI_sec_Gga1p_Gga2p_Arf1p_Apl4p_Apl2p_Apm1p_Aps1p_Chc1p_Clc1p_Pep12p_Vps45p_Vps5p_Swa2p_complex','_ALPtransport_sec_Apl6p_Aps3p_Apm3p_Apl5p_Vam3p_Clc1p_Chc1p_Arf1p_Swa2p_Vps1p_complex',...
    '_GLER_COPI_formation_sec_Arf1p_Gea1p_Gea2p_Rer1p_Erd2p_Cop1p_Sec26p_Sec27p_Sec21p_Ret2p_Sec28p_Ret3p_complex','COPII','_mature','_misfold_er'};
Components = {'DSB','NG','OG','GPI','LDSV','HDSV','CPY','ALP','COPI','COPII','Translocate','ERAD'};


for i = 1:length(Components)
    rxnidx = find(contains(model.rxns,index(i)));
    proteins = unique(extractBefore(model.rxns(rxnidx),index{i}));
    num(i) = length(proteins);
end

[num,a] = sort(num);
Components = Components(a);

figure('Name','main');
hold on;
b = bar(num,0.7,'LineWidth',0.5,'FaceColor',[253,208,162]/255);
s.LineWidth = 0.5;
xlim([0.25 12.75])
xticks([1:1:12]);
xticklabels(Components);
set(gca,'FontSize',6,'FontName','Helvetica');


xlabel('Secretory components','FontSize',7,'FontName','Helvetica');
ylabel('Protein processed','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 100 300 150]);
set(gca,'position',[0.1 0.2 0.75 0.5]);

