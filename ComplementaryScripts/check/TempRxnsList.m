% This function is to generate all temp rxns
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('pcSecYeast.mat')
load('enzymedataDummyER.mat');

% three trans
peptide = 'xxx';
SP = 1;NG = 0;DSB = 0;GPI = 0; onlyrxns = 0;
[model,rxns{1}] = translocate(model,peptide,40,SP,NG,DSB,GPI,onlyrxns); %
SP = 0;NG = 1;DSB = 0;GPI = 0; onlyrxns = 0;
[model,rxns{2}] = translocate(model,peptide,40,SP,NG,DSB,GPI,onlyrxns);
SP = 0;NG = 0;DSB = 0;GPI = 0; onlyrxns = 0;
[model,rxns{3}] = translocate(model,peptide,40,SP,NG,DSB,GPI,onlyrxns);
DSB = 1;onlyrxns = 0;
[model,~,~] = addDSB(model,peptide,peptide,40,DSB,onlyrxns);
GPI = 1; onlyrxns = 0;
[model, ~,~] = addGPI(model,peptide,peptide,40,GPI,onlyrxns);
GPI = 1; onlyrxns = 0;
[model, ~,~] = addGPI(model,peptide,peptide,40,GPI,onlyrxns);
OG=1;onlyrxns = 0;
[model, ~,rxns{4}] = addOG(model,peptide,peptide,40,OG,onlyrxns);
NG=1;onlyrxns = 0;
[model, ~,rxns{5}] = addNG(model,peptide,peptide,40,NG,onlyrxns);
NG = 1;OG = 1;DSB = 1;GPI = 1;Trans = 1;compartment = 'c';
[model, ~,rxns{6}] = addMisfold(model,peptide,peptide,40,NG,OG,DSB,GPI,Trans,compartment,onlyrxns);
NG = 0;OG = 0;DSB = 0;GPI = 0;Trans = 1;compartment = 'c';
[model, ~,rxns{6}] = addMisfold(model,peptide,peptide,40,NG,OG,DSB,GPI,Trans,compartment,onlyrxns);
NG = 0;OG = 0;DSB = 0;GPI = 0;Trans = 1;compartment = 'er';
[model, ~,rxns{6}] = addMisfold(model,peptide,peptide,40,NG,OG,DSB,GPI,Trans,compartment,onlyrxns);
NG = 0;OG = 0;DSB = 0;GPI = 0;Trans = 0;compartment = 'er';
[model, ~,rxns{6}] = addMisfold(model,peptide,peptide,40,NG,OG,DSB,GPI,Trans,compartment,onlyrxns);
GPI=1;
[model, ~,rxns{9}] = coat_GPI(model,peptide,peptide,GPI,onlyrxns);
GPI=0;Trans=1;
[model, ~,rxns{10}] = coat_trans_membrane(model,peptide,peptide,GPI,Trans,onlyrxns);
GPI=0;Trans=0;
[model, ~,rxns{11}] = coat_other(model,peptide,peptide,GPI,Trans,onlyrxns);
NG = 1;
[model, ~,rxns{12}] = golgiProcessing_N(model,peptide,peptide,40,NG,onlyrxns);
OG = 1;
[model, ~,rxns{13}] = golgiProcessing_O(model,peptide,peptide,40,OG,onlyrxns);
[model, ~,rxns{14}] = mature(model,peptide,peptide,onlyrxns);
compartment = 'c';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'er';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'erm';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'g';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'gm';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'v';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'vm';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'e';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);
compartment = 'ce';
[model, ~,rxns{15}] = transportToFinal(model,peptide,peptide,compartment,onlyrxns);

% Extract rxns
idx = find(startsWith(model.rxns,'xxx'));
rxneq = printRxnFormula(model,'rxnAbbrList',model.rxns(idx),'metNameFlag',true);
[model.rxns(idx),rxneq]