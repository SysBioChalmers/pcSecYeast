%% buildModel
% Timing: ~  s

% Before run this script, please 1) check the annotation file is correct
% for all sce proteins
addpath(genpath('../../SecYeast/'))

% Matlab will coorupt if we put all information in one sheet
% [~,~,proteins] = xlsread('TableS1.xlsx','Annotation');
% [~,~,proteins_extra] = xlsread('TableS1.xlsx','Annotation_extra');
% proteins = proteins(:,1:15);
% proteins = [proteins;proteins_extra];
%load protein information
[~,~,protein_info] = xlsread('protein_information.xlsx');
protein_info_nomodification = protein_info;
protein_info_nomodification(2:end,3:9) = {0}; % this is the critial step for the withoutSec
clear proteins_extra 
tic;
%% Import yeast8
% cd ../ComplementaryData/Yeast8/;
% org_model = readCbModel('yeastGEM.xml');
% cd ../../;
% save('Yeast8.mat','org_model');
load('Yeast8.mat');

%% Modify the original model
[model_updated,energyResults,redoxResults] = modifyYeast8(org_model);
clear org_model;

%% Reformulate the orginal model
%   1. Reactions with isozymes (i.e., 'OR' case in GPR) will be copied,
%   each isozyme will be assigned to each copy. 
%   2. Reversible reactions will be split into forward and reverse ones.

model_split = splitModel(model_updated);
model = model_split;
clear model_updated model_split;

%% reformulate the Metabolic part
load('ProteinSequence.mat');
% Add complex formation reactions based on protein stoichiometry
% promiscuous = findPromiscuous(model_split);
load('Protein_stoichiometry.mat');% obtained from pdbe see Data for detail info
model = addComplexRxns(model,Protein_stoichiometry,protein_info_nomodification,ProteinSequence);

%% reformulate the Sec part
% Add Sec complexes
[~,~,proteins]=xlsread('TableS1.xlsx','Secretory');
model = addMachineryComplex(model,proteins,protein_info_nomodification,Protein_stoichiometry,ProteinSequence);

%% reformulate the ribosome/assembly factor/proteasome part
[~,~,protein_machinery]=xlsread('TableS1.xlsx','Machinery');
model = addMachineryComplex(model,protein_machinery,protein_info_nomodification,Protein_stoichiometry,ProteinSequence);

% manually update complex formation reactions for some complexes.
model = updateComplexformation(model,protein_info);

%% Collect kcats for enzymes

% use the one for with sec parameter
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('enzymedataDummyER.mat');
% enzymedata = collectkcats(model);
% 
% % manually update kcats for some reactions.
% enzymedata = updatekcats(enzymedata);
% % CONFIDENCE SCORE 5 means maually assigned kcats, not completely correct.
% 
% % get kdeg info from petri's proteome data
% enzymedata = getkdeg(enzymedata);
% enzymedata.kdeg(1:end) = 0.1;% using 0.1 for now
% 
% % Calculate molecular weight for each enzyme
% enzymedata = calculateMW(enzymedata,ProteinSequence,protein_info);
% enzymedataSEC = [];
% 
% % calculate the enzyme molecule weight of other machinery proteins such as
% % ribosome and proteasome
% enzymedataMachine = SimulateMachineParam(model);
% enzymedataMachine = calculateMW(enzymedataMachine,ProteinSequence,protein_info);

%% Add dummy complex reactions
% Dummy complex is assumed to be a part of metabolic protein pool.
% Note that the dummy complex is synthesized or diluted in the unit of
% mmol/gCDW/h.
% [protein_info,~] = SimulateDummyERParam(model,meanprotein_info,protein_info,enzymedataSEC);
protein_info = [protein_info;[{'dummyER'},{'dummyER'},0,0,0,0,0,0,0,'c',{''},1,{''},{''}]];
model = addDummyER(model,'r_4047','s_3717[c]',protein_info); 

% r_4047 is pseudo_protein_rxn_id in the GEM
% s_3717[c] is protein id
model = addDummy(model,'r_4047','s_3717[c]'); 
% close refolding and dilute for misfolding protein those reactions are for
% later usage
refold_list = model.rxns(contains(model.rxns,'refolding_'));
model = changeRxnBounds(model,refold_list,0,'b');
misfold_dilute_list = model.rxns(contains(model.rxns,'dilution_misfolding'));
model = changeRxnBounds(model,misfold_dilute_list,0,'b');

%% Change the original biomass equation
% Estimate modeled proteome
f_modeled_protein = estimateModeledprotein(model);
f_modeled_protein = round(f_modeled_protein,2); %g/gProtein
% Change the biomass equation
model = changeBiomass(model,f_modeled_protein,'r_4041','s_3717[c]');
% r_4041 is pseudo_biomass_rxn_id in the original GEM
% s_3717[c] is protein id in the original GEM
% unmodeled_cofactor[c] is newly added metabolite for unmodeled cofactor
% s_4206[c] is ion id in the original GEM
% s_4205[c] is cofactor id in the original GEM

% change GAM and NGAM
GAM= 36.5;
NGAM = 0.5;
model = changeGAM(model,GAM,NGAM);

%% Save model
save('pcSecYeastWithoutSEC.mat','model');
%% Save model to Excel
model_excel = model;
model_excel.subSystems = cell(length(model_excel.rxns),1);
model_excel.mets = model_excel.metNames;
writeCbModel(model_excel,'xls','pcSecYeastWithoutSec.xls');
clear model_excel;
toc;