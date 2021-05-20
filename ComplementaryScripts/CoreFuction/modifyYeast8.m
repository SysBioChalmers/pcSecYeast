%% addCofactorRxns
function [model,energyResults,redoxResults] = modifyYeast8(model)

% %Load yeast model:
% initCobraToolbox
% git('clone https://github.com/SysBioChalmers/yeast-GEM.git')
% cd yeast-GEM/ComplementaryScripts/
% model = loadYeastModel;


% Add new reactions and metabolites
[~,rxnlst,~] = xlsread('Yeast8_Modification.xlsx','New reactions');
[~,~,metlst] = xlsread('Yeast8_Modification.xlsx','New metabolites');
rxnlst = rxnlst(2:end,:);
metlst =metlst(2:end,:);
% add rxns
for i = 1:length(rxnlst(:,1))
    
    
    rxnformula = rxnlst{i,3};
    [stoichCoeffList, metaboliteList, ~, revFlag]=constructS({rxnformula});
    comps = split(metaboliteList, ' [');
    if length(comps(1,:)) == 1
        comps = comps(2);
    else
        comps = comps(:,2);
    end
    comps = strrep(comps,']','');
    CONValldata = cat(2,model.compNames,model.comps);
    [~,b] = ismember(comps,CONValldata(:,1));
    comps = CONValldata(b,2);
    
    %mapping mets to model.metnames, get s_ index for new mets
    for j = 1:length(metaboliteList)
        [~,metindex] = ismember(metaboliteList(j),model.metNames);
        if metindex ~= 0
            mets(j) = model.mets(metindex);
        elseif metindex == 0
            newID = getNewIndex(model.mets);
            mets(j) = strcat('s_',newID,'[',comps(j),']');
            model = addMetabolite(model,char(mets(j)), ...
                'metName',metaboliteList{j});
        end
    end
    
    newID = getNewIndex(model.rxns);
    [model, rxnIDexists] = addReaction(model,...
        ['r_',newID],...
        'reactionName', rxnlst{i,2},...
        'metaboliteList',mets,...
        'stoichCoeffList',stoichCoeffList,...
        'reversible',revFlag,...
        'geneRule',rxnlst{i,4},...
        'notes',rxnlst{i,5},...
        'checkDuplicate',1);
    if isempty(rxnIDexists)
        newrxns(i) = {['r_',newID]};
        [energyResults(i,:),redoxResults(i,:)] = CheckEnergyProduction(model);
    else
        newrxns(i) = model.rxns(rxnIDexists);
        energyResults(i,:) = {'alreadly exists','skip','skip'};
        redoxResults(i,:) = {'alreadly exists','skip','skip'};
    end
    clear mets 
end

% add rxn notes 
[~,idx] = ismember(newrxns,model.rxns);
model.rxnNotes(idx(idx~=0)) = rxnlst((idx~=0),5);


% add new metabolite new annotation
[~,idx] = ismember(metlst(:,1),model.metNames);
model.metFormulas(idx(idx~=0)) = metlst(idx~=0,3);
model.metCharges(idx(idx~=0)) = cell2mat(metlst(idx~=0,4));
model.metChEBIID(idx(idx~=0)) = metlst(idx~=0,5);
model.metMetaNetXID(idx(idx~=0)) = metlst(idx~=0,6);

%% Add gene standard name for new genes
[~,genelst,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
genelst = genelst(2:end,:);
allgene{1} = genelst(:,1);
allgene{2} = genelst(:,2);
for i = 1: length(model.genes)
    geneIndex = strcmp(allgene{1}, model.genes{i});
    if sum(geneIndex) == 1 && ~isempty(allgene{2}{geneIndex})
        model.geneNames{i} = allgene{2}{geneIndex};
    else
        model.geneNames{i} = model.genes{i};
    end
end

% Add protein name for genes
for i = 1:length(model.genes)
    model.proteins{i} = strcat('COBRAProtein',num2str(i));
end

model = rmfield(model,'grRules');

%% check whether the rxn is mass balanceed/charge balanced
% [MassChargeresults] = CheckBalanceforSce(model,newrxns);

% print rxns
printRxnFormula(model,'rxnAbbrList',newrxns,'metNameFlag',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newID = getNewIndex(IDs)

%Find latest index and create next one:
IDs   = regexprep(IDs,'[^(\d*)]','');
IDs   = str2double(IDs);
newID = max(IDs) + 1;
newID = num2str(newID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [energyResults,redoxResults] = CheckEnergyProduction(model)
% This function is to  check whether adding new reaction will lead to infinite
% ATP and NADH production
% the input is model and new reaction ID lists
% output is a result with pass or error
%
% Feiran Li
energyResults     = {};
redoxResults = {};
    model_test = model;
    model_test = minimal_Y6(model_test);
    %Add/change ATP production reaction:
    %            ATP    +    H2O    ->  ADP     +   H+      +  PO4
    mets  = {'s_0434[c]','s_0803[c]','s_0394[c]','s_0794[c]','s_1322[c]'};
    coefs = [-1,-1,1,1,1];
    model_test = addReaction(model_test,{'GenerateATP','leaktest1'}, ...
                             mets,coefs,false,0,1000);
    model_test = changeObjective(model_test,'GenerateATP', 1);
    sol = optimizeCbModel(model_test);
    if sol.f <= 360 && sol.f > 0 %later can be changed to the experimental value
        energyResults = [energyResults; model.rxns(end),'pass',num2str(sol.f)];
    elseif sol.f > 360
        energyResults = [energyResults; model.rxns(end),'Fail',num2str(sol.f)];
    else
        energyResults = [energyResults; model.rxns(end),'error','error'];
    end
    
    model_test = model;
    model_test = minimal_Y6(model_test);
    %Add/change NADH production reaction:
    %            NADH[c] + H[c] =>  NAD[c]
    mets  = {'s_1203[c]','s_0794[c]','s_1198[c]'};
    coefs = [-1,-1,1];
    model_test = addReaction(model_test,{'GenerateNADH','leaktest12'}, ...
                             mets,coefs,false,0,1000);
    model_test = changeObjective(model_test, model_test.rxns(end), 1);
    sol = optimizeCbModel(model_test);
    if sol.f <= 120 && sol.f > 0 %later can be changed to the experimental value
        redoxResults = [redoxResults; model.rxns(end),'pass',num2str(sol.f)];
    elseif sol.f > 120
        redoxResults = [redoxResults; model.rxns(end),'Fail',num2str(sol.f)];
    elseif sol.f <= 0
        redoxResults = [redoxResults; model.rxns(end),'error','error'];
    end

end

function model = minimal_Y6(model)
% change Y6 model media to minimal - ammonium, glucose, oxygen,
% phosphate, sulphate
% Bicarbonate production is blocked to get rid of conflict due to bicarbonate and carbon dioxide
% are regarded as equivalent
% the function is from:https://doi.org/10.1371/journal.pcbi.1004530

% start with a clean slate: set all exchange reactions to upper bound = 1000
% and lower bound = 0 (ie, unconstrained excretion, no uptake)

%Feiran Li 2018.09.05

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

desiredExchanges = {'r_1654'; ... % ammonium exchange
                    'r_1992'; ... % oxygen exchange
                    'r_2005'; ... % phosphate exchange
                    'r_2060'; ... % sulphate exchange
                    'r_1861'; ... % iron exchange, for test of expanded biomass def
                    'r_1832'; ... % hydrogen exchange
                    'r_2100'; ... % water exchange
                    'r_4593'; ... % chloride exchange
                    'r_4595'; ... % Mn(2+) exchange
                    'r_4596'; ... % Zn(2+) exchange
                    'r_4597'; ... % Mg(2+) exchange
                    'r_2049'; ... % sodium exchange
                    'r_4594'; ... % Cu(2+) exchange
                    'r_4600'; ... % Ca(2+) exchange
                    'r_2020' };   % potassium exchange

blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange

glucoseExchange = {'r_1714'};     % D-glucose exchange

uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges);
glucoseExchangeIndex = findRxnIDs(model,glucoseExchange);
BlockedRxnIndex      = findRxnIDs(model,blockedExchanges);

if length(find(uptakeRxnIndexes~= 0)) ~= 15
    warning('Not all exchange reactions were found.')
end


model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0))     = -1000;
model.lb(glucoseExchangeIndex) = -1;

model.lb(BlockedRxnIndex) = 0;
model.ub(BlockedRxnIndex) = 0;

end