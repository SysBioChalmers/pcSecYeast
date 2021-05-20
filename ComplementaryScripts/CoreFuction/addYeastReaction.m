function model=addYeastReaction(model,eq,rxnID,rxnName,rule)

if nargin < 5
    rule = '';
end

if iscell(rxnID)
    rxnID = cell2mat(rxnID);
end
%find and add new metabolites and new compartements to secModel
[S, mets, ~, reversible]=constructS({eq});
CONValldata = cat(2,strcat(repmat({' ['},length(model.compNames),1),model.compNames,repmat({']'},length(model.compNames),1)),strcat(repmat({'['},length(model.comps),1),model.comps,repmat({']'},length(model.comps),1)));
metstmp = mets;
for j = 1:length(CONValldata(:,1))
metstmp = strrep(metstmp,CONValldata(j,2),CONValldata(j,1));
end

[~,b] = ismember(metstmp,model.metNames);
mets(b~=0) = model.mets(b(b~=0));
%finding reactions which are not in the model
index=find(ismember(model.rxns,rxnID));

model = addReaction(model,rxnID,'metaboliteList',mets,'stoichCoeffList',S,'reversible',reversible,'geneRule',rule);

model.rxnNames(ismember(model.rxns,rxnID)) = rxnName;
end