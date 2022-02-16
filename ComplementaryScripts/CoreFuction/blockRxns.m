%% blockRxns 
function model = blockRxns(model)
% Some reactions should be block to avoid weird flux distributions.

% block Gcy1, an alternative glycerol dissimilation pathway that is active 
% under microaerobic conditions (PMID: 22979944)
model = changeRxnBounds(model,'r_0487',0,'b');


% block newly added isozyme
model = changeRxnBounds(model,'r_0438_5',0,'b'); % ferrocytochrome-c:oxygen oxidoreductase


% block newly added isozyme for enolase as Eno1 and Eno2 show higher protein levels
model = changeRxnBounds(model,'r_0366_1_fwd',0,'b'); % enolase
model = changeRxnBounds(model,'r_0366_1_rvs',0,'b'); % enolase

model = changeRxnBounds(model,'r_0366_5_fwd',0,'b'); % enolase
model = changeRxnBounds(model,'r_0366_5_rvs',0,'b'); % enolase


model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK

model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase

model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase


% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]

model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)

model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)


% % model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
% model.ub(ismember(model.rxns,'r_0725_fwd_withoutcofactor')) = 0;
% model.lb(ismember(model.rxns,'r_0725_fwd_withoutcofactor')) = 0;
% % model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase
% model.ub(ismember(model.rxns,'r_0918_withoutcofactor')) = 0;
% model.lb(ismember(model.rxns,'r_0918_withoutcofactor')) = 0;

model = changeRxnBounds(model,'r_4235',0,'b'); % weird reaction from glc to g6p

model = changeRxnBounds(model,'r_4216_rvs',0,'b'); % block the reaction to produce FMN without ATP

model = changeRxnBounds(model,'r_4264_rvs',0,'b'); % succinate:NAD+ oxidoreductase

% The following two reactions account for the transport of pyruvate from
% [c] to [m] without enzyme cost, should be blocked.
model = changeRxnBounds(model,'r_1137',0,'b');
model = changeRxnBounds(model,'r_1138',0,'b');

% Block backward reactions for PLP synthesis.
plpbwrxns = {'r_4211_1_rvs' 'r_4211_2_rvs' 'r_4212_1_rvs' 'r_4212_2_rvs'};
model = changeRxnBounds(model,plpbwrxns,zeros(1,length(plpbwrxns)),'b');

% block the cycle of GSH and GSSG

gshrxns = {'r_0483_1_rvs' 'r_0483_2_rvs' 'r_0483_3_rvs' 'r_0483_4_rvs' 'r_0483_5_rvs' 'r_0483_6_rvs'};
model = changeRxnBounds(model,gshrxns,zeros(1,length(gshrxns)),'b');

model = changeRxnBounds(model,'r_0446_rvs',0,'b');% blocked in the recent version yeast8

