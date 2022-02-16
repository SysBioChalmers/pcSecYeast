function SimulateFakeTP(a,b)
% SimulateFakeTP
initcluster
% This is to simulate which factor would influence the productivity most
% we would create 400 fake secretory proteins with different modifications
% and amino acid composition

% generate fake proteins with 200 aa but differnet aa composition and
% differernt modification
% aa 20*5 + DSB 4 + NG 4 + OG 4 + Trans 2
% 
% [~,raw,~] = xlsread('aa_id.xlsx');
% aa_list = struct();
% aa_list.aa = raw(2:end,1);
% aa_list.subs = raw(2:end,3);
% aa_list.prod = raw(2:end,5);
% aa = aa_list.aa;
% aa_ratio = [0.2 0.4 0.6 0.8];
% modification_ratio = [1,2,3,4,5,6,7,8,9,10];
% 
% 
% % create fake proteins for different aa
% length_total = 400;
% seq_ave = cell2mat(join(aa,'')); % 20 aa;
% initial = '';
% k = 1;
% for i = 1:20
%     for j = 1:length(aa_ratio)
%         seq_fake_tmp = [initial,repmat(aa{i},1,length_total*aa_ratio(j))];
%         length_nolimit = length_total- length_total*aa_ratio(j);
%         seq_fake_tmp = [seq_fake_tmp,repmat(seq_ave,1,length_nolimit/20)];
%         seq_fake_tmp = [seq_ave,seq_fake_tmp]; % to create a SP
%         fakeProteinInfo(k,1) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,2) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,3) = {1}; % through ER
%         fakeProteinInfo(k,4) = {1}; % signal peptide
%         fakeProteinInfo(k,5) = {0}; % DSB
%         fakeProteinInfo(k,6) = {0}; %NG
%         fakeProteinInfo(k,7) = {0}; %OG
%         fakeProteinInfo(k,8) = {0}; %Trans
%         fakeProteinInfo(k,9) = {0}; %GPI
%         fakeProteinInfo(k,10) = {'e'}; % loc
%         fakeProteinInfo(k,11) = {seq_fake_tmp}; %seq
%         fakeProteinInfo(k,12) = {420}; % length
%         fakeProteinInfo(k,13) = {seq_ave}; % fake sp
%         fakeProteinInfo(k,14) = {20}; % fake sp
%         k = k+1;
%     end
% end
% % 
% % add modification proteins
% index = [5,6,7]; % DSB NG OG
% 
% for j = 1:length(index)
%     for i = 1:length(modification_ratio)
%     seq_fake_tmp = repmat(seq_ave,1,length_total/20 + 1); % +1 refers to SP
%     fakeProteinInfo(k,1) = {['fakeTP',num2str(k)]};
%     fakeProteinInfo(k,2) = {['fakeTP',num2str(k)]};
%     fakeProteinInfo(k,3:4) = {1};
%     fakeProteinInfo(k,5:9) = {0};
%     fakeProteinInfo(k,index(j)) = {modification_ratio(i)};
%     fakeProteinInfo(k,10) = {'e'}; % loc
%     fakeProteinInfo(k,11) = {seq_fake_tmp}; %seq
%     fakeProteinInfo(k,12) = {420}; % length
%     fakeProteinInfo(k,13) = {seq_ave}; % fake sp
%     fakeProteinInfo(k,14) = {20}; % fake sp
%     k = k+1;
%     end
% end
% 
% trans = [1,0];
% for i = 1:2
%         seq_fake_tmp = repmat(seq_ave,1,length_total/20 + 1); % +1 refers to SP
%          fakeProteinInfo(k,1) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,2) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,3:4) = {1};
%         fakeProteinInfo(k,5:9) = {0};
%         fakeProteinInfo(k,8) = {trans(i)};
%         fakeProteinInfo(k,10) = {'e'}; % loc
%         fakeProteinInfo(k,11) = {seq_fake_tmp}; %seq
%         fakeProteinInfo(k,12) = {420}; % length
%         fakeProteinInfo(k,13) = {seq_ave}; % fake sp
%         fakeProteinInfo(k,14) = {20}; % fake sp
%         k = k+1;
% end
% 
% % create hybrid sequence
% modification_ratio = [0,0,0,0,0,0,0,0,0,0,modification_ratio,modification_ratio,modification_ratio];
% initial = '';
% for i = 1:1000
% index_modi = randi([1 40],1,3);
% index_trans = randi([1 2],1,1);
% index_seq = floor(randfixedsum(20,1,400,1,400));
% index_seq(20) = 400-sum(index_seq(1:19));
% 
% seq(i,:) = index_seq';
% seq_fake_tmp = initial;
%         for j = 1:20
%                 seq_fake_tmp = [seq_fake_tmp,repmat(aa{j},1,index_seq(j))];
%         end
%         seq_fake_tmp = [seq_ave,seq_fake_tmp]; % to create a SP
%         
%         fakeProteinInfo(k,1) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,2) = {['fakeTP',num2str(k)]};
%         fakeProteinInfo(k,3:4) = {1};
%         fakeProteinInfo(k,5:7) = num2cell(modification_ratio(index_modi));
%         fakeProteinInfo(k,9) = {0};
%         fakeProteinInfo(k,8) = {trans(index_trans)};
%         fakeProteinInfo(k,10) = {'e'}; % loc
%         fakeProteinInfo(k,11) = {seq_fake_tmp}; %seq
%         fakeProteinInfo(k,12) = {length(seq_fake_tmp)}; % length
%         fakeProteinInfo(k,13) = {seq_ave}; % fake sp
%         fakeProteinInfo(k,14) = {20}; % fake sp
%         k = k+1;
% end
% % save('fakeProteinInfo.mat','fakeProteinInfo')
load('fakeProteinInfo.mat')
%% add fake protein to model and starts to simulate
% load model and param
load('enzymedata.mat')
load('enzymedataMachine.mat')
load('enzymedataSEC.mat')
load('enzymedataDummyER.mat');
load('pcSecYeast.mat')
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
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
f_mito = 0.1;
clear tot_protein f_modeled_protein;
factor_k = 1;
f_unmodelER = 0.046;
f_erm = 0.083;
mkdir('SimulateFakeTP')
cd SimulateFakeTP
allname = cell(0,1);
mulist = 0.1;
FakeTP = fakeProteinInfo(:,1);

for i = a:b
    [model_tmp,enzymedataTP] = addTargetProtein(model,FakeTP(i),false,fakeProteinInfo);
    [enzymedataTP] = SimulateRxnKcatCoef(model_tmp,enzymedataSEC,enzymedataTP);
    
    enzymedata_new = enzymedata;
    enzymedata_new.proteins = [enzymedata_new.proteins;enzymedataTP.proteins];
    enzymedata_new.proteinMWs = [enzymedata_new.proteinMWs;enzymedataTP.proteinMWs];
    enzymedata_new.proteinLength = [enzymedata_new.proteinLength;enzymedataTP.proteinLength];
    enzymedata_new.proteinExtraMW = [enzymedata_new.proteinExtraMW;enzymedataTP.proteinExtraMW];
    enzymedata_new.kdeg = [enzymedata_new.kdeg;enzymedataTP.kdeg];
    enzymedata_new.proteinPST = [enzymedata_new.proteinPST;enzymedataTP.proteinPST];
    enzymedata_new.rxns = [enzymedata_new.rxns;enzymedataTP.rxns];
    enzymedata_new.rxnscoef = [enzymedata_new.rxnscoef;enzymedataTP.rxnscoef];
    enzymedata_new = CombineEnzymedata(enzymedata_new,enzymedataSEC,enzymedataMachine,enzymedataDummyER);
    
    
    for k = 1:length(mulist)
        mu = mulist(k);
        model_tmp = changeRxnBounds(model_tmp,'r_2111',mu,'b');
            
            % Set optimization
            rxnID = [FakeTP{i},' exchange']; %maxmize FakeTP
            osenseStr = 'Maximize';
            
            name = FakeTP{i};
            fileName = writeLP(model_tmp,mu,f,f_unmodelER,osenseStr,rxnID,enzymedata_new,factor_k,name);
            allname = [allname;{fileName}];
    end
  
end

% write cluster file
writeclusterfileLP(allname,['sub_fakeTP_',num2str(a)])
display([num2str(a) '/' num2str(length(FakeTP))]);
cd ../
