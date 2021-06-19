%% setMedia
function model = setMedia(model,type)

% type  = 1: minimal media (Delft media) (default)
%       = 2: yeast nitrogen base without amino acids
%       = 3: YNB+CSM-Ura (DOI: 10.1039/B803529F)
%       = 4: SD_2_SCAA
%       = 5: yeast nitrogen base with all amino acids

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange
% 1: minimal media (Delft media)                
desiredExchanges_1 = {'r_1654'; ... % ammonium exchange
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
                    'r_2020'};    % potassium exchange
% 2: yeast nitrogen base without amino acids
desiredExchanges_2 = {'r_1671'; ... % biotin exchange
                      'r_2067'; ... % thiamine exchange
                      'r_2028'; ... % pyridoxine exchange >high uptake<
                      'r_1967'; ... % nicotinate exchange
                      'r_1947'; ... % myo-inositol exchange
                      'r_1604'; ... % 4-aminobenzoate exchange
                      'r_1548'; ... % (R)-pantothenate exchange
                      'r_1792'; ... % folic acid exchange
                      'r_2038'};    % riboflavin exchange
% 3: CSM-Ura
CSM_Ura_Exchanges =  {'r_1893'; ... % L-histidine exchange
                      'r_1902'; ... % L-methionine exchange
                      'r_1912'; ... % L-tryptophan exchange
                      'r_1639'; ... % adenine exchange
                      'r_1879'; ... % L-arginine exchange
                      'r_1881'; ... % L-aspartate exchange
                      'r_1897'; ... % L-isoleucine exchange
                      'r_1899'; ... % L-leucine exchange
                      'r_1900'; ... % L-lysine exchange
                      'r_1903'; ... % L-phenylalanine exchange
                      'r_1911'; ... % L-threonine exchange
                      'r_1913'; ... % L-tyrosine exchange
                      'r_1914'};    % L-valine exchange

% 4: minimal + AA(14)
SD_2_SCAA             =    {'r_1879'; ...% L-arginine exchange
                            'r_1881';...%L-aspartate exchange
                            'r_1889';...%L-glutamate exchange
                            'r_1810';...%L-glycine exchange
                            'r_1893';...%L-histidine exchange
                            'r_1897';...%L-isoleucine exchange
                            'r_1899';...%L-leucine exchange
                            'r_1900';...%L-lysine exchange
                            'r_1902';...%L-methionine exchange
                            'r_1903';...%L-phenylalanine exchange
                            'r_1911';...%L-threonine exchange
                            'r_1912';...%L-tryptophan exchange
                            'r_1913';...%L-tyrosine exchange
                            'r_1914';...%L-valine exchange
                            'r_2090'};  %uracil exchange
% 5 YNB with aa
YNB                    =    {'r_1879'; ...% L-arginine exchange
                            'r_1881';...%L-aspartate exchange
                            'r_1889';...%L-glutamate exchange
                            'r_1810';...%L-glycine exchange
                            'r_1893';...%L-histidine exchange
                            'r_1897';...%L-isoleucine exchange
                            'r_1899';...%L-leucine exchange
                            'r_1900';...%L-lysine exchange
                            'r_1902';...%L-methionine exchange
                            'r_1903';...%L-phenylalanine exchange
                            'r_1911';...%L-threonine exchange
                            'r_1912';...%L-tryptophan exchange
                            'r_1913';...%L-tyrosine exchange
                            'r_1914';...%L-valine exchange
                            'r_1873';...%L-alanine exchange
                            'r_1880';...%L-asparagine exchange
                            'r_1883';...%L-cysteine exchange
                            'r_1891';...%L-glutamine exchange
                            'r_1904';...%L-proline exchange
                            'r_1906' }; %L-serine exchange
                           


uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges_1);
uptakeRxnIndexes_2     = findRxnIDs(model,desiredExchanges_2);

blockedRxnIndex      = findRxnIDs(model,blockedExchanges);
CSM_Ura_RxnIndex     = findRxnIDs(model,CSM_Ura_Exchanges);
SD_2_SCAA_RxnIndex = findRxnIDs(model,SD_2_SCAA);
YNB_RxnIndex = findRxnIDs(model,YNB);

if length(find(uptakeRxnIndexes~= 0)) ~= length(desiredExchanges_1)
    warning('Not all exchange reactions were found.');
end

model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0)) = -1000;

model.lb(blockedRxnIndex) = 0;
model.ub(blockedRxnIndex) = 0;

if type == 2
     model.lb(uptakeRxnIndexes_2) = -2; 
elseif type == 3
     model.lb(uptakeRxnIndexes_2) = -2; 
    model.lb(CSM_Ura_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 4
     model.lb(uptakeRxnIndexes_2) = -2; 
     model.lb(SD_2_SCAA_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
elseif type == 5
     model.lb(uptakeRxnIndexes_2) = -2; 
      model.lb(YNB_RxnIndex) = -0.08; % max coefficiency of AA in biomass times max mu
end

