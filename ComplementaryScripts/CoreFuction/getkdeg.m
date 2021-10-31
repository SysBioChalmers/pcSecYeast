function enzymedata = getkdeg(enzymedata)
enzymedata.kdeg = zeros(length(enzymedata.proteins),1);
[~,~,Turnovernumber]=xlsread('TableS1.xlsx','kdeg');
Turnovernumber = Turnovernumber(~strcmp(Turnovernumber(:,2),'xxx'),:);
[~,Idx] = ismember(enzymedata.proteins,Turnovernumber(:,2));
enzymedata.kdeg(Idx~=0,1) = cell2mat(Turnovernumber(Idx(Idx~=0),3))./(0.1+ cell2mat(Turnovernumber(Idx(Idx~=0),3)));
enzymedata.kdeg(enzymedata.kdeg == 0) = 0.042/(0.042+0.1); % 0.042/(0.042+0.1) median value from petri's paper
end