function ReadFakeTPres(a,b)
initcluster
load('pcSecYeast.mat');
cd SimulateFakeTP/;
file = dir('*.out');
filename = {file.name};
maxTP = [];
for j = a:b
    [sol_obj,sol_status,~] = readSoplexResult(['Simulation_dilutionfakeTP',num2str(j),'.lp.out'],model);
    if strcmp(sol_status,'optimal')
        maxTP(j) = sol_obj;
    else
       maxTP(j) = 0;
    end
end

save(['res_maxTP',num2str(a),'.mat'],'maxTP');
