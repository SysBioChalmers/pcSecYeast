function [sol_obj,sol_status,sol_full] = readSoplexResult(fileName_out,model)

sol_full = zeros(length(model.rxns),1);

fid = fopen(fileName_out);

res = {};
line = fgets(fid);
while ischar(line)
    res = [res;{line}];
    line = fgets(fid);
end

fclose(fid);

format longE;

if find(contains(res,'problem is solved [optimal]'),1) > 0
    sol_status = 'optimal';
    obj_line = cell2mat(res(find(contains(res,'Objective value'),1)));
	sol_obj = str2double(obj_line(strfind(obj_line,':')+2:end-1));
    a = find(contains(res,'Primal solution (name, value):'));
    z = find(contains(res,'All other variables are zero'));
    flux = split(res(a+1:z-1));
    rxn_idx = flux(:,1);
    rxn_v = flux(:,2);
    for i = 1:length(sol_full)
        id = strcat('X',num2str(i));
        if ~isempty(find(strcmp(rxn_idx,id),1))
            v = str2double(cell2mat(rxn_v(strcmp(rxn_idx,id))));
        else
            v = 0;
        end
        sol_full(i) = v;
    end
else
    sol_status = 'no solution';
    sol_obj = [];
    sol_full = [];
end

