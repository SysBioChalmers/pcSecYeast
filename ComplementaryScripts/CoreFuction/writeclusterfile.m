%% Protein cost
num = 497;
k = 1:50:497;
for i = 1:length(k)
    subfileName = ['sub_PC',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');

    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:10
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(5*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(5*m-1) '\n']);
        end
        for m = 1:10
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateProteinCost($a',num2str(m),',$b',num2str(m),')" &\n']);
        fprintf(fptr,'sleep 30s\n');
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateProteinCost($a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 30s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end
%% read protein cost
num = 497;
k = 1:100:497;
for i = 1:length(k)
    subfileName = ['sub_Rd_proteincost',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:10
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(10*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(10*m-1) '\n']);
        end
        for m = 1:10
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadProteinCostResult($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(10*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(10*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadProteinCostResult($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end
%% glc con for fig 2
num = 20;
for i = 1:num
    subfileName = ['sub_glc',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1002\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');

    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateGlcCon(',num2str(i),')"']);
    fclose(fptr);
end

%% CPY % for fig 3c
TPrate = linspace(3.875E-7,1.9375e-04,6);
misP = [0.1:0.2:1,1];
order = [1,1;1,2;1,3;1,4;1,5;2,1;2,2;2,3;2,4;2,5;3,1;3,2;3,3;3,4;3,5;4,1;4,2;4,3;4,4;4,5;5,1;5,2;5,3;5,4;5,5;1,6;2,6;3,6;4,6;5,6;6,6;6,1;6,2;6,2;6,4;6,5];
for i = 1:length(order(:,1))
    subfileName = ['sub_CPY_SCE61_',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2022-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateCPY(',num2str(TPrate(order(i,1))),',',num2str(misP(order(i,2))),',5,0)"']);
    fclose(fptr);
end

%% CPY_Para for fig Sx
TPrate = linspace(3.875E-7,1.9375e-04,6);
TPrate = linspace(1.9375e-04,0.00124000000000000,6);
misP = [0.1:0.2:1,1];
order = [1,1;1,2;1,3;1,4;1,5;2,1;2,2;2,3;2,4;2,5;3,1;3,2;3,3;3,4;3,5;4,1;4,2;4,3;4,4;4,5;5,1;5,2;5,3;5,4;5,5;6,1;6,2;6,3;6,4;6,5];
for j = [0,4,5,7,8] %4 means total er constraint, 5 means total erm constraint, 6 means total retro-tranloc constraint, 7 means total ERAD constraint
    for i = 1:length(order(:,1))
    subfileName = ['sub_CPY_Para_',num2str(i),'_',num2str(j),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1002\n');
    %fprintf(fptr,'#SBATCH -A C3SE2022-1-16\n');
    fprintf(fptr,'#SBATCH -n 10\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    extraconstraint = i;
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateCPY_Para(',num2str(TPrate(order(i,1))),',',num2str(misP(order(i,2))),',5,0,',num2str(j),')"']);
  
    fclose(fptr);
    end
end

%% CPY_test % For Fig 3b
TPrate = [2e-6,5e-05,0]; % 0 is for reference
misP = [0,0.45,1];
fold = [5,10,20];
order = [1,1,1,0;1,2,1,0;1,3,1,0;1,3,2,1;1,3,3,1;2,1,1,0;2,2,1,0;2,3,1,0;2,3,2,1;2,3,3,1;3,1,1,0];
for i = 1:length(order(:,1))
    subfileName = ['sub_CPY_test',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2022-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateCPY(',num2str(TPrate(order(i,1))),',',num2str(misP(order(i,2))),',',num2str(fold(order(i,3))),',',num2str(order(i,4)),')"']);
    fclose(fptr);
end


%% FakeTP
%% Protein cost
num = 1112;
k = 1:100:1112;
for i = 1:length(k)
    subfileName = ['sub_FTP',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:5
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(20*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(20*m-1) '\n']);
        end
        for m = 1:5
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateFakeTP($a',num2str(m),',$b',num2str(m),')" &\n']);
        fprintf(fptr,'sleep 30s\n');
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateFakeTP($a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 30s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

%% Protein cost FakeTP
num = 1112;
k = 1:250:1112;
for i = 1:length(k)
    subfileName = ['sub_FTP_res',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-1004\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB/2021a GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:5
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(50*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(50*m-1) '\n']);
        end
        for m = 1:5
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadFakeTPres($a',num2str(m),',$b',num2str(m),')" &\n']);
        fprintf(fptr,'sleep 30s\n');
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadFakeTPres($a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 30s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end