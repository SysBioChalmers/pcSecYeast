%% Protein cost
num = 505;
k = 1:50:505;
for i = 1:length(k)
    subfileName = ['sub_proteincost',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB GMP CMake\n');
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
num = 505;
k = 1:100:505;
for i = 1:length(k)
    subfileName = ['sub_proteincost',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
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
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(20*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(20*m-1) '\n']);
        end
        for m = 1:10
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadProteinCostResult($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(20*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(20*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "ReadProteinCostResult($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end
%% glc con
num = 20;
for i = 1:num
    subfileName = ['sub_glc',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A snic2021-22-16\n');
    fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    %fprintf(fptr,'module load MATLAB GMP CMake\n');
    fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');

    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateGlcCon(',num2str(i),')"']);
    fclose(fptr);
end

%% CPY
TPrate = linspace(3.08E-7,1.9250e-04,6);
misP = [0.1:0.2:1,1];
order = [1,1;1,2;1,3;1,4;1,5;2,1;2,2;2,3;2,4;2,5;3,1;3,2;3,3;3,4;3,5;4,1;4,2;4,3;4,4;4,5;5,1;5,2;5,3;5,4;5,5;1,6;2,6;3,6;4,6;5,6;6,6;6,1;6,2;6,2;6,4;6,5];
for i = 1:length(order(:,1))
    subfileName = ['sub_CPY_SCE61_',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    %fprintf(fptr,'module load MATLAB GMP CMake\n');
    fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateCPY(',num2str(TPrate(order(i,1))),',',num2str(misP(order(i,2))),',5,1)"']);
    fclose(fptr);
end

%% CPY_test
TPrate = 2.85e-06;
misP = [0,0.5,1];
fold = [5,10,20];
order = [1,1,1,0;1,2,1,0;1,3,1,0;1,3,2,1;1,3,3,1];
for i = 1:length(order(:,1))
    subfileName = ['sub_CPY_test',num2str(i),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    %fprintf(fptr,'module load MATLAB GMP CMake\n');
    fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['matlab -nodesktop -singleCompThread -r "SimulateCPY(',num2str(TPrate(order(i,1))),',',num2str(misP(order(i,2))),',',num2str(fold(order(i,3))),',',num2str(order(i,4)),')"']);
    fclose(fptr);
end

%% FakeTP
%% Protein cost
num = 300;
k = 1:100:300;
for i = 1:length(k)
    subfileName = ['sub_FTP',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-16\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB GMP CMake\n');
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
num = 300;
k = 1:50:300;
for i = 1:length(k)
    subfileName = ['sub_FTP_res',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-16\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'ml MATLAB GMP CMake\n');
    %fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
            fprintf(fptr,['a=' num2str(50*(i-1)+1) '\n']);
            fprintf(fptr,['b=' num2str(50*i) '\n']);
            fprintf(fptr,'matlab -nodesktop -singleCompThread -r ReadFakeTPres(a,b)');
    else
          fprintf(fptr,['a=' num2str(50*(i-1)+1) '\n']);
            fprintf(fptr,['b=' num2str(num) '\n']);
            fprintf(fptr,'matlab -nodesktop -singleCompThread -r ReadFakeTPres(a,b)');
    end
    fclose(fptr);
end