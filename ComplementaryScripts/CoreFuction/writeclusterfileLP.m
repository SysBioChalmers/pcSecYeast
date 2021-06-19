function writeclusterfileLP(allname,subname)


subfileName = [subname,'.sh'];
fptr = fopen(subfileName,'w');
fprintf(fptr,'#!/bin/bash\n');
%fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
fprintf(fptr,'#SBATCH -n 20\n');
fprintf(fptr,'#SBATCH -o out.txt\n');
fprintf(fptr,'#SBATCH --time 0-3:00:00\n');
fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
fprintf(fptr,'#SBATCH --mail-type=end\n');
fprintf(fptr,'module load MATLAB GCCcore/10.2.0 GMP CMake/3.18.4\n');

for i = 1:length(allname)
    fprintf(fptr,['/home/f/feiranl/tools/soplex-4.0.0/build/bin/soplex  -s0 -g5 -t1000 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 ',allname{i},' > ', [allname{i},'.out'],' &\n']);
    %fprintf(fptr,['/cephyr/users/feiranl/Hebbe/tools/build/bin/soplex  -s0 -g5 -t1000 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-3 --real:fpopttol=1e-3 ',allname{i},' > ', [allname{i},'.out'],' &\n']);

end
fprintf(fptr,'wait;\n');



