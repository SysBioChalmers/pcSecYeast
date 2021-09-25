addpath(genpath('/proj/nobackup/snic2021-22-16/cplex1210/cplex/matlab'))
addpath(genpath('/home/f/feiranl/tools'))
addpath(genpath('/proj/nobackup/snic2021-22-16/pcSecYeast/'))
%addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/home/f/feiranl/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
% %addpath('/apps/Common/Core/Gurobi/8.0.0/matlab')
% workdir = pwd;
% cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
% initCobraToolbox
% savepath '~/pathdef.m'
% cd(workdir)