addpath(genpath('/proj/nobackup/snic2021-22-16/cplex1210/cplex/matlab'))
addpath(genpath('/home/f/feiranl/tools'))
addpath(genpath('/proj/nobackup/snic2021-22-16/SecYeast/'))
%addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/home/f/feiranl/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)