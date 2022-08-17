% clear;
clc;
dirName = 'tool_path\';
fileName = 'Lens array_210721_RT0.2_tilt0_rot0_SF_s6a12d0.66.txt';
pathName = fullfile(dirName,fileName);
% data = importdata(pathName,' ',52);
data = readtable(pathName);
