

close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = 'mm';
textFontSize = 14;
textFontType = 'Times New Roman';
profile on
parObj = gcp;

%% to generate the function of residual height 
% related variables: spindle direction, cutting width and arc step.





%% 
% delete(parObj);
profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));