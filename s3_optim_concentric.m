close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

msgOpts.Default = 'Cancel and quit';
msgOpts.Interpreter = 'tex';
% msgOpts.modal = 'non-modal';
profile on
parObj = gcp;

%%
% load the data of the residual function
load();





















%%
% delete(parObj);
profile off
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));