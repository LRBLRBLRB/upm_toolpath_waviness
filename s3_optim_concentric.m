close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = 'mm';
textFontSize = 14;
textFontType = 'Times New Roman';
% profile on
parObj = gcp("nocreate");

% [fileName,dirName] = uigetfile({ ...
%     '*.mat','MAT-files(*.mat)'; ...
%     '*,*','all files(*.*)'}, ...
%     'Select one tool edge data file', ...
%     'output_data\tool\tooltheo.mat', ...
%     'MultiSelect','off');
% toolName = fullfile(dirName,fileName);
toolName = 'output_data\tool\toolTheo_3D.mat';
toolData = load(toolName);

%% to generate the function of residual height 
% related variables: spindle direction, cutting width and arc step.

tool2D = toolData.toolBform;
tool2D.coefs(1,:) = [];
tool2D.dim = 2;
r = toolData.radius;

travNum = 150;
c1 = [0;0];
c2 = [linspace(0,150,travNum);zeros(1,travNum)];
vec1 = [0;-1];
vec2 = [zeros(1,travNum);ones(1,travNum)];
R1 = rotz(pi);
s1 = tool2D;
s1.coefs = R1(1:2,1:2)*s1.coefs;
R2 = rotz(pi);
R2 = R2(1:2,1:2);
res = zeros(1,travNum);
for ii = 1:travNum
    s2 = tool2D;
    s2.coefs = R2*s2.coefs + c2(:,ii);
    [res(ii),~] = residual2D_numeric(s1,s2,1e-3, ...
        c1,c2(:,ii),vec1,vec2(:,ii),r,'DSearchn');
end

figure('Name','Residual height function');
plot(c2(1,:),res);


%% 
% delete(parObj);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));