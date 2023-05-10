% The file is aimed to deal with the data measured by the Mitutoyo formtracer, 
% and regenerate the measured surface.

% close all;
clear;
% clc;
addpath(genpath('funcs'));
if ~(exist('workspaceDir','var'))
    workspaceDir = fullfile('..','workspace','\20230417';
    unit = 'mm';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    unitList = {'m','mm','\mum','nm'};
    aimUnit = find(strcmp(unitList,unit),1);
end

%% data load
[lineFileName,lineDirName] = uigetfile({ ...
    '*.DAT','data files(*.DAT)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one measured result file', ...
    fullfile(workspaceDir), ...
    'MultiSelect','off');

lineFilePath = fullfile(lineDirName,lineFileName);
lineFileID = fopen(lineFilePath,'r');
lineUnit = fgetl(lineFileID);
if strcmp(lineUnit,'um')
    lineUnit = '\mum';
end
presUnit = find(strcmp(unitList,lineUnit),1);
lineNum = str2double(fgetl(lineFileID));
lineData0 = readmatrix(lineFilePath,"NumHeaderLines",2);

%% surface data import
if ~(exist('surfFunc','var'))
    % A = tand(20)/(2*2000);
    A = 0.091;
    C = 0;
    syms x;
    curveFunc = matlabFunction(A.*x.^2./2 + C);
    curveFx = diff(curveFunc,x);
    surfDomain = [-1,1;-1,1];
    surfDomain = 1.05*surfDomain;
end


if false
    %% 2D line data calculation

    % eliminate the outer-range segments
    lineDomain = surfDomain(1,:)/1.05;
    ind = lineData0(:,1) > lineDomain(1,1) & lineData0(:,1) < lineDomain(1,2);
    lineData = lineData0(ind,:);
    
    figure('Name','1 - Formtracer original 2D result');
    plot(lineData(:,1),lineData(:,2));
    hold on;
    xlabel('r(\mum)');
    ylabel('z(\mum)');
    
    % reference line
    line([min(lineData(:,1)),max(lineData(:,1))],[0,0],'LineStyle','--');
    line([0,0],[min(lineData(:,2)),max(lineData(:,2))],'LineStyle','--');
    
    % theoretical curve
    lineTheo = lineData(:,1);
    lineTheo(:,2) = curveFunc(lineTheo(:,1));
    plot(lineTheo(:,1),lineTheo(:,2));
else
    %% 2D substraction data calculation

    % eliminate the outer-range segments
    lineDomain = surfDomain(1,:)/1.05*0.9;
    ind = lineData0(:,1) > lineDomain(1,1) & lineData0(:,1) < lineDomain(1,2);
    lineData = lineData0(ind,:);
    
    figure('Name','1 - Formtracer original 2D result');
    plot(lineData(:,1),lineData(:,2),'Color',[0 0.4470 0.7410]);
    hold on;
    xlabel('r(\mum)');
    ylabel('z(\mum)');
    
    % reference line
    line([min(lineData(:,1)),max(lineData(:,1))],[0,0],'LineStyle','--','Color',[0.8500 0.3250 0.0980]);
    line([0,0],[min(lineData(:,2)),max(lineData(:,2))],'LineStyle','--','Color',[0.8500 0.3250 0.0980]);

    set(gca,'XLim',[min(lineData(:,1)),max(lineData(:,1))], ...
        'YLim',[min(lineData(:,2)),max(lineData(:,2))]);
end


%%
% rmpath(genpath('funcs'));