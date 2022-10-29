% to deal with the tool tip measurement data
% and to get the 3D point cloud of the tool tip

if true
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20221020-tooltip';
    fitMethod = 'Levenberg-Marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
end

%% 3D curve results
debug = 2;
switch debug
    case 1
        pathName = fullfile(workspaceDir,"tooltip result/20221019-strategy-2+40-5.csv");
    case 2
        optsInput.Resize = 'on';
        optsInput.WindowStyle = 'normal';
        optsInput.Interpreter = 'tex';
        toolInput = inputdlg({'Tool fitting method','Tool parameterization method', ...
            'Unit','Font type in figure','Font size in figure'}, ...
            'Tool Processing Input', ...
            [1 50; 1 50; 1 50; 1 50; 1 50], ...
            {'Levenberg-Marquardt','centripetal','\mu m','Times New Roman','14'},optsInput);
        fitMethod = toolInput{1};
        paramMethod = toolInput{2};
        unit = toolInput{3};
        textFontType = toolInput{4};
        textFontSize = str2double(toolInput{5});
        % tool measurement file loading
        [fileName,dirName] = uigetfile({ ...
            '*.csv','Comma-Separated Values-files(*.csv)'; ...
            '*.mat','MAT-files(*.mat)'; ...
            '*.txt','text-files(*.txt)'; ...
            '*.*','all files(*.*)'...
            }, ...
            'Select one tool tip measurement data', ...
            workspaceDir, ...
            'MultiSelect','off');
        pathName = fullfile(dirName,fileName);
end

% get rid of the header of the csv file
numHeader = 0;
tooltipFile = fopen(pathName);
while ~feof(tooltipFile)
    tmpLine = fgets(tooltipFile);
    if ~isnan(str2double(tmpLine(1:2)))
        break;
    end
    numHeader = numHeader + 1;
end
fclose(tooltipFile);
oriPts = importdata(pathName,',',numHeader);
oriPts = oriPts.data;
figure('Name','Original tool data');
plot(oriPts(:,1),oriPts(:,2),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);

% outliners removed
% rmoutliers(oriPts,2,"median");
% 如果多边形的直线段和圆弧段重合部分的数据差别比较大，说明目前的测量模式有问题。


%% ransac line fitting
% sampleSz = 2; % number of points to sample per trial
% maxDist = 2; % max allowable distance for inliers
% 
% fitLineFcn = @(pts) polyfit(pts(:,1),pts(:,2),1); % fit function using polyfit
% evalLineFcn = ...   % distance evaluation function
%   @(mdl, pts) sum((pts(:, 2) - polyval(mdl, pts(:,1))).^2,2);
% 
% [modelRANSAC,inlierIdx] = ransac(oriPts,fitLineFcn,evalLineFcn, ...
%   sampleSz,maxDist);
% 
% inlierPts = oriPts(inlierIdx,:);
% x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
% y = modelRANSAC(1)*x + modelRANSAC(2);
% plot(x, y, 'm-');

%% tool tip arc fitting and extraction

% ransac
sampleSz = 3; % number of points to sample per trial
maxDist = 0.05; % max allowable distance for inliers

fitLineFcn = @(pts) circleFit2D(pts');  % fit function
evalLineFcn = ...   % distance evaluation function
  @(mdl, pts) abs((pts(:,1) - x0).^2 + (pts(:,2) - y0).^2 - r^2);

[modelRANSAC,inlierIdx] = ransac(oriPts,fitLineFcn,evalLineFcn, ...
  sampleSz,maxDist);

















%%

% rmpath(genpath('.')