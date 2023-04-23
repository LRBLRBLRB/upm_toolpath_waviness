% The file is aimed to deal with the data measured by the Zygo white light interferometer, 
% and regenerate the measured surface

close all;
clear;
% clc;
addpath(genpath('funcs'));
if ~(exist('workspaceDir','var'))
    workspaceDir = '..\workspace\20230417';
    surfUnit = '\mum';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    unitList = {'m','mm','\mum','nm'};
    presUnit = find(strcmp(unitList,surfUnit),1);
    aimUnit = find(strcmp(unitList,unit),1);
end

% surface data import
if ~(exist('surfFunc','var'))
    % A = tand(20)/(2*2000);
    A = 0.091/1000/(1000^(aimUnit - presUnit));
    C = 0;
    syms x y;
    surfSym = A.*(x.^2 + y.^2)./2 + C;
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-1000,1000;-1000,1000];
    surfDomain = 1.05*surfDomain;
    rMax = max(surfDomain(1,2),surfDomain(2,2));
    % sampling density
    spar = 501;
    conR = linspace(0,rMax,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
end

%% 3D surface data calculation

% load the 3D original data
[surfFileName,surfDirName,surfFileType] = uigetfile({ ...
    '*.datx','data-files(*.datx)'; ...
    '*.xyz','xyz-files(*.xyz)'; ...
    '*.txt','text files(*.txt)'; ...
    '*.mat','MATLAB data files(*.mat)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one zygo data file', ...
    fullfile(workspaceDir), ...
    'MultiSelect','off');

surfFilePath = fullfile(surfDirName,surfFileName);

% extract the original data
switch surfFileType
    case 1
        % case for the .datx data
        surfData
    case 2
        % case for the .xyz data
        fid = fopen(surfFilePath);

        surfData
    case 3
        % case for the mat data exported from datx files
        surfData = load(surfFilePath);
    case 4
end


% plot the original result
figure('Name','1 Zygo original result');





%% 2D line data calculation
[lineFileName,lineDirName] = uigetfile({ ...
    '*.csv','comma-seperated-value files(*.dcsv)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one tool path data file', ...
    fullfile(workspaceDir), ...
    'MultiSelect','off');

lineFilePath = fullfile(lineDirName,lineFileName);
lineData = readmatrix(lineFilePath,"NumHeaderLines",1,'Range','A:B');

figure('Name','Zygo original 2D result');
plot(lineData(:,1),lineData(:,2));
hold on;
xlabel('r(\mum)');
ylabel('z(\mum)');

% 需要首先确定测量值的原点在哪里然后再减去

lineRange(1) = min(lineData(:,1));
lineRange(2) = max(lineData(:,2));
lineTheo = lineData(:,1);
lineTheo(:,2) = surfFunc(lineTheo(:,1),0);
plot(lineTheo(:,1),lineTheo(:,2));



%%
% rmpath(genpath('funcs'));