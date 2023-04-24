% The file is aimed to deal with the data measured by the Zygo white light interferometer, 
% and regenerate the measured surface

% close all;
clear;
% clc;
addpath(genpath('funcs'));
if ~(exist('workspaceDir','var'))
    workspaceDir = '..\workspace\20230417';
    surfUnit = 'um';
    unit = 'um';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    unitList = {'m','mm','um','nm'};
    presUnit = find(strcmp(unitList,surfUnit),1);
    aimUnit = find(strcmp(unitList,unit),1);
end

%% surface data import
if ~(exist('surfFunc','var'))
    % A = tand(20)/(2*2000);
    A = 0.091/(1000^(aimUnit - presUnit));
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
    spar = 1000;
    conR = linspace(0,rMax,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
end

%% 3D surface data load
[surfFileName,surfDirName] = uigetfile({ ...
    '*.mat','MATLAB data files(*.mat)'; ...
    '*.datx','data-files(*.datx)'; ...
    '*.xyz','xyz-files(*.xyz)'; ...
    '*.txt','text files(*.txt)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one zygo data file', ...
    fullfile(workspaceDir), ...
    'MultiSelect','off');

surfFilePath = fullfile(surfDirName,surfFileName);
[~,~,surfFileType] = fileparts(surfFilePath);

% extract the original data
switch surfFileType
    case '.mat'
        % case for the mat data exported from datx files
        surfDataStruct = load(surfFilePath);
        surfData = zeros([size(surfDataStruct.mh),3]);
        [surfData(:,:,2),surfData(:,:,1)] = meshgrid(surfDataStruct.vx,surfDataStruct.vy);
        surfData(:,:,3) = surfDataStruct.mh;
        % unit convertion
        xyUnit0 = find(strcmp(unitList,surfDataStruct.XYunit),1); % default unit is milimeter
        zUnit0 = find(strcmp(unitList,surfDataStruct.Zunit),1); % default unit is meter
        aimUnit = find(strcmp(unitList,unit),1);
        surfData(:,:,1:2) = 1000^(aimUnit - xyUnit0)*surfData(:,:,1:2);
        surfData(:,:,3) = 1000^(aimUnit - zUnit0)*surfData(:,:,3);
    case '.datx'
        % case for the .datx data
        surfData
    case '.xyz'
        % case for the .xyz data
        surfData0 = loadXYZ(surfFilePath,1000,1000,8.70983e-07,'um','um');
        surfData(:,:,1) = surfData0.X;
        surfData(:,:,2) = surfData0.Y;
        surfData(:,:,3) = surfData0.Z;
%         surfData = zygo_xyz(surfFilePath);
    case '.txt'
        % case for the txt data files
        surfData = load(surfFilePath);
    otherwise
        return;
end

%         [surfFileName,surfDirName,surfFileType] = uigetfile({ ...
%             '*.datx','data-files(*.datx)'; ...
%             '*.xyz','xyz-files(*.xyz)'; ...
%             '*.txt','text files(*.txt)'; ...
%             '*.mat','MATLAB data files(*.mat)'; ...
%             '*.*','all files(*.*)'}, ...
%             'Select one zygo data file', ...
%             fullfile(workspaceDir), ...
%             'MultiSelect','off');
%         surfFilePath1 = fullfile(surfDirName,surfFileName);
%         surfDataStruct1 = load(surfFilePath1);
%         [surfData1(:,:,1),surfData1(:,:,2)] = meshgrid(surfDataStruct1.vx,surfDataStruct1.vy);
%         surfData1(:,:,3) = surfDataStruct1.mh;


% plot the original result
figure('Name','1 Zygo original result');
surf( ...
    surfData(:,:,1),surfData(:,:,2),surfData(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
% surf( ...
%     surfData1(:,:,1),surfData1(:,:,2),surfData1(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

%% point cloud registration
% using zyq's icp function
% surfDataPt = reshape(surfData,[],3);
% surfMeshPt = reshape(surfData,[],3);
% [surfRot,surfTran] = icp_yq(surfDataPt,surfMeshPt);
% surfDataPt = (surfRot*(surfDataPt.') + surfTran).';
% surfData = reshape(surfDataPt,size(surfData,1),size(surfData,2),3);

% using computer vision toolbox
% pcregisterndt
surfDataPC = pointCloud(surfData);
surfMeshPC = pointCloud(surfMesh);
% low-pass filter 
gridStep = 0.1; % down sampling
surfDataPC = pcdownsample(surfDataPC,'gridAverage',gridStep);
surfMeshPC = pcdownsample(surfMeshPC,'gridAverage',gridStep);
[surfRot,surfDataPC,rmse] = pcregistericp(surfDataPC,surfMeshPC);


% manual registration
% surfDataPt = reshape(surfData,[],3);
% surfMeshPt = reshape(surfData,[],3);
% surfRot = eye(3,3);
% tmpZ = 0 - min(surfDataPt(:,3));
% surfTran = [0;0;tmpZ];
% surfDataPt = (surfRot*(surfDataPt.') + surfTran).';
% surfData = reshape(surfDataPt,size(surfData,1),size(surfData,2),3);


figure('Name','2 measured surface and designed surface');
surf(surfData(:,:,1),surfData(:,:,2),surfData(:,:,3), ...
    'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','FaceAlpha',0.1);
hold on;
surf(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor',[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.1);



%% fit the theoretical surface
% z = A.*((x - x0).^2 + (y - y0).^2)./2 + C + z0;




%% surface error calculation
deltaZ = surfData(:,:,3) - surfFunc(surfData(:,:,1),surfData(:,:,2));
Sa = mean(abs(deltaZ));
Sz
Sq = max(deltaZ,[],'all') - min(deltaZ,[],'all');



%%
% rmpath(genpath('funcs'));