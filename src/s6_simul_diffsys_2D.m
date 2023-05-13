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

% workspaceDir = fullfile('..','workspace','\20220925-contrast\nagayama_concentric';
% workspaceDir = fullfile('..','workspace','\20221020-tooltip\tooltip fitting result';
workspaceDir = fullfile('..','workspace','20230504 D906');
diaryFile = fullfile('..','workspace','diary',['diary',datestr(now,'yyyymmddTHHMMSS')]);
% diary diaryFile;
% diary on;

%% load tool file
[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one tool edge data file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');
toolName = fullfile(dirName,fileName);
% toolName = 'output_data\tool\toolTheo_3D.mat';
% tool data unit convertion
toolData = load(toolName);
unitList = {'m','mm','\mum','nm'};
presUnit = find(strcmp(unitList,toolData.unit),1);
aimUnit = find(strcmp(unitList,unit),1);
toolData.center = 1000^(aimUnit - presUnit)*toolData.center;
toolData.radius = 1000^(aimUnit - presUnit)*toolData.radius;
toolData.toolBform.coefs = 1000^(aimUnit - presUnit)*toolData.toolBform.coefs;
toolData.toolCpts = 1000^(aimUnit - presUnit)*toolData.toolCpts;
toolData.toolEdgePt = 1000^(aimUnit - presUnit)*toolData.toolEdgePt;
toolData.toolFit = 1000^(aimUnit - presUnit)*toolData.toolFit;

%% load nc file
[fileName,dirName] = uigetfile({ ...
    '*.nc;.pgm','CNC-files(*.nc,*pgm)'; ...
    '*.csv','Comma-separated-value files(*.csv)'; ...
    '*,*','all files(*.*)'}, ...
    'Select one cnc file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');

cncName = fullfile(dirName,fileName);
[~,~,cncFileType] = fileparts(cncName);

switch cncFileType
    case {'.nc'}
        cncFid = fopen(cncName,'r');
        % get rid of the header of the nc file
        numHeader = 1;
        while ~feof(cncFid)
            tmpLine = fgetl(cncFid);
            if strncmp(tmpLine,'#105',4)
                feedVel = sscanf(tmpLine,'#105=%d%.s');
                continue;
            end
            if strncmp(tmpLine,'#201',4)
                spindleVel = sscanf(tmpLine,'#201=%d%.s');
                continue;
            end
            % if the line begins with %d%d or -%d, then break
            if strcmp(tmpLine,'( CUTTING BLOCK )')
                break;
            end
            numHeader = numHeader + 1;
        end
        % load the X Z data
        cncData = zeros(2,0);
        while ~feof(cncFid)
            tmpLine = fgetl(cncFid);
            % if the line begins with %d%d or -%d, then break
            if strcmp(tmpLine,'(linking block)') || strcmp(tmpLine,'(LINKING BLOCK)')
                break;
            end
            cncData = [cncData,sscanf(tmpLine,'X%fZ%f')];
        end
        % ncData = fscanf(fid,'X%fZ%f');
        fclose(cncFid);
    case '.pgm'
    case '.csv'
end


%% surface import
A = tand(20)/(2*2000);
C = 0;
syms x y;
surfSym = A*(x.^2 + y.^2) + C;
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

%% parameter settings
spindleDirection = 'Clockwise';
if strcmp(spindleDirection,'Clockwise') % 'Counterclockwise'
    conThetaBound = [0,-2*pi];
else
    conThetaBound = [0,2*pi];
end

%% 2-axes convertion to cartesian
spiralPtNum = size(cncData,2);
spiralPath = zeros(3,spiralPtNum);
spiralAngle = zeros(1,spiralPtNum);
spiralQuat = zeros(spiralPtNum,4);
spiralPath(1,1) = cncData(1);
spiralQuat(:,1) = [1,0,0,0];    
for ii = 2:spiralPtNum
    spiralAngle(ii) = spiralAngle(ii - 1) + spindleVel/feedVel*(cncData(1,ii) - cncData(1,ii - 1));
    spiralQuat(:,ii) = rotm2quat(rotz(spiralAngle(ii)));
    spiralPath(1,ii) = cncData(1,ii)*cos(spiralAngle(ii));
    spiralPath(2,ii) = cncData(1,ii)*sin(spiralAngle(ii));
end
spiralPath = 1000*spiralPath;

figure;
% surf(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
plot(spiralPath(1,:),spiralPath(2,:),'-.','Color',[0.8500 0.3250 0.0980]);
% for ii = 1:spiralPtNum
%     toolSp1 = toolData.toolBform;
%     toolSp1.coefs = toolSp1.coefs + spiralPath(:,ii);
%     toolPt1 = fnval(toolSp1,0:0.01:1);
%     plot3(toolPt1(1,:),toolPt1(2,:),toolPt1(3,:),'Color',[0,0.4450,0.7410]);
% end

%% spiral tool path simulation and residual height calculation
spiralRes = nan(2,spiralPtNum);
spiralPeakPt = zeros(10,spiralPtNum);
spiralInterPtIn = cell(1,spiralPtNum);
spiralInterPtOut = cell(1,spiralPtNum);
spiralULim = cell(1,spiralPtNum);
tSpiralRes0 = tic;

parfor ind1 = 1:spiralPtNum
    % inner ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + conThetaBound(end),1,'first');
    ind3 = find(spiralAngle < spiralAngle(ind1) + conThetaBound(end),1,'last');
    if isempty(ind2) || isempty(ind3)
%         ind2 = find(spiralAngle >= spiralAngle(ind1) + pi,1,'first');
%         ind3 = find(spiralAngle < spiralAngle(ind1) + pi,1,'last');
%         [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
%             spiralPath,spiralNorm,spiralCut,spiralContactU, ...
%             toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
        spiralULim{ind1} = [0;1];
        tmpRes1 = nan;
        tmpPeak1 = zeros(5,1);
    else
        [tmpRes1,tmpPeak1,spiralInterPtIn{ind1},spiralULim{ind1}] = ...
            residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolData,toolRadius,spiralULim{ind1},ind1,ind2,ind3);
    end

    % outer ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + 2*pi,1,'first');
    ind3 = find(spiralAngle < spiralAngle(ind1) + 2*pi,1,'last');
    if isempty(ind2) || isempty(ind3)
        tmpRes2 = nan;
        tmpPeak2 = zeros(5,1);
    else

%         scatter3(spiralPath(1,ind1),spiralPath(2,ind1),spiralPath(3,ind1));
%         quiver3(spiralPath(1,ind1),spiralPath(2,ind1),spiralPath(3,ind1), ...
%             spiralCut(1,ind1),spiralCut(2,ind1),spiralCut(3,ind1));
%         toolSp1 = toolData.toolBform;
%         toolSp1.coefs = quat2rotm(spiralQuat(ind1,:))*toolSp1.coefs + spiralPath(:,ind1);
%         fnplt(toolSp1);
%         scatter3(spiralPath(1,ind2),spiralPath(2,ind2),spiralPath(3,ind2));
%         quiver3(spiralPath(1,ind2),spiralPath(2,ind2),spiralPath(3,ind2), ...
%             spiralCut(1,ind2),spiralCut(2,ind2),spiralCut(3,ind2),'AutoScale','on');
%         scatter3(spiralPath(1,ind3),spiralPath(2,ind3),spiralPath(3,ind3));
%         quiver3(spiralPath(1,ind3),spiralPath(2,ind3),spiralPath(3,ind3), ...
%             spiralCut(1,ind3),spiralCut(2,ind3),spiralCut(3,ind3));

        [tmpRes2,tmpPeak2,spiralInterPtOut{ind1},spiralULim{ind1}] = ...
            residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolData,toolRadius,spiralULim{ind1},ind1,ind2,ind3);
    end
    spiralRes(:,ind1) = [tmpRes1;tmpRes2];
    spiralPeakPt(:,ind1) = [tmpPeak1;tmpPeak2];

%     if spiralRes(:,ind1) > 5
%         1;
%     end

    % debug
    % plot3(toolPathPtRes(1,ii),toolPathPtRes(2,ii),toolPathPtRes(3,ii), ...
    %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
    % toolSp1 = toolSp;
    % R1 = axesRot([0;0;1],[1;0;0],toolNormDirectRes(:,ii),toolCutDirectRes(:,ii),'zx');
    % toolSp1.coefs = R1*toolSp.coefs + toolPathPtRes(:,ii);
    % Q = fnval(toolSp1,uLimTmp(1,ii):0.01:uLimTmp(2,ii));
    % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1);
end

tSpiralRes = toc(tSpiralRes0);
fprintf('The time spent in the residual height calculation for spiral toolpath process is %fs.\n',tSpiralRes);






























