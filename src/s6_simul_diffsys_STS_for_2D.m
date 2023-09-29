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
workspaceDir = fullfile('..','workspace','\20230510\0-toolSelect');
% diaryFile = fullfile('..','workspace','\diary',['diary',datestr(now,'yyyymmddTHHMMSS')]);
% diary diaryFile;
% diary on;

%% load tool file
[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)'; ...
    '*,*','all files(*.*)'}, ...
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
    '*,*','all files(*.*)'}, ...
    'Select one cnc file', ...
    fullfile(['D:\OneDrive - sjtu.edu.cn\Research\Projects' ...
    '\202111-考虑刀具几何的路径规划\experiment\非球面加工\20230508-15\cnc'],'tooltheo.mat'), ...
    'MultiSelect','off');

cncName = fullfile(dirName,fileName);
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
cncData = zeros(3,0);
tmpLine = fgetl(cncFid);
while ~feof(cncFid)
    tmpLine = fgetl(cncFid);
    % if the line begins with %d%d or -%d, then break
    if strcmp(tmpLine,'(linking block)') || strcmp(tmpLine,'(LINKING BLOCK)')
        break;
    end
    cncData = [cncData,sscanf(tmpLine,'C%f X%f Z%f')];
end
% ncData = fscanf(fid,'X%fZ%f');
fclose(cncFid);

% diffsys - nanocam convertion
cncData(3,:) = cncData(3,:) - cncData(3,end);
cncData(2,:) = -1.*cncData(2,:);
cncData(1,:) = wrapTo360(-1.*cncData(1,:));
cncData(1,find(abs(cncData(1,:) - 360) < 1e-3)) = 0;
cncData(2:3,:) = 1000*cncData(2:3,:);

%% surface import
c = 0.69/1000/(1000^(aimUnit - presUnit));
syms x y;
surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
surfFunc = matlabFunction(surfSym);
surfFx = diff(surfFunc,x);
surfFy = diff(surfFunc,y);
surfNormFunc = matlabFunction([surfFx;surfFy;-1],'Vars',{'x','y'});
rMax = max(abs(cncData(2,:)));
% sampling density
spar = 501;
conR = linspace(0,rMax,spar); % concentric radius vector
conTheta = linspace(0,2*pi,spar);
[rMesh,thetaMesh] = meshgrid(conR,conTheta);
surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));

% other parameters
cutDirection = 'Edge to Center'; % 'Center to Edge'
startDirection = 'X Minus'; % 'X Minus'
angularIncrement = 'Constant Arc'; % 'Constant Angle'
aimRes = 1; % um
rStep = toolData.radius/2; % rStep can also be determined by the axial differentiate of the surface

% related parameters
isUIncrease = toolData.toolBform.coefs(end,1) - toolData.toolBform.coefs(1,1);
switch startDirection
    case 'X Plus' % plus in this program, but minus in moore
        rStep = -1*rStep;
    case 'X Minus' % minus in this program, but plus in moore
        rStep = 1*rStep;
end
if isUIncrease*rStep < 0
    % ([1,0] & X minus) or ([0,1] & X plus)
    uDirection = 'U Minus';
else
    % ([1,0] & X plus) or ([0,1] & X minus)
    uDirection = 'U Plus';
end

if strcmp(startDirection,'X Plus') % 'X Minus'
    conThetaBound = [0,2*pi];
    uLimOrder2 = 'first';
    uLimOrder3 = 'last';
else
    conThetaBound = [0,-2*pi];
    uLimOrder2 = 'last';
    uLimOrder3 = 'first';
end

%% 2-axes convertion to cartesian
spiralPtNum = size(cncData,2);
spiralPath = zeros(3,spiralPtNum);
spiralAngle = zeros(1,spiralPtNum);
spiralQuat = zeros(spiralPtNum,4);
spiralNorm = zeros(3,spiralPtNum);
spiralCut = zeros(3,spiralPtNum);
spiralPath(:,1) = [cncData(2,1);0;cncData(3,1)];
spiralQuat(1,:) = [1,0,0,0];    
for ii = 2:spiralPtNum
    spiralAngle(ii) = cncData(1,ii);
    spiralQuat(ii,:) = rotm2quat(rotz(spiralAngle(ii)));
    spiralNorm(:,ii) = quat2rotm(spiralQuat(ii,:))*toolData.toolEdgeNorm;
    spiralCut(:,ii) = quat2rotm(spiralQuat(ii,:))*toolData.cutDirect;
    spiralPath(1,ii) = cncData(2,ii)*cosd(spiralAngle(ii));
    spiralPath(2,ii) = cncData(2,ii)*sind(spiralAngle(ii));
    spiralPath(3,ii) = cncData(3,ii);
end

% shift
toolSp = toolData.toolBform;
toolSp.coefs = quat2rotm(spiralQuat(1,:))*toolSp.coefs + spiralPath(:,1);
toolSpPt = fnval(toolSp,0:0.0001:1);
spiralPath(3,:) = spiralPath(3,:) - min(toolSpPt(3,:));

figure;
% surf(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
plot3(spiralPath(1,:),spiralPath(2,:),spiralPath(3,:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');
% for ii = 1:spiralPtNum
%     toolSp1 = toolData.toolBform;
%     toolSp1.coefs = toolSp1.coefs + spiralPath(:,ii);
%     toolPt1 = fnval(toolSp1,0:0.01:1);
%     plot3(toolPt1(1,:),toolPt1(2,:),toolPt1(3,:),'Color',[0,0.4450,0.7410]);
% end

%% spiral tool path simulation and residual height calculation
spiralRes = 5*aimRes*ones(2,spiralPtNum);
spiralPeakPt = zeros(10,spiralPtNum);
spiralInterPtIn = cell(1,spiralPtNum);
spiralInterPtOut = cell(1,spiralPtNum);
spiralULim = cell(1,spiralPtNum);
tSpiralRes0 = tic;

switch uDirection
    case 'U Plus'
        uLimIni = [0;1];
    case 'U Minus'
        uLimIni = [1;0];
end

spiralContactU = zeros(1,spiralPtNum);
surfPt = zeros(3,spiralPtNum);
parfor ii = 1:spiralPtNum
    surfNorm = surfNormFunc(spiralPath(1,ii),spiralPath(2,ii));
    surfNorm = transpose(quat2rotm(spiralQuat(ii,:)))*surfNorm;
    [spiralContactU(ii),surfPt(:,ii),~] = toolPtInv(toolData.toolBform,surfNorm,1e-3, ...
        "Type",'TangentPlane',"Radius",toolData.radius);
end
toolRadius = toolData.radius;
for ind1 = 1:spiralPtNum
    % find the outer side of ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) - conThetaBound(end),1,uLimOrder2);
    ind3 = find(spiralAngle < spiralAngle(ind1) - conThetaBound(end),1,uLimOrder3);
    if isempty(ind2) || isempty(ind3)
        spiralULim{ind1} = uLimIni;
        tmpRes1 = 5*aimRes;
        tmpPeak1 = zeros(5,1);
    else
        [tmpRes1,tmpPeak1,spiralInterPtIn{ind1},spiralULim{ind1}] = ...
            residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolData,toolRadius,spiralULim{ind1},aimRes,uLimIni,ind1,ind2,ind3);
    end

    % find the inner side of ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + conThetaBound(end),1,uLimOrder2);
    ind3 = find(spiralAngle < spiralAngle(ind1) + conThetaBound(end),1,uLimOrder3);
    if isempty(ind2) || isempty(ind3)
        tmpRes2 = 5*aimRes;
        tmpPeak2 = zeros(5,1);
    else
        [tmpRes2,tmpPeak2,spiralInterPtOut{ind1},spiralULim{ind1}] = ...
            residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolData,toolRadius,spiralULim{ind1},aimRes,uLimIni,ind1,ind2,ind3);
    end
    
    spiralRes(:,ind1) = [tmpRes1;tmpRes2];
    spiralPeakPt(:,ind1) = [tmpPeak1;tmpPeak2];
    disp(ind1);
end

% get rid of the redundant part of the uLim of the innermost circle
rdomain = rMax;
if true
    parfor ii = 1:spiralPtNum
        if spiralULim{ii}(end) == 0
            spiralULim{ii}(end) = innermostU;
        end
        if spiralULim{ii}(1) == 1
            tmpU = spiralULim{ii}(2):abs(curvePlotSpar):1;
            toolSp1 = toolSp;
            toolSp1.coefs = quat2rotm(spiralQuat(ii,:))*toolSp1.coefs + spiralPath(:,ii);
            toolSp1Pt = fnval(toolSp1,tmpU);
            tmpPt = vecnorm(toolSp1Pt(1:2,:),2,1);
            [~,tmpUInd] = find(tmpPt == min(tmpPt(tmpPt > rdomain)));
            spiralULim{ii}(1) = tmpU(tmpUInd);
        end
    end
end
% delete(waitBar);
tSpiralRes = toc(tSpiralRes0);
fprintf('The time spent in the residual height calculation for spiral toolpath process is %fs.\n',tSpiralRes);
% warningTone = load('handel');
% sound(warningTone.y,warningTone.Fs);

%% spiral tool path result
% plot the result
figure('Name','Spiral tool path result');
tPlot0 = tic;
plotSpar = 1;
plot3(spiralPath(1,1:plotSpar:end), ...
    spiralPath(2,1:plotSpar:end), ...
    spiralPath(3,1:plotSpar:end), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none');
axis equal;
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
toolCoefs = toolSp.coefs;
stepNum = abs(log10(abs(curvePlotSpar)));
waitBar = waitbar(0,'Drawing ...','Name','Concentric Results Drawing', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for ii = 1:size(spiralPath,2)
    % Check for clicked Cancel button
    if getappdata(waitBar,'canceling')
        break;
    end
    displayData = ii/accumPtNum(end); % Calculate percentage
    waitbar(displayData,waitBar,['Figure Plotting ... ', ...
        num2str(roundn(displayData*100,-2),'%.2f'),'%']); % Progress bar dynamic display
    toolSp1 = toolSp;
    toolSp1.coefs = quat2rotm(spiralQuat(ii,:))*toolCoefs + spiralPath(:,ii);
    for jj = 1:size(spiralULim{ii},2)
        uLimRound = round(spiralULim{ii},stepNum);
        Q = fnval(toolSp1,uLimRound(1,jj):curvePlotSpar:uLimRound(2,jj));
        plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
        drawnow;
    end
end
delete(waitBar);

% axis equal;
grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('tool center point','','Location','northeast'); % 'tool edge',
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);

% sprial tool path error
s6_visualize_spiral_multi;

msgfig = helpdlg({sprintf(['\\fontsize{%d}\\fontname{%s}', ...
    'Spiral tool path was generated successfully!'], ...
    textFontSize,textFontType)},'Success');
