% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

% 实际上，我打算把concentric和freeform的方案给合并。目前aspheric文件中的是旧的刀位点计算方案，freeform中是新的

isAPP = false;
if isAPP
    workspaceDir = app.workspaceDir;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;

    % machining paramters
    cutDirection = app.cutDirection;
    spindleDirection = app.spindleDirection;
    angularDiscrete = app.angularDiscrete;
    aimRes = app.aimRes;
    toolData = app.toolData;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = app.maxIter;
    arcLength = app.arcLength;
    maxAngPtDist = app.maxAngPtDist;
    angularLength = app.angularLength;

    surfFunc = app.surfFuncs;
    surfFx = app.surfFx;
    surfFy = app.surfFy;
    surfDomain = app.surfDomain;
    surfMesh = app.surfMesh;
    Geometry2DCell = app.Geometry2DCell;
    surfType = app.surfType;
    rMax = app.rMax;
    
    msgOpts.Default = 'Cancel and quit';
    msgOpts.Interpreter = 'tex';
    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);
else
    close all;
    clear;
    clc;
    addpath(genpath('funcs'));
    % global variables
    % global textFontSize textFontType;
    % workspaceDir = '../workspace/20220925-contrast/nagayama_concentric';
    workspaceDir = '..\workspace\20220925-contrast\nagayama_concentric';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    
    msgOpts.Default = 'Cancel and quit';
    msgOpts.Interpreter = 'tex';
    % msgOpts.modal = 'non-modal';
    % profile on
    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);
    
    %% concentric surface generation / import
    % tool data import
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select one tool edge data file', ...
        fullfile(workspaceDir,'tooltheo.mat'), ...
        'MultiSelect','off');
    toolName = fullfile(dirName,fileName);
    % toolName = 'output_data\tool\toolTheo_3D.mat';
    toolData = load(toolName);
%     toolData.radius = 1000*toolData.radius;
%     toolData.toolBform.coefs = 1000*toolData.toolBform.coefs;

    % surface data import
    default = false;
    if default
        [fileName,dirName] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','all files(*.*)'}, ...
            'Select the surface edge data file', ...
            '..\workspace\input_data\surface\ellipsoidAray.mat', ...
            'MultiSelect','off');
        surfName = fullfile(dirName,fileName);
        load(surfName);
        D = 0;
        R = 10/2;
        A = 3.5/2;
        B = 4/2;
        C = 5/2;
        % machining surface
        syms x y;
        surfSym = D - C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
        surfFunc = matlabFunction(surfSym);
        surfFx = diff(surfFunc,x);
        surfFy = diff(surfFunc,y);
        rMax = R/2;
    else % ellipsoid
        D = 0;
        R = 4*1000;
        A = 3/2;
        B = 4/2;
        C = 5/2;
        % machining surface
        syms x y;
        surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
        surfFunc = matlabFunction(surfSym);
        surfFx = diff(surfFunc,x);
        surfFy = diff(surfFunc,y);
        surfDomain = [-A*R/4,A*R/4;-B*R/4,B*R/4];
        rMax = max(surfDomain(1,2),surfDomain(2,2));
        % rMax = 5;
        % sampling density
        spar = 501;
        % surfCenter = [0;0;sqrt(C^2*(R.^2-R/4.^2))]; % concentric circle center
        conR = linspace(0,rMax,spar); % concentric radius vector
        conTheta0 = linspace(0,2*pi,spar);
        [rMesh,thetaMesh] = meshgrid(conR,conTheta0);
        % surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
        % surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
        surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
        surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
        surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
        % save('input_data/surface/ellipsoidAray.mat', ...
        %    "surfMesh","surfNorm","surfCenter");
    end
    
    surfType = '3D';
    Geometry2DCell = {'Rotating Paraboloid','Aspheric'};

    % machining paramters
    cutDirection = 'Edge to Center'; % 'Center to Edge'
    spindleDirection = 'Counterclockwise'; % 'Clockwise'
    angularDiscrete = 'Constant Arc'; % 'Constant Angle'
    aimRes = 1;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = 10;
    arcLength = 30;
    maxAngPtDist = 6*pi/180;
    angularLength = 6*pi/180;
end

% plot the initial surface
[surfNormIni(:,:,1),surfNormIni(:,:,2),surfNormIni(:,:,3)] = surfnorm( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));

fig1 = figure('Name','original xyz scatters of the surface (sparsely)');
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
%     quiver3(surfMesh(1:10:end,1:10:end,1), ...
%         surfMesh(1:10:end,1:10:end,2), ...
%         surfMesh(1:10:end,1:10:end,3), ...
%         surfNormIni(1:10:end,1:10:end,1), ...
%         surfNormIni(1:10:end,1:10:end,2), ...
%         surfNormIni(1:10:end,1:10:end,3), ...
%         'AutoScale','on','Color',[0.85,0.33,0.10], ...
%         'DisplayName','Normal Vectors');
legend('Original Points','Orthogonal direction','Location','northeast');
% axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

% interaction
msgfig = msgbox('Surface was generated successfully!','Exit','help','non-modal');
uiwait(msgfig);
msgfig = questdlg({'Surface was generated successfully!', ...
    'Ready to continue?'}, ...
    'Surface Generation','OK & continue','Cancel & quit','OK & continue');
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

% pre-processing
surfNormSym = [surfFx;surfFy;-1];
surfNormSym = surfNormSym./norm(surfNormSym);
syms x y;
surfNormFunc = matlabFunction(surfNormSym,'Vars',{x,y});

% switch spindleDirection
%     case 'Clockwise'
%         conThetaBound = [2*pi,0];
%     case 'Counterclockwise'
%         conThetaBound = [0,2*pi];
% end
% 
% switch cutDirection
%     case 'Edge to Center'
%         interpR = fnval(Fr,accumPtNum(end):-1:accumPtNum(1) + 1);
%         interpR = [interpR,zeros(1,accumPtNum(1))];
%         accumPtNumlength = [0,accumPtNum];
%     case 'Center to Edge'
%     otherwise
%         errordlg('Invalid cut-direction value!','Parameter Error');
% end

% the spiral toolpath of clockwise & center-edge and that of
% counterclockwise & edge-center remains the same except for the direction
if (strcmp(spindleDirection,'Clockwise') && strcmp(cutDirection,'Center to Edge')) || ...
    (strcmp(spindleDirection,'Counterclockwise') && strcmp(cutDirection,'Edge to Center'))
    conThetaBound = [0,2*pi];
else
    conThetaBound = [2*pi,0];
end

% slCharacterEncoding('UTF-8');

%% load the data of the residual function, and get the appropriate cutting width
% load(fullfile(workspaceDir,"widthRes.mat")); % widthRes
% wid = interp1(widthRes(2,:),widthRes(1,:),aimRes);



%% Tool path adjusting for the first loop
t0 = tic;
tRes0 = tic;
t1 = tic;

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius

switch angularDiscrete
    case 'Constant Arc'
        conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
            ceil(2*pi/maxAngPtDist) + 1);
        conTheta0(end) = [];
    case 'Constant Angle'
        conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
            ceil(2*pi/angularLength) + 1);
        conTheta0(end) = [];
    otherwise
        error('Invalid angular discretization type.');
end

clear global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision toolPathAngle ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision toolPathAngle ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

% initialize the cutting pts: the rotation center is seen as the first loop
loopPtNum = [length(conTheta0)]; % the number of points in each loop
accumPtNum = [length(conTheta0)]; % the No. of the last points in each loop
toolREach = 0; % the radius R of the current loop
toolRAccum = zeros(1,loopPtNum(end));
toolNAccum = zeros(1,loopPtNum(end)); % the number of loops of each point
toolPathAngle(1:accumPtNum(1)) = conTheta0; % the tool path angle of each point
surfPt(1,:) = zeros(1,loopPtNum(end)); % x coordinates of surface cutting points
surfPt(2,:) = zeros(1,loopPtNum(end)); % y coordinates of surface cutting points
surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:)); % z coordinates of surface cutting points
surfNorm = surfNormFunc(surfPt(1,:),surfPt(2,:)); % normolization of the normal vector
surfDirect = cutdirection(surfPt,[0;0;0],conTheta0,'method','vertical'); % cutting direction of each points
% surfDirect = [[0;1;0],cutdirection(surfPt(:,2:end),[0;0;0])]; % cutting direction of each points
% 和[0;0;1]垂直、且与surfDirect(:,end)夹角最小的向量作为surDirect(:,1)，可视作surfDirect(:,end)到[0;0;1]平面的投影
% surfDirect(:,1) = vecOnPlane(surfDirect(:,end),surfPt(:,1),[0;0;1]);

% tool path initialization
toolQuat = zeros(accumPtNum(end),4); % the quaternions for each points
toolVec = zeros(3,accumPtNum(end)); % the translation vectors for each points
toolContactU = zeros(1,accumPtNum(end)); % the parameter u of the contact point on the tool tip
isCollision = false(1,accumPtNum(end));
toolPathPt = zeros(3,accumPtNum(end));
toolCutDirect = zeros(3,accumPtNum(end));
toolNormDirect = zeros(3,accumPtNum(end));
uLim = [zeros(1,accumPtNum(end));ones(1,accumPtNum(end))]; % the interval of each toolpath
peakPt = zeros(6,accumPtNum(end));
res = 5*aimRes*ones(2,accumPtNum(end)); % the residual height, initialized with 5 times the standard aimRes

% figure;
% calculate the 1st loop of the tool path
for ii = 1:accumPtNum(end)
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = tooltippos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),[0;0;-1],surfDirect(:,ii));
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
    % plot3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii)); hold on;
    % quiver3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
    %     toolCutDirect(1,ii),toolCutDirect(2,ii),toolCutDirect(3,ii), ...
    %     'AutoScale','on');
    
%     toolSp1 = toolSp;
%     toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolSp.coefs + toolVec(:,ii);
%     Q = fnval(toolSp1,uLim(1,ii):0.01:uLim(2,ii));
%     plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
%     xlabel('x'); ylabel('y'); axis equal; hold on;
end

fprintf('No.1\tElapsed time is %f seconds.\n-----\n',toc(t1));

%% Tool path adjusting for the rest
r = rStep; % 可以通过widthRes确定迭代初值
delta = rStep;

% specially for the 2nd loop 
% ---------------- (I have no idea how to put it together with the other
% loops!!! My god. Damn it!!!) ----------------
r = iterfunc_surfpath01(r,toolData,delta,surfFunc,surfNormFunc, ...
    angularDiscrete,conThetaBound,maxAngPtDist,arcLength,angularLength, ...
    aimRes,rMax,conTheta0);

iterfunc_surfpath(r,toolData,delta,surfFunc,surfNormFunc, ...
    angularDiscrete,conThetaBound,maxAngPtDist,arcLength,angularLength, ...
    aimRes,rMax,tRes0);

% debug
% figure('Name','concentric tool path debug');
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',1,'LineStyle','none'); hold on;
% colormap('summer');
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel('x'); ylabel('y'); zlabel('z');
% view(0,90);



%% plot the result above
figure('Name','tool path optimization');
plotSpar = 1;
tPlot0 = tic;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none'); hold on;
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
ptNum = sum(loopPtNum);
toolCoefs = toolSp.coefs;
% figure;
if size(surfDomain,1) == 1
    for jj = 1:ptNum
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolVec(:,jj);
        Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
        % isQDomain = (Q(1,:).^2 + Q(2,:).^2 - surfDomain(1,2)^2) >= 0;
        % Q(:,isQDomain) = [];
        plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1); hold on;
    end
else
    for jj = 1:ptNum
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolVec(:,jj);
        Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
%         isQXDomain = abs(Q(1,:)) - surfDomain(1,2) >= 0;
%         Q(:,isQXDomain) = [];
%         isQYDomain = abs(Q(2,:)) - surfDomain(2,2) >= 0;
%         Q(:,isQYDomain) = [];
%         isQDomain = Q(1,:).^2/A^2 + Q(2,:).^2/B^2 - ;
        plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1); hold on;
    end
end

% axis equal;
grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
axis equal;
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
legend('designed surface','tool center point','tool edge','Location','northeast');
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);

% save ..\workspace\20220925-contrast\nagayama_concentric\loopR.mat toolREach

s4_visualize_concentric;
% 这里还有问题：中间一圈的残高比外圈的小很多，而且中间的一圈残高的计算是有问题的。

msgfig = questdlg({'Concentric tool path was generated successfully!', ...
    'Ready to continue?'}, ...
    'Concentric tool path Generation','OK & continue','Cancel & quit','OK & continue');
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% Feed rate smoothing
% to smooth the loopR to get the real tool path

% cubic spline approximation
% the function between the numeric label of tool path and surf radius R
Fr = csape(accumPtNum,toolREach,[1,1]);

toolNoTheta = linspace(2*pi*1,2*pi*length(accumPtNum),length(accumPtNum));
rTheta = csape(toolNoTheta,toolREach,[1,1]);

figure('Name','Feed Rate Smoothing');
tiledlayout(2,1);
nexttile;
scatter(accumPtNum,toolREach);
hold on;
fnplt(Fr,'r',[accumPtNum(1) + 1,accumPtNum(end)]);
plot(1:accumPtNum(end),toolRAccum);
% line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
%     'Color',[0.929,0.694,0.1250]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel('Loop Accumulating Point Number');
ylabel(['Radius of the Loop (',unit,')']);
legend('No.-R scatters','csape result','Concentric result');
nexttile;
scatter(toolNoTheta,toolREach);
hold on;
fnplt(rTheta,'r',[toolPathAngle(accumPtNum(1) + 1),toolPathAngle(end)+2*pi*toolNAccum(end)]);
plot(toolPathAngle + 2*pi*toolNAccum,toolRAccum);
% line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
%     'Color',[0.929,0.694,0.1250]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel('Accumulating Toolpath Angle');
ylabel(['Radius of the Loop (',unit,')']);
legend('\theta-R scatters','csape result','Concentric result');
% save the feed rate curve
[smoothFileName,smoothDirName,smoothFileType] = uiputfile({ ...
    '*.mat','MAT-file(*.mat)'; ...
    '*.txt','text-file(.txt)';...
    '*.*','all file(*.*)';...
    }, ...
    'Select the directory and filename to save the surface concentric tool path', ...
    fullfile(workspaceDir,['toolPath',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
smoothName = fullfile(smoothDirName,smoothFileName);
switch smoothFileType
    case 0
        msgfig = msgbox("No approximation saved","Warning","warn","non-modal");
        uiwait(msgfig);
    case 1
        Comments = cell2mat(inputdlg( ...
            'Enter Comment of the feed rate smoothing processing:', ...
            'Saving Comments', ...
            [5 60], ...
            string(datestr(now))));
        save(smoothName,"Comments","Fr","rTheta");
end


%% spiral tool path generation with the smoothing result
% initialize the spiral path
spiralPath0 = zeros(3,accumPtNum(end)); % the spiral tool path & orientation
spiralNorm0 = zeros(3,accumPtNum(end));
spiralCut0 = zeros(3,accumPtNum(end));
spiralQuat0 = zeros(accumPtNum(end),4);
spiralContactU0 = zeros(1,accumPtNum(end));

interpR = fnval(Fr,1:accumPtNum(end));
interpR(1:accumPtNum(1)) = 0;
accumPtNumlength = [0,accumPtNum];

numLoop = length(accumPtNum);

% video capture & debug
% figure('Name','concentric-to-spiral debug & video');
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',1,'LineStyle','none'); hold on;
% colormap('summer');
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel('x'); ylabel('y'); view(0,90);
% plot3(spiralPath(1,1:accumPtNum(1)),spiralPath(2,1:accumPtNum(1)), ...
%     spiralPath(3,1:accumPtNum(1)),'.','MarkerSize',6,'Color',[0,0.4470,0.7410]);

tSpiral0 = tic;

% initialize for the first loop
spiralPath0(:,1:accumPtNum(1)) = toolPathPt(:,1:accumPtNum(1));
spiralNorm0(:,1:accumPtNum(1)) = toolNormDirect(:,1:accumPtNum(1));
spiralCut0(:,1:accumPtNum(1)) = toolCutDirect(:,1:accumPtNum(1));
spiralQuat0(1:accumPtNum(1),:) = toolQuat(1:accumPtNum(1),:);
spiralContactU0(1:accumPtNum(1)) = toolContactU(1:accumPtNum(1));
% for each loop, shift the tool path point by decreasing the radius
for kk = 2:numLoop % begin with the 2nd loop
    angleN = toolPathAngle(accumPtNumlength(kk - 1) + 1:accumPtNumlength(kk));
    for indInterp = accumPtNum(kk - 1) + 1:accumPtNum(kk)
        % Method 1: get the (x,y) by interpolation and use residual3D to get z
        % tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
        % tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));

        % Method 2: get the inner closest point and linearly interpolate them
        angleDel = angleN - toolPathAngle(indInterp);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel); % get the closest tol point in the inner loop
        ind1 = ind1 + accumPtNumlength(kk - 1); % get the index of the closest in the whole list
        ind2 = ind2 + accumPtNumlength(kk - 1);
        [spiralPath0(:,indInterp),spiralContactU0(indInterp),spiralQuat0(indInterp,:)] = toolInterp( ...
            interpR(indInterp),toolRAccum,indInterp,ind1,ind2,toolPathPt,toolContactU,toolQuat,toolCutDirect(:,indInterp));
        spiralNorm0(:,indInterp) = quat2rotm(spiralQuat0(indInterp,:))*toolData.toolEdgeNorm;
        spiralCut0(:,indInterp) = quat2rotm(spiralQuat0(indInterp,:))*toolData.cutDirect;
        % test & debug & video
        % plot3(spiralPath(1,indInterp),spiralPath(2,indInterp),spiralPath(3,indInterp), ...
        %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
        % toolSp1 = toolSp;
        % toolSp1.coefs = quat2rotm(toolQuat(indInterp,:))*toolSp.coefs + toolVec(:,indInterp);
        % Q = fnval(toolSp1,uLim(1,indInterp):0.01:uLim(2,indInterp));
        % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
        % xlabel('x'); ylabel('y'); 
    end
end

% elliminate the first loop, which focuses on the z=0 point
spiralPath0(:,1:accumPtNum(1) - 1) = [];
spiralNorm0(:,1:accumPtNum(1) - 1) = [];
spiralCut0(:,1:accumPtNum(1) - 1) = [];
spiralQuat0(1:accumPtNum(1) - 1,:) = [];
spiralContactU0(1:accumPtNum(1) - 1) = [];

% change the spiral tool path to constant arc length one
spiralAngle0 = toolPathAngle + 2*pi*toolNAccum;
spiralAngle0(:,1:accumPtNum(1) - 1) = [];
if strcmp(angularDiscrete,'Constant Arc')
    [spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
        spiralAngle0,spiralPath0,spiralContactU0,{spiralNorm0;spiralCut0},toolData,'interpType','linear'); % spiralQuat0,
    spiralNorm = spiralVec{1};
    spiralCut = spiralVec{2};
else
    spiralAngle = spiralAngle0;
    spiralPath = spiralPath0;
    spiralNorm = spiralNorm0;
    spiralCut = spiralCut0;
    spiralQuat = spiralQuat0;
    spiralContactU = spiralContactU0;
end
spiralPtNum = length(spiralAngle);

% [tEq,fEq,vecEq] = arclengthparam(arcLength,toolThetaEach,toolREach,{},'algorithm','lq-fitting');

tSpiral = toc(tSpiral0);
fprintf('The time spent in the spiral toolpath generation process is %fs.\n',tSpiral);

%% Spiral Residual height calculation of the spiral tool path
spiralRes = 5*aimRes*ones(2,spiralPtNum);
spiralPeakPt = zeros(6,spiralPtNum);
spiralULim = [zeros(1,spiralPtNum);ones(1,spiralPtNum)];
tSpiralRes0 = tic;
parfor ind1 = 1:spiralPtNum
    % inner ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) - 2*pi,1,'first');
    ind3 = find(spiralAngle < spiralAngle(ind1) - 2*pi,1,'last');
    if isempty(ind2) || isempty(ind3)
%         ind2 = find(spiralAngle >= spiralAngle(ind1) + pi,1,'first');
%         ind3 = find(spiralAngle < spiralAngle(ind1) + pi,1,'last');
%         [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
%             spiralPath,spiralNorm,spiralCut,spiralContactU, ...
%             toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
        tmpres1 = 5*aimRes;
        ptres1 = zeros(3,1);
    else
        [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
            spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
    end

    % outer ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + 2*pi,1,'first');
    ind3 = find(spiralAngle < spiralAngle(ind1) + 2*pi,1,'last');
    if isempty(ind2) || isempty(ind3)
        tmpres2 = 5*aimRes;
        ptres2 = zeros(3,1);
    else
        [tmpres2,ptres2,spiralULim(:,ind1)] = residual3D( ...
            spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
    end
    
    spiralRes(:,ind1) = [tmpres1;tmpres2];
    spiralPeakPt(:,ind1) = [ptres1;ptres2];
    
    if spiralRes(:,ind1) > 5
        1;
    end
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

%% spiral tool path result

% get the correct direction of the spiral tool path
if strcmp(cutDirection,'Edge to Center')
    spiralPtNum = fliplr(spiralPtNum);
    spiralAngle = fliplr(spiralAngle);
    spiralPath = fliplr(spiralPath);
    spiralQuat = fliplr(spiralQuat);
    spiralNorm = fliplr(spiralNorm);
    spiralCut = fliplr(spiralCut);
    spiralContactU = fliplr(spiralContactU);
    spiralRes = fliplr(spiralRes);
    spiralULim = fliplr(spiralULim);
    spiralPeakPt = fliplr(spiralPeakPt);
end

% plot the result
figure('Name','Spiral tool path result');
tPlot0 = tic;
plotSpar = 1;
plot3(spiralPath(1,1:plotSpar:end), ...
    spiralPath(2,1:plotSpar:end), ...
    spiralPath(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
% toolCoefs = toolSp.coefs;
% parfor ii = 1:ptNum
%     toolSp1 = toolSp;
%     toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
%     Q = fnval(toolSp1,uLim(1,ii):0.05:uLim(2,ii));
%     plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
% end
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
s4_visualize_spiral;

msgfig = msgbox('Spiral tool path was generated successfully!','Success','help','non-modal');

%%
% delete(parObj);
% clear global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
%     toolQuat toolVec toolContactU isCollision ...
%     toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));
