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

    toolData = app.toolData;

    % machining paramters
    cutDirection = app.cutDirection;
    spindleDirection = app.spindleDirection;
    angularDiscrete = app.angularDiscrete;
    aimRes = app.aimRes;
    toolData = app.toolData;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = app.maxIter;
    rMax = app.rMax;
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
    % workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
    workspaceDir = 'workspace\20220925-contrast\nagayama_concentric';
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
    
    % surface data import
    default = false;
    if default
        [fileName,dirName] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','all files(*.*)'}, ...
            'Select the surface edge data file', ...
            'workspace\input_data\surface\ellipsoidAray.mat', ...
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
        % sampling density
        spar = 501;
        surfCenter = [0;0;sqrt(C^2*(R.^2-R/4.^2))]; % concentric circle center
        conR = linspace(0,R/4,spar); % concentric radius vector
        conTheta = linspace(0,2*pi,spar);
        [rMesh,thetaMesh] = meshgrid(conR,conTheta);
        surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
        surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
        surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
        % save('input_data/surface/ellipsoidAray.mat', ...
        %    "surfMesh","surfNorm","surfCenter");
    end

    % normal direction & cutting direction of the surface mesh

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

% pre-processing
surfNormSym = [surfFx;surfFy;-1];
surfNormSym = surfNormSym./norm(surfNormSym);
surfNormFunc = matlabFunction(surfNormSym,'Vars',{x,y});

switch spindleDirection
    case 'Clockwise'
        conThetaBound = [2*pi,0];
    case 'Counterclockwise'
        conThetaBound = [0,2*pi];
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
        conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
            ceil(2*pi/maxAngPtDist) + 1);
        conTheta(end) = [];
    case 'Constant Angle'
        conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
            ceil(2*pi/angularLength) + 1);
        conTheta(end) = [];
    otherwise
        error('Invalid angular discretization type.');
end

clear global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

% initialize the cutting pts: the rotation center is seen as the first loop
loopPtNum = [length(conTheta)]; % the number of points in each loop
accumPtNum = [length(conTheta)]; % the No. of the last points in each loop
toolREach = 0; % the radius R of the current loop
toolRAccum = zeros(1,loopPtNum(end));
toolNAccum = zeros(1,loopPtNum(end));
surfPt(1,:) = zeros(1,loopPtNum(end)); % x coordinates of surface cutting points
surfPt(2,:) = zeros(1,loopPtNum(end)); % y coordinates of surface cutting points
surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:)); % z coordinates of surface cutting points
surfNorm = surfNormFunc(surfPt(1,:),surfPt(2,:)); % normolization of the normal vector
surfDirect = cutdirection(surfPt,[0;0;0],conTheta,'method','vertical'); % cutting direction of each points
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
parfor ii = 1:accumPtNum(end)
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
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
    % toolSp1 = toolSp;
    % toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolSp.coefs + toolVec(:,ii);
    % Q = fnval(toolSp1,uLim(1,ii):0.01:uLim(2,ii));
    % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
    % xlabel('x'); ylabel('y'); axis equal; hold on;
end

fprintf('No.1\tElapsed time is %f seconds.\n-----\n',toc(t1));

%% Tool path adjusting for the rest
r = rStep; % 可以通过widthRes确定迭代初值
delta = rStep;
circletoolpath(r,toolData,delta,surfFunc,surfNormFunc, ...
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

save workspace\20220925-contrast\nagayama_concentric\loopR.mat toolREach

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

% cut direction change doesn't affect the concentric tool path
switch cutDirection
    case 'Edge to Center'
    case 'Center to Edge'
end

% cubic spline approximation
% the function between the numeric label of tool path and surf radius R
Fr = csape(accumPtNum,toolREach,[1,1]);

toolThetaEach = linspace(0,2*pi*(length(accumPtNum) - 1),length(accumPtNum));
rTheta = csape(toolThetaEach,toolREach,[1,1]);

figure('Name','Feed Rate Smoothing');
scatter(accumPtNum,toolREach);
hold on; grid on;
fnplt(Fr,'r',[accumPtNum(1),accumPtNum(end)]);
plot(1:accumPtNum(end),toolRAccum);
% line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
%     'Color',[0.929,0.694,0.1250]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel('Loop Accumulating Point Number');
ylabel(['Radius of the Loop (',unit,')']);
legend('No.-R scatters','csape result','Concentric result');

%% spiral tool path generation with the smoothing result
interpR = fnval(Fr,1:accumPtNum(end));
interpR(1:accumPtNum(1)) = 0;
spiralPtNum = loopPtNum;
accumPtNumlength = [0,accumPtNum];
numLoop = length(accumPtNum);

% the concentric angle of each tool path
angle0 = atan2(toolPathPt(2,:),toolPathPt(1,:));

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

% for each loop, shift the tool path point by decreasing the radius
% initialize the spiral path
spiralPath0 = zeros(3,accumPtNum(end)); % the spiral tool path
spiralNorm0 = zeros(3,accumPtNum(end));
spiralCut0 = zeros(3,accumPtNum(end));
spiralPath0(:,1:accumPtNum(1)) = toolPathPt(:,1:accumPtNum(1));
spiralNorm0(:,1:accumPtNum(1)) = toolNormDirect(:,1:accumPtNum(1));
spiralCut0(:,1:accumPtNum(1)) = toolCutDirect(:,1:accumPtNum(1));
for kk = 2:numLoop % begin with the 2nd loop
    angleN = angle0(accumPtNumlength(kk - 1) + 1:accumPtNumlength(kk));
    for indInterp = accumPtNum(kk - 1) + 1:accumPtNum(kk)
        % Method 1: get the (x,y) by interpolation and use residual3D to get z
        % tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
        % tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));

        % Method 2: get the inner closest point and linearly interpolate them
        angleDel = angleN - angle0(indInterp);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel); % get the closest tol point in the inner loop
        ind1 = ind1 + accumPtNumlength(kk - 1); % get the index of the closest in the whole list
        ind2 = ind2 + accumPtNumlength(kk - 1);
        [spiralPath0(:,indInterp),spiralNorm0(:,indInterp),spiralCut0(:,indInterp)] = toolInterp( ...
            interpR(indInterp),toolPathPt,toolNormDirect,toolCutDirect,toolRAccum,indInterp,ind1,ind2);

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

% change the spiral tool path to constant arc length one
spiralAngle = angle0 + 2*pi*toolNAccum;
if strcmp(angularDiscrete,'Constant Arc')
    [angle,spiralPath,spiralVec] = arclengthparam(arcLength, ...
        spiralAngle,spiralPath0,{spiralNorm0;spiralCut0},'interpType','linear');
    spiralNorm = spiralVec{1};
    spiralCut = spiralVec{2};
else
    spiralPath = spiralPath0;
    spiralNorm = spiralNorm0;
    spiralCut = spiralCut0;
end

% numLoop = numLoop + 1;
% accumPtNumlength = [0,accumPtNum];
% for kk = 3:numLoop % begin with the 2nd loop
%     for jj = 1:accumPtNumlength(kk)
%         % Method 1: get the (x,y) by interpolation and use residual3D to get z
%         % tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
%         % tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));
% 
%         % Method 2: get the inner closest point and linearly interpolate them
%         indInterp = accumPtNumlength(kk - 1) + jj;
%         angleN = angle(accumPtNumlength(kk - 2) + 1:accumPtNumlength(kk - 1));
%         angleDel = angleN - angleN(indInterp);
%         [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel);
%         [spiralPath(:,indInterp),spiralNorm(:,indInterp),spiralCutDir(:,indInterp)] = toolInterp( ...
%             toolPathPt,toolNormDirect,toolCutDirect,toolR,ind1,ind2,indInterp);
%     end
% end

% spiral tool path result
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
toolCoefs = toolSp.coefs;
parfor ii = 1:ptNum
    toolSp1 = toolSp;
    toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp1,uLim(1,ii):0.05:uLim(2,ii));

    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
end
% axis equal;
grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('tool center point','','tool edge','Location','northeast');
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);

% sprial tool path error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%
% delete(parObj);
clear global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));
