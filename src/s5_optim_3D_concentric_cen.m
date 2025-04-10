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
    % workspaceDir = fullfile('..','workspace','/20220925-contrast/nagayama_concentric';
    workspaceDir = fullfile('..','workspace','\20220925-contrast\nagayama_concentric';
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
            fullfile('..','workspace','\input_data\surface\ellipsoidAray.mat', ...
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
    aimRes = 50;
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

%% load the data of the residual function, and get the appropriate cutting width
% load(fullfile(workspaceDir,"widthRes.mat")); % widthRes
% wid = interp1(widthRes(2,:),widthRes(1,:),aimRes);



%% Tool path adjusting for the first loop
t0 = tic;
tRes0 = tic;
t1 = tic;

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius
toolCoefs = toolSp.coefs;

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

% initialize the cutting pts: the rotation center is seen as the first loop
loopPtNum = [length(conTheta)]; % the number of points in each loop
accumPtNum = [length(conTheta)]; % the No. of the last points in each loop
toolREach = 0; % the radius R of the current loop
toolRAccum = zeros(1,loopPtNum(end));

% tool path initialization
toolQuat = zeros(accumPtNum(end),4); % the quaternions for each points
toolContactU = zeros(1,accumPtNum(end)); % the parameter u of the contact point on the tool tip
% the xOy projection of tool path points should be arranged in loops
toolPathPt = zeros(3,accumPtNum(end));
toolCutDirect = zeros(3,accumPtNum(end));
toolNormDirect = zeros(3,accumPtNum(end));
uLim = [zeros(1,accumPtNum(end));ones(1,accumPtNum(end))]; % the interval of each toolpath
peakPt = zeros(6,accumPtNum(end));
res = 5*aimRes*ones(2,accumPtNum(end)); % the residual height, initialized with 5 times the standard aimRes

% tool contact point values
surfPt = zeros(3,loopPtNum(end)); % x coordinates of surface cutting points
surfDirect = cutdirection(surfPt,[0;0;0],conTheta,'method','vertical'); % cutting direction of each points
% surfDirect = [[0;1;0],cutdirection(surfPt(:,2:end),[0;0;0])]; % cutting direction of each points
% 和[0;0;1]垂直、且与surfDirect(:,end)夹角最小的向量作为surDirect(:,1)，可视作surfDirect(:,end)到[0;0;1]平面的投影
% surfDirect(:,1) = vecOnPlane(surfDirect(:,end),surfPt(:,1),[0;0;1]);

% debug
% figure;
% calculate the 1st loop of the tool path
% 这里可以直接只算一个然后直接旋转得到，不用每一个都算一遍的。而且为什么第一个这么慢？
parfor ii = 1:accumPtNum(end)
        [toolQuat(ii,:),toolPathPt(:,ii),toolContactU(ii),surfPt(:,ii)] = toolpathpos( ...
            toolData,toolPathPt(:,ii),[0;0;-1],surfDirect(:,ii),surfFunc,surfNormFunc);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;

        % debug
        % scatter3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii)); hold on;
        % scatter3(surfPt(1,ii),surfPt(2,ii),surfPt(3,ii));
        % qui1 = quiver3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
        %     toolCutDirect(1,ii),toolCutDirect(2,ii),toolCutDirect(3,ii), ...
        %     'AutoScale','on');
        % toolSp1 = toolSp;
        % toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolSp.coefs + toolPathPt(:,ii);
        % Q = fnval(toolSp1,0:0.01:1);
        % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
        % xlabel('x'); ylabel('y'); axis equal;
        % delete(qui1);
end

angleTmp = conTheta;
fprintf('No.1\tElapsed time is %f seconds.\n-----\n',toc(t1));

%% Tool path adjusting for the rest
r = rStep/2; % 可以通过widthRes确定迭代初值
delta = rStep;
iter = 1;

% debug
% figure('Name','concentric tool path debug');
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',1,'LineStyle','none'); hold on;
% colormap('summer');
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel('x'); ylabel('y'); zlabel('z');
% view(0,90);

while true
    tic
    while iter <= maxIter
        % calculate the discretization scatters
        switch angularDiscrete
            case 'Constant Arc'
                % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
                conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
                    ceil(2*pi/min(maxAngPtDist,arcLength/r)) + 1);
                conTheta(end) = [];
            case 'Constant Angle'
                conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
                    ceil(2*pi/angularLength) + 1);
                conTheta(end) = [];
            otherwise
                error('Invalid angular discretization type.');
        end
    
        % calculate the xOy projection of tool path points
        loopPtNumTmp = length(conTheta); % the number of points in the loop
        loopPtNumLast = loopPtNum(end);  % the number of points in the former loop
        
        % calculate the tool path and residual height
        toolPathPtTmp = zeros(3,loopPtNumTmp);
        toolPathPtTmp(1,:) = r*cos(conTheta); % x coordinates of surface cutting points in the loop
        toolPathPtTmp(2,:) = r*sin(conTheta); % y coordinates of surface cutting points
        toolQuatTmp = zeros(loopPtNumTmp,4);
        toolContactUTmp = zeros(1,loopPtNumTmp);
        toolCutDirectTmp = zeros(3,loopPtNumTmp);
        toolNormDirectTmp = zeros(3,loopPtNumTmp);

        % initialize the tool contact point values
        surfPtTmp = zeros(3,loopPtNumTmp); % the surface contact points
        surfDirectTmp = cutdirection(toolPathPtTmp,[0;0;0],conTheta,'method','vertical'); % cutting direction of each points
        for ii = 1:loopPtNumTmp
            [toolQuatTmp(ii,:),toolPathPtTmp(:,ii),toolContactUTmp(ii),surfPtTmp(:,ii)] = surfPos( ...
                toolData,toolPathPtTmp(:,ii),[0;0;-1],surfDirectTmp(:,ii),surfFunc,surfNormFunc);
            toolCutDirectTmp(:,ii) = quat2rotm(toolQuatTmp(ii,:))*toolData.cutDirect;
            toolNormDirectTmp(:,ii) = quat2rotm(toolQuatTmp(ii,:))*toolData.toolEdgeNorm;

            % debug
            plot3(toolPathPtTmp(1,ii),toolPathPtTmp(2,ii),toolPathPtTmp(3,ii), ...
                '.','MarkerSize',6,'Color',[0,0.4470,0.7410]); hold on;
            plot3(surfPtTmp(1,ii),surfPtTmp(2,ii),surfPtTmp(3,ii), ...
                '.','MarkerSize',12,'Color',[0.85,0.325,0.098]); hold on;
            qui1 = quiver3(toolPathPtTmp(1,ii),toolPathPtTmp(2,ii),toolPathPtTmp(3,ii), ...
                toolNormDirectTmp(1,ii),toolNormDirectTmp(2,ii),toolNormDirectTmp(3,ii), ...
                300,'r');
            qui2 = quiver3(toolPathPtTmp(1,ii),toolPathPtTmp(2,ii),toolPathPtTmp(3,ii), ...
                toolCutDirectTmp(1,ii),toolCutDirectTmp(2,ii),toolCutDirectTmp(3,ii), ...
                300,'b');
            toolSp1 = toolData.toolBform;
            toolSp1.coefs = quat2rotm(toolQuatTmp(ii,:))*toolData.toolBform.coefs + toolPathPtTmp(:,ii);
            Q = fnval(toolSp1,0:0.01:1);
            plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
            delete(qui1);
            delete(qui2);
        end

%         plot3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
%             '.','MarkerSize',6,'Color',[0,0.4470,0.7410]); hold on;
%         plot3(surfPtTmp(1,:),surfPtTmp(2,:),surfPtTmp(3,:), ...
%             '.','MarkerSize',12,'Color',[0.85,0.325,0.098]); hold on;
%         qui1 = quiver3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
%             toolNormDirectTmp(1,:),toolNormDirectTmp(2,:),toolNormDirectTmp(3,:), ...
%             'off','r','AutoScale','off');
%         qui2 = quiver3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
%             toolCutDirectTmp(1,:),toolCutDirectTmp(2,:),toolCutDirectTmp(3,:), ...
%             'off','b','AutoScale','off');
    
        % restore the data that will be used in the residual height calculation
        angleTmp = [angleTmp,conTheta];
        toolPathPtRes = [toolPathPt(:,end - loopPtNumLast + 1:end),toolPathPtTmp];
        toolNormDirectRes = [toolNormDirect(:,end + 1 - loopPtNumLast:end),toolNormDirectTmp];
        toolCutDirectRes = [toolCutDirect(:,end + 1 - loopPtNumLast:end),toolCutDirectTmp];
        toolContactURes = [toolContactU(end + 1 - loopPtNumLast:end),toolContactUTmp];
        uLimTmp = [uLim(:,end - loopPtNumLast + 1:end), ...
            [zeros(1,loopPtNumTmp);ones(1,loopPtNumTmp)]]; % the interval of each toolpath
        resTmp = [res(:,end - loopPtNumLast + 1:end), ...
            5*aimRes*ones(2,loopPtNumTmp)];
        peakPtInTmp = [peakPt(1:3,end - loopPtNumLast + 1:end), ...
            zeros(3,loopPtNumTmp)];
        peakPtOutTmp = zeros(3,loopPtNumLast + loopPtNumTmp);

        % calculate the residual height of the loop and the inner nearest loop
        % inner side of each point on the tool path
        angleN = angleTmp(1:loopPtNumLast);
        parfor ii = loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp
            angleDel = angleN - angleTmp(ii);
            [ind2,ind3] = getInnerLoopToolPathIndex(angleN,angleDel);
            [resTmp(1,ii),peakPtInTmp(:,ii),uLimTmp(:,ii)] = residual3D( ...
                toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
                toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
        end

        % if residual height does not satisfy the reqiurement
        maxRes = max(resTmp(1,loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp),[],"all");
        if maxRes < aimRes
            break;
        else
            fprintf('\tIter %d: The maximum residual height is %f %cm.\n',iter,maxRes,char(956));
            delta = delta/3;
            r = r - delta;
            iter = iter + 1;
        end
    end

    if iter == maxIter
        fprintf('The iteration ends with the maxIter is reached. \n');
        fprintf('There is somewhere that has no intersection with the inside loop of toolpath.\n');
        return
    end

    % outer side of each point in the tool path3
    % 这里，内圈的计算，出了问题。为什么算来算去都是250？？？？？？？？？
    angleN = angleTmp(loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp);
    parfor ii = 1:loopPtNumLast
        angleDel = angleN - angleTmp(ii);
        [ind2,ind3] = getOuterLoopToolPathIndex(angleN,angleDel,loopPtNumLast);
        [resTmp(2,ii),peakPtOutTmp(:,ii),uLimTmp(:,ii)] = residual3D( ...
            toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
            toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
    end
    
    % then store the data of this loop
    loopPtNum = [loopPtNum,loopPtNumTmp];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNumTmp];
    toolREach = [toolREach,r];
    toolRAccum = [toolRAccum,r*ones(1,loopPtNumTmp)];
    surfPt = [surfPt,surfPtTmp];
    % surfDirect = [surfDirect,surfDirectTmp];
    toolQuat = [toolQuat;toolQuatTmp];
    toolContactU = [toolContactU,toolContactUTmp];
    toolPathPt = [toolPathPt,toolPathPtTmp];
    toolNormDirect = [toolNormDirect,toolNormDirectTmp];
    toolCutDirect = [toolCutDirect,toolCutDirectTmp];

    uLim(:,end - loopPtNumLast + 1:end) = [];
    uLim = [uLim,uLimTmp]; % the interval of each toolpath
    peakPt(:,end - loopPtNumLast + 1:end) = [];
    peakPt = [peakPt,[peakPtInTmp;peakPtOutTmp]];
    res(:,end - loopPtNumLast + 1:end) = [];
    res = [res,resTmp];

    % debug
%     plot3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
%         '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
%     quiver3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
%         toolNormDirectTmp(1,:),toolNormDirectTmp(2,:),toolNormDirectTmp(3,:));
%     for jj = accumPtNum(end - 1) + 1:accumPtNum(end)
%         toolSp1 = toolSp;
%         toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolPathPt(:,jj);
%         Q = fnval(toolSp1,0:0.01:1);
%         plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
%     end

    % clear the temporary variables
    angleTmp = conTheta;
    clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp ...
        toolContactUTmp toolPathPtTmp toolNormDirectTmp ...
        toolCutDirectTmp;
    clear uLimTmp peakPtOutTmp resTmp loopPtNumLast;
    fprintf('No.%d\tElapsed time is %f seconds.\n-----\n',length(loopPtNum),toc);
    if r > rMax, break; end
    delta = rStep;
    r = r + delta;
    iter = 1;
end
tRes = toc(tRes0);

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
        toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolPathPt(:,jj);
        Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
        % isQDomain = (Q(1,:).^2 + Q(2,:).^2 - surfDomain(1,2)^2) >= 0;
        % Q(:,isQDomain) = [];
        plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1); hold on;
    end
else
    for jj = 1:ptNum
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolPathPt(:,jj);
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

save ..\workspace\20220925-contrast\nagayama_concentric\loopR.mat toolREach

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
% accumPtNumlength = [0,accumPtNum];
% loopRLength = [0,loopR];
% Fr = csape(accumPtNumlength,loopRLength); 

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

% initialize the spiral path
spiralPath = zeros(3,accumPtNum(end)); % the spiral tool path
spiralNorm = zeros(3,accumPtNum(end));
spiralCutDir = zeros(3,accumPtNum(end));
spiralPath(:,1:accumPtNum(1)) = toolPathPt(:,1:accumPtNum(1));
spiralNorm(:,1:accumPtNum(1)) = toolNormDirect(:,1:accumPtNum(1));
spiralCutDir(:,1:accumPtNum(1)) = toolCutDirect(:,1:accumPtNum(1));

% the concentric angle of each tool path
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));

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
for kk = 2:numLoop % begin with the 2nd loop
    angleN = angle(accumPtNumlength(kk - 1) + 1:accumPtNumlength(kk));
    for indInterp = accumPtNum(kk - 1) + 1:accumPtNum(kk)
        % Method 1: get the (x,y) by interpolation and use residual3D to get z
        % tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
        % tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));

        % Method 2: get the inner closest point and linearly interpolate them
        angleDel = angleN - angle(indInterp);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel); % get the closest tol point in the inner loop
        ind1 = ind1 + accumPtNumlength(kk - 1); % get the index of the closest in the whole list
        ind2 = ind2 + accumPtNumlength(kk - 1);
        [spiralPath(:,indInterp),spiralNorm(:,indInterp),spiralCutDir(:,indInterp)] = toolInterp( ...
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

%% Comparison: directly generate the spiral tool path
% 实际上，这种显然更好。用上面的那种方法，会导致不是真正的等弧长，而且在交接段会突然减速，动力学应该会影响表面质量
% parfor ii = (sparTheta + 1):ptNum
%     % 如果是沿同一个极径的，就可以直接不用投影；否则还是需要这样子找
%     nLoop = floor((ii - 1)/sparTheta) - 1;
%     angleN = angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1));
%     angleDel = angleN - angle(ii);
%     % ind2(ii) remains the index of angle nearest to angle(ii) within 
%     % those which is larger than the angle(ii) and in angleN
%     if isempty(angleN(angleDel >= 0))
%         % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
%         angleDel = angleDel + 2*pi;
%     end
%     ind2 = sparTheta*nLoop + find(angleN == min(angleN(angleDel >= 0)));
%     % ind3(ii) remains the index of angle nearest to angle(ii) within 
%     % those which is smaller than the angle(ii) and in angleN
%     if isempty(angleN(angleDel < 0))
%         angleDel = angleDel - 2*pi;
%     end
%     ind3 = sparTheta*nLoop + find(angleN == max(angleN(angleDel < 0)));
%     [res(resNum + ii),peakPt(:,resNum + ii),uLim(:,ii)] = residual3D( ...
%         toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp,toolRadius, ...
%         uLim(:,ii),ii,ind2,ind3);
% end
% 
% 
% 
% 
%         spiralPath(:,accumPtNum(kk) + jj) = tmpSpiral;
%         spiralNorm(:,accumPtNum(kk) + jj) = ;
%     end
% end




%%
% delete(parObj);
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));