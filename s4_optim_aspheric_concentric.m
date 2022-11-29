% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

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
    rStep = toolData.radius; % 每步步长可通过曲面轴向偏导数确定
    rMax = app.rMax;
    arcLength = app.arcLength;
    maxAngPtDist = app.maxAngPtDist;
    angularLength = app.angularLength;

    surfFunc = app.surfFuncs;
    surfFx = app.surfFx;
    surfFy = app.surfFy;
    surfDomain = app.surfDomain;
    surfMesh = app.surfMesh;
    
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
    workspaceDir = 'workspace\20221020-tooltip';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    
    msgOpts.Default = 'Cancel and quit';
    msgOpts.Interpreter = 'tex';
    % msgOpts.modal = 'non-modal';
    profile on
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
        R = 10/2*1000;
        A = 3.5/2;
        B = 4/2;
        C = 5/2;
        % machining surface
        syms x y;
        surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
        surfFunc = matlabFunction(surfSym);
        surfFx = matlabFunction(diff(surfFunc,x));
        surfFy = matlabFunction(diff(surfFunc,y));
        rMax = R/2;
    else % ellipsoid
        R = 10/2*1000;
        A = 3.5/2;
        B = 4/2;
        C = 5/2;
        % machining surface
        syms x y;
        surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
        surfFunc = matlabFunction(surfSym);
        surfFx = matlabFunction(diff(surfFunc,x));
        surfFy = matlabFunction(diff(surfFunc,y));
        rMax = R/2;
        % sampling density
        spar = 501;
        surfCenter = [0;0;sqrt(C^2*(R.^2-R/4.^2))]; % concentric circle center
        conR = linspace(0,R/4,spar); % concentric radius vector
        conTheta = linspace(0,2*pi,spar);
        [rMesh,thetaMesh] = meshgrid(conR,conTheta);
        surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
        surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
        surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
    %     [surfNorm(:,:,1),surfNorm(:,:,2),surfNorm(:,:,3)] = surfnorm( ...
    %         surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));
        % save('input_data/surface/ellipsoidAray.mat', ...
        %    "surfMesh","surfNorm","surfCenter");
    end
    
    figure('Name','original xyz scatters of the surface (sparsely)');
    surf( ...
        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
    hold on; axis equal;
    set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
    xlabel(['x (',unit,')']);
    ylabel(['y (',unit,')']);
    zlabel(['z (',unit,')']);
    msgfig = questdlg({'Surface was generated successfully!', ...
        'Ready to continue?'}, ...
        'Surface Generation','OK & continue','Cancel & quit','OK & continue');
    if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
        msgbox('Exit for the program','Exit','error','modal');
        uiwait(msgbox);
        return;
    end
    
    % machining paramters    
    cutDirection = 'Center to Edge'; % 'Edge to Center'
    spindleDirection = 'Counterclockwise'; % 'Clockwise'
    angularDiscrete = 'Constant Arc'; % 'Constant Angle'
    aimRes = 50;
    rStep = toolData.radius; % 每步步长可通过曲面轴向偏导数确定
    arcLength = 30;
    maxAngPtDist = 6*pi/180;
    angularLength = 6*pi/180;
end

%% load the data of the residual function, and get the appropriate cutting width
load(fullfile(workspaceDir,"widthRes.mat")); % widthRes
% wid = interp1(widthRes(2,:),widthRes(1,:),aimRes);



%% Tool path adjusting for the first loop
t0 = tic;
tRes0 = tic;
t1 = tic;
r = rStep/2; % 可以通过widthRes确定迭代初值
delta = rStep;

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius

while true
    % calculate the discretization scatters
    % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
    switch angularDiscrete
        case 'Constant Arc'
            conTheta = linspace(0,2*pi, ...
                ceil(2*pi/maxAngPtDist) + 1);
            conTheta(1) = [];
        case 'Constant Angle'
            conTheta = linspace(0,2*pi, ...
                ceil(2*pi/angularLength) + 1);
            conTheta(1) = [];
        otherwise
            error('Invalid angular discretization type.');
    end

%     % initialize the cutting pts
%     loopPtNum = 1; % the number of points in each loop
%     loopR = 0; % the radius R of the current loop
%     toolR = 0; % the radius R of each tool path points
%     surfPt = [0;0;surfFunc(0,0)]; % coordinates of surface cutting points
%     surfNorm = [surfFx(0,0);surfFy(0,0);-1]; % coordinates of the normal vectors on surface points
%     surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1)); % normolization of the normal vector
%     surfDirect = [1;0;0]; % cutting direction of each points


    % initialize the cutting pts: the rotation center is seen as the first loop
    loopPtNum = [1,length(conTheta)]; % the number of points in each loop
    accumPtNum = [1,1 + length(conTheta)]; % the No. of the last points in each loop
    loopR = [0,r]; % the radius R of the current loop
    toolR = [0,r*ones(1,loopPtNum(end))];
    surfPt(1,:) = [0,r*cos(conTheta)]; % x coordinates of surface cutting points
    surfPt(2,:) = [0,r*sin(conTheta)]; % y coordinates of surface cutting points
    surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:)); % z coordinates of surface cutting points
    surfNorm(1,:) = surfFx(surfPt(1,:),surfPt(2,:)); % x coordinates of the normal vectors on surface points
    surfNorm(2,:) = surfFy(surfPt(1,:),surfPt(2,:)); % y coordinates of the normal vectors on surface points
    surfNorm(3,:) = -1*ones(1,accumPtNum(end)); % z coordinates of the normal vectors on surface points
    surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1)); % normolization of the normal vector
    surfDirect = [[0;1;0],cutdirection(surfPt(:,2:end),[0;0;0])]; % cutting direction of each points
    % 和[0;0;1]垂直、且与surfDirect(:,end)夹角最小的向量作为surDirect(:,1)，可视作surfDirect(:,end)到[0;0;1]平面的投影
    surfDirect(:,1) = vecOnPlane(surfDirect(:,end),surfPt(:,1),[0;0;1]);

    % calculate the tool path and residual height
    [tQuat(1,:),tVec(:,1),tContactU(1),isColl(1)] = toolPos( ...
        toolData,surfPt(:,1),surfNorm(:,1),[0;0;1],surfDirect(:,1));
    if isColl(1) == false
        tPathPt(:,1) = quat2rotm(tQuat(1,:))*toolData.center + tVec(:,1);
        tCutDir(:,1) = quat2rotm(tQuat(1,:))*toolData.cutDirect;
        tNorm(:,1) = quat2rotm(tQuat(1,:))*toolData.toolEdgeNorm;
    end
    [tQuat(2,:),tVec(:,2),tContactU(2),isColl(2)] = toolPos( ...
        toolData,surfPt(:,end),surfNorm(:,end),[0;0;1],surfDirect(:,end));
    if isColl(2) == false
        tPathPt(:,2) = quat2rotm(tQuat(2,:))*toolData.center + tVec(:,2);
        tCutDir(:,2) = quat2rotm(tQuat(2,:))*toolData.cutDirect;
        tNorm(:,2) = quat2rotm(tQuat(2,:))*toolData.toolEdgeNorm;
    end

    tSp1 = toolSp;
    tSp1.coefs = axesRot([0;0;1],[1;0;0],tNorm(:,1),tCutDir(:,1),'zx')*toolSp.coefs + tPathPt(:,1);
    tSp2 = toolSp;
    tSp2.coefs = axesRot([0;0;1],[1;0;0],tNorm(:,2),tCutDir(:,2),'zx')*toolSp.coefs + tPathPt(:,2);
    % figure;
    % plot3(tSp1.coefs(1,:),tSp1.coefs(2,:),tSp1.coefs(3,:),'.'); hold on;
    % plot3(tSp2.coefs(1,:),tSp2.coefs(2,:),tSp2.coefs(3,:),'.');
    
    [res,~] = residual2D_numeric(tSp1,tSp2,1e-3,fnval(tSp1,tContactU(1)),fnval(tSp2,tContactU(2)),'DSearchn');

    if res < aimRes
        clear tQuat tVec tContactU isColl tPathPt tCutDir tNorm res
        break;
    else
        fprintf('The residual height of No.%d is beyond the expected range.\n-----\n',length(loopPtNum));
        delta = delta/3;
        r = r - delta;
    end
end

% tool path initialization
toolQuat = zeros(accumPtNum(end),4); % the quaternions for each points
toolVec = zeros(3,accumPtNum(end)); % the translation vectors for each points
toolContactU = zeros(1,accumPtNum(end)); % the parameter u of the contact point on the tool tip
isCollision = false(1,accumPtNum(end));
toolPathPt = zeros(3,accumPtNum(end));
toolCutDirect = zeros(3,accumPtNum(end));
toolNormDirect = zeros(3,accumPtNum(end));
uLim = [zeros(1,accumPtNum(end));ones(1,accumPtNum(end))]; % the interval of each toolpath
peakPtIn = zeros(3,accumPtNum(end));
res = 5*aimRes*ones(2,accumPtNum(end)); % the residual height, initialized with 5 times the standard aimRes

% calculate the 1st & 2nd loop of the tool path
parfor ii = 1:accumPtNum(end)
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),[0;0;1],surfDirect(:,ii));
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
end

% calculate the residual height together with the rotation center
angle = atan2(toolPathPt(2,2:accumPtNum(end)),toolPathPt(1,2:accumPtNum(end)));
for ii = 2:accumPtNum(end)
    angleDel = (angle - wrapToPi(angle(ii - 1) + pi));
    [ind2,ind3] = getInnerLoopToolPathIndex(angle,angleDel);
    [res(1,ii),peakPtIn(:,ii),uLim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU, ...
        toolSp,toolRadius,uLim(:,ii),ii,ind2,ind3);
end
peakPt = [peakPtIn;zeros(3,accumPtNum(end))];


r = r + delta;
fprintf('No.2\tElapsed time is %f seconds.\n',toc(t1));

%% Tool path adjusting for the rest
while true
    tic
    while true
        % calculate the discretization scatters
        switch angularDiscrete
            case 'Constant Arc'
                % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
                conTheta = linspace(0,2*pi, ...
                    ceil(2*pi/min(maxAngPtDist,arcLength/r)) + 1);
                conTheta(end) = [];
            case 'Constant Angle'
                conTheta = linspace(0,2*pi, ...
                    ceil(2*pi/angularLength) + 1);
                conTheta(1) = [];
            otherwise
                error('Invalid angular discretization type.');
        end
    
        % calculate the cutting pts
        loopPtNumTmp = length(conTheta); % the number of points in the loop
        loopPtNumLast = loopPtNum(end);  % the number of points in the former loop
        surfPtTmp(1,:) = r*cos(conTheta); % x coordinates of surface cutting points in the loop
        surfPtTmp(2,:) = r*sin(conTheta); % y coordinates of surface cutting points
        surfPtTmp(3,:) = surfFunc(surfPtTmp(1,:),surfPtTmp(2,:)); % z coordinates of surface cutting points
        surfNormTmp(1,:) = surfFx(surfPtTmp(1,:),surfPtTmp(2,:)); % x coordinates of the normal vectors on surface points
        surfNormTmp(2,:) = surfFy(surfPtTmp(1,:),surfPtTmp(2,:)); % y coordinates of the normal vectors on surface points
        surfNormTmp(3,:) = -1*ones(1,loopPtNumTmp); % z coordinates of the normal vectors on surface points
        surfNormTmp = -1*(surfNormTmp./vecnorm(surfNormTmp,2,1)); % normolization of the normal vector
        surfDirectTmp = cutdirection(surfPtTmp,[0;0;0]); % cutting direction of each points
        
        % calculate the tool path and residual height
        toolQuatTmp = zeros(loopPtNumTmp,4);
        toolVecTmp = zeros(3,loopPtNumTmp);
        toolContactUTmp = zeros(1,loopPtNumTmp);
        isCollisionTmp = zeros(1,loopPtNumTmp);
        toolPathPtTmp = zeros(3,loopPtNumTmp);
        toolCutDirectTmp = zeros(3,loopPtNumTmp);
        toolNormDirectTmp = zeros(3,loopPtNumTmp);
        for ii = 1:loopPtNumTmp
            [toolQuatTmp(ii,:),toolVecTmp(:,ii),toolContactUTmp(ii)] = toolPos( ...
                toolData,surfPtTmp(:,ii),surfNormTmp(:,ii),[0;0;1],surfDirectTmp(:,ii));
            toolPathPtTmp(:,ii) = quat2rotm(toolQuatTmp(ii,:))*toolData.center + toolVecTmp(:,ii);
            toolCutDirectTmp(:,ii) = quat2rotm(toolQuatTmp(ii,:))*toolData.cutDirect;
            toolNormDirectTmp(:,ii) = quat2rotm(toolQuatTmp(ii,:))*toolData.toolEdgeNorm;
        end
    
        % restore the data that will be used in the residual height calculation
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
        angle = atan2(toolPathPtRes(2,:),toolPathPtRes(1,:));
        % inner side of each point on the tool path
        angleN = angle(1:loopPtNumLast);
        parfor ii = loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp
            angleDel = angleN - angle(ii);
            [ind2,ind3] = getInnerLoopToolPathIndex(angleN,angleDel);
            % if isempty(angleN(angleDel >= 0))
            %     % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
            %     angleDel = angleDel + 2*pi;
            % end
            % ind2 = find(angleN == min(angleN(angleDel >= 0)));
            % if isempty(angleN(angleDel < 0))
            %     angleDel = angleDel - 2*pi;
            % end
            % ind3 = find(angleN == max(angleN(angleDel < 0)));
            [resTmp(1,ii),peakPtInTmp(:,ii),uLimTmp(:,ii)] = residual3D( ...
                toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
                toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
        end

        % if residual height does not satisfy the reqiurement,
        if max(resTmp(1,loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp),[],"all") < aimRes
            break;
        else
            fprintf('The residual height of No.%d is beyond the expected range.\n',length(loopPtNum) + 1);
            delta = delta/3;
            r = r - delta;
                clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp toolVecTmp ...
                    toolContactUTmp isCollisionTmp toolPathPtTmp toolNormDirectTmp ...
                    toolCutDirectTmp;
                clear uLimTmp peakPtOutTmp resTmp loopPtNumLast;
        end
    end

    % outer side of each point in the tool path
    angleN = angle(loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp);
    parfor ii = 1:loopPtNumLast
        angleDel = angleN - angle(ii);
        [ind2,ind3] = getOuterLoopToolPathIndex(angleN,angleDel,loopPtNumLast);
        % if isempty(angleN(angleDel >= 0))
        %     % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
        %     angleDel = angleDel + 2*pi;
        % end
        % ind2 = loopPtNumLast + find(angleN == min(angleN(angleDel >= 0)));
        % if isempty(angleN(angleDel < 0))
        %     angleDel = angleDel - 2*pi;
        % end
        % ind3 = loopPtNumLast + find(angleN == max(angleN(angleDel < 0)));
        [resTmp(2,ii),peakPtOutTmp(:,ii),uLimTmp(:,ii)] = residual3D( ...
            toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
            toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
    end
    
    % then store the data of this loop
    loopPtNum = [loopPtNum,loopPtNumTmp];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNumTmp];
    loopR = [loopR,r];
    toolR = [toolR,r*ones(1,loopPtNumTmp)];
    surfPt = [surfPt,surfPtTmp];
    % surfNorm = [surfNorm,surfNormTmp];
    % surfDirect = [surfDirect,surfDirectTmp];
    toolQuat = [toolQuat;toolQuatTmp];
    toolVec = [toolVec,toolVecTmp];
    toolContactU = [toolContactU,toolContactUTmp];
    isCollision = [isCollision,isCollisionTmp];
    toolPathPt = [toolPathPt,toolPathPtTmp];
    toolNormDirect = [toolNormDirect,toolNormDirectTmp];
    toolCutDirect = [toolCutDirect,toolCutDirectTmp];

    uLim(:,end - loopPtNumLast + 1:end) = [];
    uLim = [uLim,uLimTmp]; % the interval of each toolpath
    peakPt(:,end - loopPtNumLast + 1:end) = [];
    peakPt = [peakPt,[peakPtInTmp;peakPtOutTmp]];
    res(:,end - loopPtNumLast + 1:end) = [];
    res = [res,resTmp];
    
    clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp toolVecTmp ...
        toolContactUTmp isCollisionTmp toolPathPtTmp toolNormDirectTmp ...
        toolCutDirectTmp;
    clear uLimTmp peakPtOutTmp resTmp loopPtNumLast;
    fprintf('No.%d\tElapsed time is %f seconds.\n-----\n',length(loopPtNum),toc);
    if r > rMax, break; end
    delta = rStep;
    r = r + delta;
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
for ii = 1:ptNum
    toolSp = toolData.toolBform;
    toolCoefs = toolData.toolBform.coefs;
    toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp,uLim(1,ii):0.01:uLim(2,ii));
    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
end
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
legend('designed surface','tool center point','tool edge','Location','northeast');
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);

save workspace\20220925-contrast\nagayama_concentric\loopR.mat loopR

s4_visualize_concentric;
% 这里还有问题：中间一圈的残高比外圈的小很多，而且中间的一圈残高的计算是有问题的。

%% Feed rate smoothing
% to smooth the loopR to get the real tool path

% cubic spline approximation
% the function between the numeric label of tool path and surf radius R
Fr = csape(accumPtNum,loopR); 

figure('Name','Feed Rate Smoothing');
scatter(accumPtNum,loopR);
hold on; grid on;
fnplt(Fr,'r',[0,accumPtNum(end)]);
plot(1:accumPtNum(end),toolR);
% line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
%     'Color',[0.929,0.694,0.1250]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel('Loop Accumulating Point Number');
ylabel(['Radius of the Loop (',unit,')']);
legend('No.-R scatters','csape result','Concentric result');

% tool path generation with the smoothing result
interpR = fnval(Fr,1:accumPtNum(end));
spiralPtNum = loopPtNum;
numLoop = length(accumPtNum);
spiralPath = zeros(3,size(toolPathPt,2)); % the spiral tool path
spiralNorm = zeros(3,size(toolPathPt,2));
spiralCutDir = zeros(3,size(toolPathPt,2));
spiralPath(:,1) = toolPathPt(:,1);
spiralNorm(:,1) = toolNormDirect(:,1);
spiralCutDir(:,1) = toolCutDirect(:,1);
angle = atan2(toolPathPt(2,:),toolPathPt(1,:)); % the concentric angle of each tool path
% for each loop, shift the tool path point by decreasing the radius
for kk = 2:numLoop % begin with the 2nd loop
    for jj = 1:spiralPtNum(kk)
        % Method 1: get the (x,y) by interpolation and use residual3D to get z
        % tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
        % tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));

        % Method 2: get the inner closest point and linearly interpolate them
        angleN = angle(accumPtNum(kk - 1) + 1:accumPtNum(kk));
        angleDel = angleN - angleN(jj);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel);
        ind1 = ind1 + accumPtNum(kk - 1);
        ind2 = ind2 + accumPtNum(kk - 1);
        indInterp = jj + accumPtNum(kk - 1);
        [spiralPath(:,indInterp),spiralNorm(:,indInterp),spiralCutDir(:,indInterp)] = toolInterp( ...
            toolPathPt,toolNormDirect,toolCutDirect,interpR,ind1,ind2,indInterp);
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
plotSpar = 5;
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
for ii = 1:ptNum
    toolSp = toolData.toolBform;
    toolCoefs = toolData.toolBform.coefs;
    toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp,uLim(1,ii):0.05:uLim(2,ii));
    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
end
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('tool center point','','tool edge','Location','northeast');
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);

% sprial tool path error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Comparison: directly generate the spiral tool path
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