% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath
% close all;
clear;
clc;
addpath(genpath('funcs'));
t0 = tic;

% global variables
% global textFontSize textFontType;
workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
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
else % ellipsoid
    R = 10/2*1000;
    A = 3.5/2;
    B = 4/2;
    C = 5/2;
    % sampling density
    spar = 101;
    surfCenter = [0;0;sqrt(C^2*(R.^2-R/4.^2))]; % concentric circle center
    conR = linspace(0,R/4,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = sqrt(C^2*(R.^2 - rMesh.^2));
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

% machining surface
syms x y;
surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
surfFunc = matlabFunction(surfSym);
surfFx = matlabFunction(diff(surfFunc,x));
surfFy = matlabFunction(diff(surfFunc,y));
rMax = R/2;

% machining paramters
aimRes = 10;
rStep = toolData.radius; % 每步步长可通过曲面轴向偏导数确定
arcLength = 30;
maxAngPtDist = 6*pi/180;

%% load the data of the residual function, and get the appropriate cutting width
load(fullfile(workspaceDir,"widthRes.mat")); % widthRes
% wid = interp1(widthRes(2,:),widthRes(1,:),aimRes);



%% Tool path adjusting for the first loop
tRes0 = tic;
t1 = tic;
r = rStep/2; % 可以通过widthRes确定迭代初值
delta = rStep;
while true
    % calculate the discretization scatters
    % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
    conTheta = linspace(0,2*pi, ...
        ceil(2*pi/min(maxAngPtDist,arcLength/r)) + 1);
    conTheta(end) = [];
    
    % calculate the cutting pts
    loopPtNum = length(conTheta); % the number of points in each loop
    loopR = r; % the radius R of the current loop
    surfPt(1,:) = r*cos(conTheta); % x coordinates of surface cutting points
    surfPt(2,:) = r*sin(conTheta); % y coordinates of surface cutting points
    surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:)); % z coordinates of surface cutting points
    surfNorm(1,:) = surfFx(surfPt(1,:),surfPt(2,:)); % x coordinates of the normal vectors on surface points
    surfNorm(2,:) = surfFy(surfPt(1,:),surfPt(2,:)); % y coordinates of the normal vectors on surface points
    surfNorm(3,:) = -1*ones(1,loopPtNum); % z coordinates of the normal vectors on surface points
    surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1)); % normolization of the normal vector
    surfDirect = cutDirection(surfPt,[0;0;0]); % cutting direction of each points
    
    % tool path initialization
    toolSp = toolData.toolBform; % B-form tool tip arc
    toolRadius = toolData.radius; % fitted tool radius
    toolQuat = zeros(loopPtNum,4); % the quaternions for each points
    toolVec = zeros(3,loopPtNum); % the translation vectors for each points
    toolContactU = zeros(1,loopPtNum); % the parameter u of the contact point on the tool tip
    isCollision = false(1,loopPtNum);
    toolPathPt = zeros(3,loopPtNum);
    toolCutDirect = zeros(3,loopPtNum);
    toolNormDirect = zeros(3,loopPtNum);
    uLim = [zeros(1,loopPtNum);ones(1,loopPtNum)]; % the interval of each toolpath
    peakPtIn = zeros(3,loopPtNum);
    res = 5*aimRes*ones(2,loopPtNum);
    
    % calculate the tool path and residual height
    parfor ii = 1:loopPtNum
        [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
            toolData,surfPt(:,ii),surfNorm(:,ii),[0;0;1],surfDirect(:,ii));
        if isCollision(ii) == false
            toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
            toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
            toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
        end
    end
    
    angle = atan2(toolPathPt(2,:),toolPathPt(1,:));
    for ii = 1:loopPtNum
        angleDel = (angle - wrapToPi(angle(ii) + pi));
        if isempty(angle(angleDel >= 0))
            % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
            angleDel = angleDel + 2*pi;
        end
        ind2 = find(angle == min(angle(angleDel >= 0)));
        if isempty(angle(angleDel < 0))
            angleDel = angleDel - 2*pi;
        end
        ind3 = find(angle == max(angle(angleDel < 0)));
        [res(1,ii),peakPtIn(:,ii),uLim(:,ii)] = residual3D( ...
            toolPathPt,toolNormDirect,toolCutDirect,toolContactU, ...
            toolSp,toolRadius,uLim(:,ii),ii,ind2,ind3);
    end
    
    if max(res(1,:)) < aimRes
        break;
    else
        warning('The residual height of No.%d is beyond the expected range.',length(loopPtNum));
        delta = delta/3;
        r = r - delta;
    end
end
peakPt = [peakPtIn;zeros(3,loopPtNum)];
r = r + delta;
fprintf('No.1\tElapsed time is %f seconds.\n',toc(t1));

%% Tool path adjusting for the rest
while true
    tic
    while true
        % calculate the discretization scatters
        % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
        conTheta = linspace(0,2*pi, ...
            ceil(2*pi/min(maxAngPtDist,arcLength/r)) + 1);
        conTheta(end) = [];
    
        % calculate the cutting pts
        loopPtNumTmp = length(conTheta);
        loopPtNumLast = loopPtNum(end);
        surfPtTmp(1,:) = r*cos(conTheta);
        surfPtTmp(2,:) = r*sin(conTheta);
        surfPtTmp(3,:) = surfFunc(surfPtTmp(1,:),surfPtTmp(2,:));
        surfNormTmp(1,:) = surfFx(surfPtTmp(1,:),surfPtTmp(2,:));
        surfNormTmp(2,:) = surfFy(surfPtTmp(1,:),surfPtTmp(2,:));
        surfNormTmp(3,:) = -1*ones(1,loopPtNumTmp);
        surfNormTmp = -1*(surfNormTmp./vecnorm(surfNormTmp,2,1));
        surfDirectTmp = cutDirection(surfPtTmp,[0;0;0]);
        
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
        peakPtTmpIn = [peakPt(1:3,end - loopPtNumLast + 1:end), ...
            zeros(3,loopPtNumTmp)];
        peakPtTmpOut = zeros(3,loopPtNumLast + loopPtNumTmp);
    
        % calculate the residual height of the loop and the inner nearest loop
        angle = atan2(toolPathPtRes(2,:),toolPathPtRes(1,:));
        % inner side of each point on the tool path
        angleN = angle(1:loopPtNumLast);
        for ii = loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp
            angleDel = angleN - angle(ii);
            if isempty(angleN(angleDel >= 0))
                % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
                angleDel = angleDel + 2*pi;
            end
            ind2 = find(angleN == min(angleN(angleDel >= 0)));
            if isempty(angleN(angleDel < 0))
                angleDel = angleDel - 2*pi;
            end
            ind3 = find(angleN == max(angleN(angleDel < 0)));
            [resTmp(1,ii),peakPtTmpIn(:,ii),uLimTmp(:,ii)] = residual3D( ...
                toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
                toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
        end

        % if residual height satisfies the reqiurement,
        if max(resTmp(1,loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp),[],"all") < aimRes
            break;
        else
            warning('The residual height of No.%d is beyond the expected range.',length(loopPtNum) + 1);
            delta = delta/3;
            r = r - delta;
                clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp toolVecTmp ...
                    toolContactUTmp isCollisionTmp toolPathPtTmp toolNormDirectTmp ...
                    toolCutDirectTmp;
                clear uLimTmp peakPtTmpOut resTmp loopPtNumLast;
        end
    end

    % outer side of each point in the tool path
    angleN = angle(loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp);
    for ii = 1:loopPtNumLast
        angleDel = angleN - angle(ii);
        if isempty(angleN(angleDel >= 0))
            % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
            angleDel = angleDel + 2*pi;
        end
        ind2 = loopPtNumLast + find(angleN == min(angleN(angleDel >= 0)));
        if isempty(angleN(angleDel < 0))
            angleDel = angleDel - 2*pi;
        end
        ind3 = loopPtNumLast + find(angleN == max(angleN(angleDel < 0)));
        [resTmp(2,ii),peakPtTmpOut(:,ii),uLimTmp(:,ii)] = residual3D( ...
            toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
            toolSp,toolRadius,uLimTmp(:,ii),ii,ind2,ind3);
    end
    
    % then store the data of this loop
    loopPtNum = [loopPtNum,loopPtNumTmp];
    loopR = [loopR,r];
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
    peakPt = [peakPt,[peakPtTmpIn;peakPtTmpOut]];
    res(:,end - loopPtNumLast + 1:end) = [];
    res = [res,resTmp];
    
    clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp toolVecTmp ...
        toolContactUTmp isCollisionTmp toolPathPtTmp toolNormDirectTmp ...
        toolCutDirectTmp;
    clear uLimTmp peakPtTmpOut resTmp loopPtNumLast;
    fprintf('No.%d\tElapsed time is %f seconds.\n',length(loopPtNum),toc);
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

s5_visualize_process;

%% Feed rate smoothing
% to smooth the loopR to get the real tool path

% cubic spline approximation
accumPtNum = loopPtNum;
for kk = 2:length(loopPtNum)
    accumPtNum(kk) = loopPtNum(kk) + accumPtNum(kk - 1);
end
accumPtNum = [0,accumPtNum];
loopRcsape = [0,loopR];
Fr = csape(accumPtNum,loopRcsape);

figure('Name','Feed Rate Smoothing');
scatter(accumPtNum,loopRcsape);
hold on; grid on;
fnplt(Fr,'r',[0,accumPtNum(end)]);
% line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
%     'Color',[0.929,0.694,0.1250]);
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel('Loop Accumulating Point Number');
ylabel(['Radius of the Loop (',unit,')']);

% tool path generation with the smoothing result
spiralPtNum = loopPtNum;
spiralPath = zeros(3,size(toolPathPt,2)); % the spiral tool path
spiralNorm = zeros(3,size(toolPathPt,2));
for kk = 1:length(loopPtNum) % each loop
    for jj = 1:loopPtNum(kk)
        tmpPt = toolPathPt(1:2,accumPtNum(kk) + jj);
        tmpSpiral = tmpPt + tmpPt/norm(tmpPt)*(loopR(kk) - fnval(Fr,accumPtNum(kk) + jj));
        
        spiralPath(:,accumPtNum(kk) + jj) = tmpSpiral;
        spiralNorm(:,accumPtNum(kk) + jj) = ;
    end
end




%%
% delete(parObj);
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));