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
%     tPar0 = tic;
%     parObj = gcp;
%     tPar = toc(tPar0);
%     fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);
    
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
    A = tand(20)/(2*2000);
    C = -2*2000*A;
    syms x y;
    surfSym = A*(x.^2 + y.^2) + C;
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-1000,1000;-1000,1000];
    rMax = max(surfDomain(1,2),surfDomain(2,2));
    % sampling density
    spar = 501;
    conR = linspace(0,rMax,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
    % save('input_data/surface/ellipsoidAray.mat', ...
    %    "surfMesh","surfNorm","surfCenter");
    
    surfType = '2D';
    Geometry2DCell = {'Rotating Paraboloid','Aspheric'};

    % machining paramters
    cutDirection = 'Edge to Center'; % 'Center to Edge'
    spindleDirection = 'Counterclockwise'; % 'Clockwise'
    angularDiscrete = 'Constant Arc'; % 'Constant Angle'
    aimRes = 500;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = 100;
    arcLength = 30;
    maxAngPtDist = 6*pi/180;
    angularLength = 6*pi/180;
end

% plot the importing result
[surfNormIni(:,:,1),surfNormIni(:,:,2),surfNormIni(:,:,3)] = surfnorm( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));

fig1 = figure('Name','original xyz scatters of the surface (sparsely)');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile;
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))));
hold on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);
% nexttile;
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
% hold on;
%     quiver3(surfMesh(1:10:end,1:10:end,1), ...
%         surfMesh(1:10:end,1:10:end,2), ...
%         surfMesh(1:10:end,1:10:end,3), ...
%         surfNormIni(1:10:end,1:10:end,1), ...
%         surfNormIni(1:10:end,1:10:end,2), ...
%         surfNormIni(1:10:end,1:10:end,3), ...
%         'AutoScale','on','Color',[0.85,0.33,0.10], ...
%         'DisplayName','Normal Vectors');
%     legend('Original Points','Orthogonal direction','Location','northeast');
% axis equal;
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel(['x (',unit,')']);
% ylabel(['y (',unit,')']);
% zlabel(['z (',unit,')']);

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
% surfNormSym = [surfFx;surfFy;-1];
% surfNormSym = surfNormSym./norm(surfNormSym);
% surfNormFunc = matlabFunction(surfNormSym,'Vars',{x,y});
% 
% switch spindleDirection
%     case 'Clockwise'
%         conThetaBound = [2*pi,0];
%     case 'Counterclockwise'
%         conThetaBound = [0,2*pi];
% end

% switch angularDiscrete
%     case 'Constant Arc'
%         conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
%             ceil(2*pi/maxAngPtDist) + 1);
%         conTheta(end) = [];
%     case 'Constant Angle'
%         conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
%             ceil(2*pi/angularLength) + 1);
%         conTheta(end) = [];
%     otherwise
%         error('Invalid angular discretization type.');
% end

%% Tool path adjusting
t0 = tic;
tRes0 = tic;
t1 = tic;

% get the f(r) function of the surface
surfFuncr = matlabFunction(subs(surfSym,x,0),'Vars',y);
surfFyr = matlabFunction(subs(surfFy,x,0),'Vars',y);

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius
toolCoefs = toolSp.coefs;

% initialize the cutting pts: the rotation center is seen as the first loop
toolREach = 0; % the radius R of the current loop
toolPathPt = zeros(3,1);
toolPathPt(3,1) = surfFuncr(toolPathPt(2));

uLim = [0;1]; % the interval of each toolpath
peakPt = zeros(3,1);
res = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes

% the first toolpath pt
[toolPathPt(:,1),toolQuat,toolContactU,surfPt] = curvepos( ...
    surfFuncr,surfFyr,toolData,toolPathPt(:,1),[0;0;-1]);
toolNormDirect = quat2rotm(toolQuat)*toolData.toolEdgeNorm;
fprintf('No.1\t toolpath point is calculated.\n-----\n');

% the rest
r = rStep;
ind = 1;
iter = 1;
delta = rStep;
while r <= rMax
    ind = ind + 1;
%     toolPathPt(:,ind) = [0;r;surfFuncr(r)];
    surfPt(:,ind) = [0;r;surfFuncr(r)];
    surfNorm = [0;surfFyr(r);-1];
    surfNorm = surfNorm./norm(surfNorm);
    while iter <= maxIter
        % calculate the surfPt and toolpathPt from center to edge
        [toolPathPt(:,ind),toolQuat(ind,:),toolContactU(ind)] = ...
            curvetippos(toolData,surfPt(:,ind),surfNorm,[0;0;-1]);
%         [toolPathPt(:,ind),toolQuat(ind,:),toolContactU(ind),surfPt(:,ind)] = curvepos( ...
%             surfFuncr,surfFyr,toolData,toolPathPt(:,ind),[0;0;-1]);
        toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

        % calculate the residual height of the loop and the inner nearest loop
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(ind,:))*toolSp1.coefs + toolPathPt(:,ind);
        toolContactPt1 = fnval(toolSp1,toolContactU(ind));
        toolSp2 = toolSp;
        toolSp2.coefs = quat2rotm(toolQuat(ind - 1,:))*toolSp2.coefs + toolPathPt(:,ind - 1);
        toolContactPt2 = fnval(toolSp2,toolContactU(ind - 1));
        [res(ind),peakPt(:,ind),uLim(1,ind),uLim(2,ind)] = residual2D_numeric(toolSp1,toolSp2,1e-3, ...
            toolContactPt1,toolContactPt2,'DSearchn');

        % if residual height does not satisfy the reqiurement
        if res(ind) < aimRes
            break;
        else
            fprintf('\tIter %d: The maximum residual height is %f %cm.\n',iter,max(res(ind)),char([956]));
            delta = delta/3;
            r = r - delta;
            iter = iter + 1;
        end
    end

    if iter == maxIter
        fprintf('The iteration ends with the maxIter is reached. \n');
    end
    delta = rStep;
    r = r + delta;
    iter = 1;
fprintf('No.%d\t toolpath point is calculated.\n-----\n',ind);
end

fprintf('The toolpath concentric optimization process causes %f seconds.\n',toc(tRes0));


%% plot the result above
figure('Name','tool path optimization');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile;
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))),'Color',[0,0.4470,0.7410]);
hold on;
plot(toolPathPt(2,:),toolPathPt(3,:),'.','Color',[0.8500,0.3250,0.0980]);
for jj = 1:size(toolPathPt,2)
    toolSp1 = toolSp;
    toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolSp1.coefs + toolPathPt(:,jj);
    toolSp1Pt = fnval(toolSp1,0:0.01:1);
    plot(toolSp1Pt(2,:),toolSp1Pt(3,:),'Color',[0.6,0.6,0.6],'LineWidth',0.1);
end
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);


nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none'); hold on;
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
plotSpar = 5;
% plot3(toolPathPt(1,1:plotSpar:end), ...
%     toolPathPt(2,1:plotSpar:end), ...
%     toolPathPt(3,1:plotSpar:end), ...
%     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
% if size(surfDomain,1) == 1
%     for jj = 1:ptNum
%         toolSp1 = toolSp;
%         toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolPathPt(:,jj);
%         Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
%         % isQDomain = (Q(1,:).^2 + Q(2,:).^2 - surfDomain(1,2)^2) >= 0;
%         % Q(:,isQDomain) = [];
%         plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1); hold on;
%     end
% else
%     for jj = 1:ptNum
%         toolSp1 = toolSp;
%         toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolPathPt(:,jj);
%         Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
%         plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1); hold on;
%     end
% end
grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
axis equal;
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
% legend('designed surface','tool center point','tool edge','Location','northeast');



% s4_visualize_concentric;

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