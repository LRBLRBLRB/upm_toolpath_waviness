% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

% 实际上，我打算把concentric和freeform的方案给合并。目前aspheric文件中的是旧的刀位点计算方案，freeform中是新的

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

workspaceDir = 'workspace\20220925-contrast\nagayama_concentric';
diaryFile = fullfile('workspace\diary',['diary',datestr(now,'yyyymmddTHHMMSS')]);
% diary diaryFile;
% diary on;

% profile on
% tPar0 = tic;
% parObj = gcp;
% tPar = toc(tPar0);
% fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);

radius = 190;

% surface data import
A = tand(20)/(2*2000);
C = 0;
syms x y;
surfSym = A*(x.^2 + y.^2) + C;
surfFunc = matlabFunction(surfSym);
surfFx = diff(surfFunc,x);
surfFy = diff(surfFunc,y);
surfDomain = [-2000,2000;-2000,2000];
rMax = max(surfDomain(1,2),surfDomain(2,2));
% sampling density
spar = 501;
conR = linspace(0,rMax,spar); % concentric radius vector
conTheta = linspace(0,2*pi,spar);
[rMesh,thetaMesh] = meshgrid(conR,conTheta);
surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));

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
% % quiver3(surfMesh(1:10:end,1:10:end,1), ...
% %     surfMesh(1:10:end,1:10:end,2), ...
% %     surfMesh(1:10:end,1:10:end,3), ...
% %     surfNormIni(1:10:end,1:10:end,1), ...
% %     surfNormIni(1:10:end,1:10:end,2), ...
% %     surfNormIni(1:10:end,1:10:end,3), ...
% %     'AutoScale','on','Color',[0.85,0.33,0.10], ...
% %     'DisplayName','Normal Vectors');
% % legend('Original Points','Orthogonal direction','Location','northeast');
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

%% machining paramters
cutDirection = 'Edge to Center'; % 'Center to Edge'
spindleDirection = 'Counterclockwise'; % 'Clockwise'
angularDiscrete = 'Constant Arc'; % 'Constant Angle'
conThetaBound = [2*pi,0];
aimRes = 0.2;
rStep = radius/2; % 每步步长可通过曲面轴向偏导数确定
maxIter = 100;
arcLength = 30;
maxAngPtDist = 6*pi/180;
angularLength = 6*pi/180;

%% Tool path adjusting
t0 = tic;
tRes0 = tic;
t1 = tic;

% get the f(r) function of the surface
curveFunc = matlabFunction(subs(surfSym,y,0),'Vars',x);
curveFx = matlabFunction(subs(surfFx,y,0),'Vars',x);

rRange = [rMax,0];
rStep = -1*rStep;

% initialize the cutting pts: the outest loop
curvePt = [rRange(1);0;curveFunc(rRange(1))];
curveNorm = [curveFx(rRange(1));0;-1];
curveNorm = curveNorm./norm(curveNorm);
% the first toolpath pt
% ---!!! the tool should be fixed just as the picture !!!---
[curvePathPt,curveQuat,curveContactU] = radiustippos(radius,curvePt,curveNorm, ...
    [0;0;-1],[0;-1;0],'directionType','norm-feed');
curvePathNorm = quat2rotm(curveQuat)*[0;0;-1];
fprintf('No.1\t toolpath point is calculated.\n-----\n');

%     scatter(curvePathPt(1,1),curvePathPt(3,1),36,[0.4940,0.1840,0.5560]);
%     toolThe0 = 0:0.01:2*pi;
%     toolPt0(1,:) = curvePathPt(1,1) + radius*cos(toolThe0);
%     toolPt0(3,:) = curvePathPt(3,1) + radius*sin(toolThe0);
%     plot(toolPt0(1,:),toolPt0(3,:),'Color',[0.7,.7,.7]);
%     axis equal

% the rest
% opt = optimset('Display','iter','MaxIter',50,'TolFun',1e-6,'TolX',1e-6); % fzero options

% opt = optimoptions(@fsolve,'Display','iter','useParallel',true, ...
%     'MaxIterations',50,'FunctionTolerance',1e-6,'StepTolerance',1e-6); % fsolve options
% % 'PlotFcn',{@optimplotx,@optimplotfval},

opt = optimoptions('ga','UseParallel',true,'Display','diagnose');

% opt = optimoptions('particleswarm','UseParallel',true,'Display','iter');

[curvePathPt,curveQuat,curveContactU,curvePt,curveRes,curveULim] = ...
    iterfunc_curvepath_radius_res(curveFunc,curveFx,radius,curvePathPt,curveQuat, ...
    curveContactU,curvePt,rStep,aimRes,rRange,'algorithm','genetic','optimopt',opt);

% calculate the center point
if curvePathPt(1,end) < 0
    tic;
    curvePathPt(1,end) = 0;
    [curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind),curvePt(:,ind)] = ...
    radiuspos(curveFunc,curveFx,radius,curvePathPt(:,ind),[0;0;-1],[0;-1;0]);

    % calculate the residual height of the loop and the inner nearest loop
    a = norm(curvePathPt(:,ind) - curvePathPt(:,ind - 1));
    curveRes(ind) = radius - sqrt(radius^2 - a^2/4);
    Norm = norm(curveResPathPt(:,ind) - curveResPathPt(:,ind - 1));
    peakPt = 1/2*(curveResPathPt(:,ind) + curveResPathPt(:,ind - 1)) + sqrt(radius^2 - Norm^2/4) ...
        *roty(-90,'deg')*(curveResPathPt(:,ind) - curveResPathPt(:,ind - 1))/Norm;    
    vec1 = curveResPathPt(:,ind) - peakPt;
    curveRes(1,ind) = atan2(vec1(end),vec1(1));
    vec2 = curveResPathPt(:,ind - 1) - peakPt;
    curveRes(2,ind - 1) = atan2(vec2(end),vec2(1));
    fprintf('-----\nNo.%d\t toolpath point at [r = 0] is calculated within %fs.\n-----\n',ind,toc);
end

fprintf('The toolpath concentric optimization process causes %f seconds.\n',toc(tRes0));

% diary off;

%% plot the result above
figure('Name','tool path optimization');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile;
plot(curvePathPt(1,:),curvePathPt(3,:),'.','Color',[0.4940,0.1840,0.5560]);
hold on;
plot(curvePt(1,:),curvePt(3,:),'.','Color',[0.9290,0.6940,0.1250]);
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))),'Color',[0.7,0.7,0.7],'LineWidth',0.3);
theta = 0:0.01:2*pi;
theta = [theta,theta(1)];
for jj = 1:size(curvePathPt,2)
    patch('XData',curvePathPt(1,jj) + radius*cos(theta), ...
        'YData',[curvePathPt(3,jj) + radius*sin(theta(1:end - 1)),NaN], ...
    'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
    'LineWidth',0.5,'LineStyle','-');
%     plot(curvePathPt(1,jj) + radius*cos(theta),curvePathPt(3,jj) + radius*sin(theta), ...
%         'Color',[0,0.4470,0.7410]);
end
legend('tool path point','tool contact point','ideal surface','actual surface');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);

nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;

% concentric toolpath for each loop
toolPathAngle = [];
loopPtNum = [];
accumPtNum = 0;
toolNAccum = [];
toolQuat = [];
toolNormDirect = [];
toolCutDirect = [];
res = [];
uLim = [];
peakPt = [];
for ii = 1:size(curvePathPt,2) - 1
    conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
        ceil(2*pi/min(maxAngPtDist,arcLength/curvePathPt(1,ii))) + 1);
    conTheta0(end) = [];
    toolPathAngle = [toolPathAngle,conTheta0];
    loopPtNum = [loopPtNum,length(conTheta0)];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)];
    toolNAccum = [toolNAccum,ii*ones(1,loopPtNum(end))];
    res = [res,curveRes(ii)*ones(1,loopPtNum(end))];
    uLim = [uLim,ndgrid(curveULim(:,ii),1:loopPtNum(end))];
end
conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
    ceil(2*pi/maxAngPtDist) + 1);
conTheta0(end) = [];
toolPathAngle = [toolPathAngle,conTheta0];
loopPtNum = [loopPtNum,length(conTheta0)];
accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)];
toolNAccum = [toolNAccum,(ii + 1)*ones(1,loopPtNum(end))];
res = [res,curveRes(ii + 1)*ones(1,loopPtNum(end))];
uLim = [uLim,ndgrid(curveULim(:,ii + 1),1:loopPtNum(end))];
accumPtNum(1) = [];
for ii = 1:length(toolPathAngle)
    R = rotz(toolPathAngle(ii));
    kk = find(ii <= accumPtNum,1);
    loopQuat(ii,:) = rotm2quat(R);
    toolPathPt(:,ii) = R*curvePathPt(:,kk);
    peakPt(ii,:) = R*curvePeakPt(:,kk);
    toolQuat(ii,:) = quatmul(curveQuat(kk,:),loopQuat(ii,:));
%     toolNormDirect(ii,:) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
%     toolCutDirect(ii,:) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
end
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0.850,0.3250,0.0980]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% spiralAngle0 = toolPathAngle;

s6_visualize_concentric; % res for two lines

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
toolNoTheta = 0:2*pi:2*pi*(size(curvePathPt,2) - 1);
rTheta = csape(toolNoTheta,curvePathPt(2,:),[1,1]);

figure('Name','Feed Rate Smoothing');
scatter(toolNoTheta,curvePathPt(2,:));
hold on;
fnplt(rTheta,'r',[toolNoTheta(1),toolNoTheta(end)]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel('Accumulating Toolpath Angle');
ylabel(['Radius of the Loop (',unit,')']);
legend('\theta-R scatters','csape result');
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
        save(smoothName,"Comments","rTheta");
end

%% spiral tool path generation with the smoothing result

spiralPath0 = zeros(3,length(toolPathAngle)); % the spiral tool path
spiralQuat0 = zeros(length(toolPathAngle),4);
spiralContactU0 = zeros(1,length(toolPathAngle));

tSpiral0 = tic;

% for each loop, shift the tool path point by decreasing the radius
for kk = 2:size(curvePathPt,2) % begin with the 2nd loop
    % the concentric tool path of the current loop
    angleN = toolPathAngle(accumPtNumlength(kk - 1) + 1:accumPtNumlength(kk));
    for indInterp = accumPtNum(kk - 1) + 1:accumPtNum(kk)
        angleDel = angleN - toolPathAngle(indInterp);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel); % get the closest tol point in the inner loop
        ind1 = ind1 + accumPtNumlength(kk - 1); % get the index of the closest in the whole list
        ind2 = ind2 + accumPtNumlength(kk - 1);
        [spiralPath0(:,indInterp),spiralContactU0(indInterp),spiralQuat0(indInterp,:)] = toolInterp( ...
            interpR(indInterp),toolRAccum,indInterp,ind1,ind2,curvePathPt,curveContactU,curveQuat,toolCutDirect(:,indInterp));
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
toolPathAngle = toolPathAngle + 2*pi*toolNAccum;
toolPathAngle(:,1:accumPtNum(1) - 1) = [];
[spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
    toolPathAngle,spiralPath0,spiralContactU0,{spiralNorm0;spiralCut0},toolData,'interpType','linear'); % spiralQuat0,
spiralNorm = spiralVec{1};
spiralCut = spiralVec{2};

spiralPtNum = length(spiralAngle);

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
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));