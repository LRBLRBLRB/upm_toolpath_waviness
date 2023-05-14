% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

% constant-residual-height based spiral tool path generation, with the
% radius tooltip instead of B-spline, the optimiation of which is to
% directly solve.
% - tool: tooldata
% - path: constant x increment
% no residual & ulim calculation: still have bugs!!!

isAPP = false;
if isAPP
    %% app-used
    workspaceDir = app.workspaceDir;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;

    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);

    % tool data import
    toolData = app.toolData;

    % surface import
    surfFunc = app.surfFuncs;
    surfFx = app.surfFx;
    surfFy = app.surfFy;
    surfDomain = app.surfDomain;
    surfMesh = app.surfMesh;
    rMax = app.rMax;

    % machining paramters
    cutDirection = app.cutDirection;
    startDirection = app.spindleDirection;
    angularIncrement = app.angularDiscrete;
    xIncrement = app.aimRes;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = app.maxIter;
    arcLength = app.arcLength;
    maxAngPtDist = app.maxAngPtDist;
    angularLength = app.angularLength;
else
    %% function-used
    close all;
    clear;
    clc;
    addpath(genpath('funcs'));
    % global variables
    % workspaceDir = fullfile('..','workspace','\20220925-contrast\nagayama_concentric';
    % workspaceDir = fullfile('..','workspace','\20221020-tooltip\tooltip fitting result';
    workspaceDir = uigetdir( ...
        fullfile('..','workspace'), ...
        'select the workspace directory');
    if ~workspaceDir
        workspaceDir = fullfile('..','workspace');
    end
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';

    questOpt.Interpreter = 'tex';
    questOpt.Default = 'OK & Continue';

    diaryFile = fullfile(workspaceDir,['diary',datestr(now,'yyyymmddTHHMMSS'),'.log']);
    diary(diaryFile);
    diary on;

    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);

    % tool data import
    [toolFileName,toolDirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select one tool edge data file', ...
        fullfile(workspaceDir,'tooltheo.mat'), ...
        'MultiSelect','off');
    if ~toolFileName
        fprintf('No tool data file loaded.\n');
        return;
    end
    toolName = fullfile(toolDirName,toolFileName);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % machining paramters
    cutDirection = 'Edge to Center'; % 'Center to Edge'
    startDirection = 'X Minus'; % 'X Minus'
    angularIncrement = 'Constant Arc'; % 'Constant Angle'
    arcLength = 20; % um
    maxAngPtDist = 1*pi/180;
    angularLength = 1*pi/180;
    radialIncrement = 'On-Axis'; % 'Surface'
    xIncrement = 28; % um
    maxIter = 100;
    spiralMethod = 'Radius-Number'; % Radius-Angle
    frMethodDefault = 'Approximation'; % 'Approximation'
    frParamDefault = 1-1e-5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % concentric surface generation / import
    % A = tand(20)/(2*2000);
    c = 0.69/1000/(1000^(aimUnit - presUnit));
    syms x y;
    surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-500,500;-500,500];
    surfDomain = 1.2*surfDomain;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% related parameters
switch startDirection
    case 'X Plus' % plus in this program, but minus in moore
        rMax = max(surfDomain(1,2),surfDomain(2,2));
        rStep = -1*xIncrement;
    case 'X Minus' % minus in this program, but plus in moore
        rMax = min(surfDomain(1,1),surfDomain(2,1)); % reverse
        rStep = 1*xIncrement;
end

switch cutDirection
    case 'Edge to Center'
        cutDirect = [0;-1;0];
        rRange = [rMax,0];
    case 'Center to Edge'
        cutDirect = [0;1;0];
        rRange = [0,rMax];
end

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

% plot the importing result
[surfNormIni(:,:,1),surfNormIni(:,:,2),surfNormIni(:,:,3)] = surfnorm( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));

fig1 = figure('Name','original xyz scatters of the surface (sparsely)');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile;
plot(toolData.toolFit(2,:),toolData.toolFit(3,:),'Color',[0,0.4470,0.7410]);
hold on;
patch('XData',toolData.toolFit(2,:),'YData',toolData.toolFit(3,:),...
    'EdgeColor','none','FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.3);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
title('Tooltip Geometry');
nexttile;
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))));
hold on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);
title('2D-Surface Geometry');

msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s}', ...
    'Surface was generated successfully!\n'],textFontSize,textFontType), ...
    'The workspace directory name is: ', ...
    sprintf('%s\n',getlastfoldername(workspaceDir)), ...
    sprintf('The parameters are listed below:'), ...
    sprintf('1. Surface radius: %f%s',abs(rMax),unit), ...
    sprintf('2. Surface curvature: %f%s^{-1}',c,unit), ...
    sprintf('3. Tool file: %s (radius: %f%s)',toolFileName,toolData.radius,unit), ...
    '4. X increment: ', ...
    sprintf('\tX direction (in program): %s',startDirection), ...
    sprintf('\tX increment: %f%s',xIncrement,unit), ...
    '5. C increment: ', ...
    sprintf('\tIncrement type: %s',angularIncrement), ...
    sprintf('\tArc length: %f%s',arcLength,unit), ...
    sprintf(['\tMax angle: %f',char(176),')\n'],maxAngPtDist*180/pi), ...
    'Ready to continue?'}, ...
    'Surface Generation','OK & Continue','Cancel & quit',questOpt);
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% Tool path adjusting
t0 = tic;
tRes0 = tic;
t1 = tic;

% get the f(r) function of the surface
curveFunc = matlabFunction(subs(surfSym,y,0),'Vars',x);
curveFx = matlabFunction(subs(surfFx,y,0),'Vars',x);

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius

% initialize the cutting pts: the outest loop
curvePt = [rRange(1);0;curveFunc(rRange(1))];
curveNorm = [curveFx(rRange(1));0;-1];
curveNorm = curveNorm./norm(curveNorm);
[curvePathPt,curveQuat,curveContactU] = curvetippos(toolData,curvePt,curveNorm, ...
    [0;0;-1],cutDirect,'directionType','norm-cut');
curveNorm = quat2rotm(curveQuat)*toolData.toolEdgeNorm;
% [curvePathPt,curveQuat,curveContactU,curvePt] = ...
%     curvepos(curveFunc,curveFx,toolData,curvePathPt,[0;0;-1],[0;-1;0]);
% rRange(1) = curvePt(1);
fprintf('-----\nNo.1\t toolpath %f\t[r = %f] is calculated within %fs.\n-----\n',curvePathPt(1,1),curvePt(1,1),toc);

% num of loops
curveNum = ceil(abs(curvePt(1,1))/xIncrement) + 1;
rStep = abs(curvePt(1,1))/(curveNum - 1);
switch startDirection
    case 'X Plus' % plus in this program, but minus in moore
        rStep = -1*rStep;
    case 'X Minus' % minus in this program, but plus in moore
        rStep = 1*rStep;
end

% assign space
curveRes = 5*xIncrement*ones(1,curveNum);
curvePeakPt = zeros(5,curveNum);
curveInterPt = cell(1,curveNum);
curveInterPt{1} = zeros(3,1);
curveULim = cell(1,curveNum);
curveULim{1} = [0;1];
for ind = 2:curveNum
    tic
    curvePt(1,ind) = curvePt(1,ind - 1) + rStep;
    curvePt(3,ind) = curveFunc(curvePt(1,ind));
    curveNorm = [curveFx(curvePt(1,ind));0;-1];
    curveNorm = curveNorm./norm(curveNorm);
    [curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind)] = curvetippos( ...
        toolData,curvePt(:,ind),curveNorm,[0;0;-1],cutDirect, ...
        "directionType",'norm-cut');

    % calculate the residual height of the loop and the inner nearest loop
    toolSp1 = toolSp;
    toolSp1.coefs = quat2rotm(curveQuat(ind,:))*toolSp1.coefs + curvePathPt(:,ind);
    % toolContactPt1 = fnval(toolSp1,curveContactU(ind));
    toolSp2 = toolSp;
    toolSp2.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp2.coefs + curvePathPt(:,ind - 1);
    % toolContactPt2 = fnval(toolSp2,curveContactU(ind - 1));

    [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim{ind}, ...
        curveULim{ind - 1}] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
        curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1});
    fprintf('No.%d\t toolpath %f\t[r = %f] is calculated within %fs.\n-----\n',ind,curvePathPt(1,ind),curvePt(1,ind),toc);
end

if strcmp(startDirection,'X Minus')
    % uLim reverse
    curveULim{1}(1) = 1;
    curveULim{end}(end) = 0;
    for ii = 1:length(curveULim)
        tmp = curveULim{ii};
        tmpSorted = sort(tmp(:),'ascend');
        curveULim{ii} = reshape(tmpSorted,2,[]);
    end
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
for jj = 1:size(curvePathPt,2)
    plot(curveInterPt{jj}(1,:),curveInterPt{jj}(3,:),'.','Color',[0.850,0.3250,0.0980]);
    toolSp1 = toolData.toolBform;
    toolSp1.coefs = quat2rotm(curveQuat(jj,:))*toolSp1.coefs + curvePathPt(:,jj);
    toolSp1Pt = fnval(toolSp1,0:0.0001:curveULim{jj}(1,1));
    toolSp1Pt(3,end) = NaN;
    patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
        'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
        'LineWidth',0.5,'LineStyle','-');
    for ii = 1:size(curveULim{jj},2) - 1
        toolSp1Pt = fnval(toolSp1,curveULim{jj}(1,ii):0.0001:curveULim{jj}(2,ii));
        plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410], ...
            'LineWidth',0.5);
        toolSp1Pt = fnval(toolSp1,curveULim{jj}(2,ii):0.0001:curveULim{jj}(1,ii + 1));
        toolSp1Pt(3,end) = NaN;
        patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
            'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
            'LineWidth',0.5,'LineStyle','-');
    end
    toolSp1Pt = fnval(toolSp1,curveULim{jj}(1,end):0.0001:curveULim{jj}(2,end));
    plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410], ...
        'LineWidth',0.5);
    toolSp1Pt = fnval(toolSp1,curveULim{jj}(2,end):0.0001:1);
    toolSp1Pt(3,end) = NaN;
    patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
        'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
        'LineWidth',0.5,'LineStyle','-');
end
legend('tool path point','tool contact point','ideal surface','peak point','','actual surface');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);

%% concentric toolpath generation for each loop
toolPathAngle = [];
loopPtNum = [];
accumPtNum = 0;
toolNAccum = [];
toolRAccum = [];
toolQuat = [];
toolNormDirect = [];
toolCutDirect = [];
toolContactU = [];
res = [];
uLim = [];
peakPt = [];

if strcmp(startDirection,'X Plus') % 'X Minus'
    conThetaBound = [0,2*pi];
else
    conThetaBound = [2*pi,0];
end

for ii = 1:size(curvePathPt,2)
    conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
        ceil(2*pi/min(maxAngPtDist,arcLength/abs(curvePathPt(1,ii)))) + 1);
    conTheta0(end) = [];
    loopPtNum = [loopPtNum,length(conTheta0)]; % the tool path number of each loop
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)]; % the tool path number of the existing loops
    toolNAccum = [toolNAccum,(ii - 1)*ones(1,loopPtNum(end))]; % the loop No. of each tool path
    toolPathAngle = [toolPathAngle,conTheta0]; % tool path angle
%     if toolPathAngle(accumPtNum(end - 1) + 1) > 6
%         toolPathAngle(accumPtNum(end - 1) + 1) = toolPathAngle(accumPtNum(end - 1) + 1) - 2*pi;
%     end
    toolRAccum = [toolRAccum,curvePt(1,ii)*ones(1,loopPtNum(end))]; % the loop radius of each tool path
    toolContactU = [toolContactU,curveContactU(ii)*ones(1,loopPtNum(end))]; % the B-spline parameter of the contact point
    res = [res,curveRes(ii)*ones(1,loopPtNum(end))]; % the residual height of each tool point and its closest point on the previous loop 
    uLim = [uLim,ndgrid(curveULim(:,ii),1:loopPtNum(end))];
end
accumPtNum(1) = [];
for ii = 1:length(toolPathAngle)
    R = rotz(toolPathAngle(ii)/pi*180);
    kk = find(ii <= accumPtNum,1);
    loopQuat(ii,:) = rotm2quat(R);
    toolPathPt(:,ii) = R*curvePathPt(:,kk); % the concentric tool path point
%     peakPt(ii,:) = R*curvePeakPt(:,kk);
    toolQuat(ii,:) = quatmul(curveQuat(kk,:),loopQuat(ii,:)); 
    toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
end

nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
% colormap('summer');
% cb = colorbar;
% cb.Label.String = ['Height (',unit,')'];
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0.850,0.3250,0.0980]);
grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
axis equal;
legend('designed surface','tool center point','Location','northeast');

warningTone = load('handel');
sound(warningTone.y,warningTone.Fs);
s6_visualize_concentric; % res for two lines

msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
    'Concentric tool path was generated successfully!'],textFontSize,textFontType), ...
    'Ready to continue?'}, ...
    'Concentric tool path Generation','OK & Continue','Cancel & quit',questOpt);
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% Feed rate smoothing
% to smooth the loopR to get the real tool path

% cubic spline approximation
% the function between the numeric label of tool path and surf radius R
toolREach = curvePt(1,:);
diffR = abs(diff(toolREach));
% Fr = csape(accumPtNum,toolREach,[1,1]);
% toolNoTheta = linspace(2*pi*1,2*pi*length(accumPtNum),length(accumPtNum));
% rTheta = csape(toolNoTheta,toolREach,[1,1]);

switch spiralMethod
    case 'Radius-Number'
        SurfEach = accumPtNum;
    case 'Radius-Angle'
        SurfEach = linspace(2*pi*1,2*pi*length(accumPtNum),length(accumPtNum));
end

questOpt1 = questOpt;
questOpt1.Default = 'Refit';

hFeedrate = figure('Name','Feed Rate Smoothing');
% while true
%     approxOut = selectfr(textFontType,textFontSize);
%     Fr = approxOut.fittedmodel;
approxMethod = 'pchip';
switch approxMethod
    case 'csape'
        Fr = csape(SurfEach,toolREach,[1 1]);
    case 'pchip'
        Fr = pchip(SurfEach,toolREach);
end

    yyaxis left;
    scatter(SurfEach,toolREach);
    hold on;
    fnplt(Fr,'r',[SurfEach(1) + 1,SurfEach(end)]);
    plot(1:SurfEach(end),toolRAccum);
    ylim1 = [min(toolREach),max(toolREach)];
    set(gca,'YLim',[2*ylim1(1) - ylim1(2),ylim1(2)]);
    ylabel(['Radius of the Loop (',unit,')']);
    yyaxis right;
    bar(SurfEach(1:end - 1),diffR);
    ylim2 = [min(diffR),max(diffR)];
    set(gca,'YLim',[ylim2(1),2*ylim2(2) - ylim2(1)]);
    ylabel(['Cutting width of each Loop (',unit,')']);
    % line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
    %     'Color',[0.929,0.694,0.1250]);
    set(gca,'FontSize',textFontSize,'FontName',textFontType);
    xlabel('Concentric Angle');
    legend('No.-R scatters','Approx result','Concentric result','Pitch');
    hold off;

%     msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
%         'Feed rate is fittted successfully!'],textFontSize,textFontType), ...
%         'Ready to continue?'}, ...
%         'Feed rate fit','OK & Continue','Refit','Cancel & quit',questOpt1);
%     switch msgfig
%         case {'','Cancel & quit'}
%             return;
%         case 'Refit'
%             continue;
%         case 'OK & Continue'
%             break;
%     end
% end

if exist(fullfile(workspaceDir,['feedrate-',approxMethod,'.fig']),'file')
    savefig(hFeedrate,fullfile(workspaceDir,['feedrate-',approxMethod,'-',datestr(now,'yyyymmddTHHMMSS'),'.fig']));
else
    savefig(hFeedrate,fullfile(workspaceDir,['feedrate-',approxMethod,'.fig']));
end

% nexttile;
% scatter(toolNoTheta,toolREach);
% hold on;
% fnplt(rTheta,'r',[toolPathAngle(accumPtNum(1) + 1),toolPathAngle(end)+2*pi*toolNAccum(end)]);
% plot(toolPathAngle + 2*pi*toolNAccum,toolRAccum);
% % line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
% %     'Color',[0.929,0.694,0.1250]);
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel('Accumulating Toolpath Angle');
% ylabel(['Radius of the Loop (',unit,')']);
% legend('\theta-R scatters','csape result','Concentric result');
% save the feed rate curve
% [smoothFileName,smoothDirName,smoothFileType] = uiputfile({ ...
%     '*.mat','MAT-file(*.mat)'; ...
%     '*.txt','text-file(.txt)';...
%     '*.*','all file(*.*)';...
%     }, ...
%     'Select the directory and filename to save the surface concentric tool path', ...
%     fullfile(workspaceDir,['feedRate',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
% smoothName = fullfile(smoothDirName,smoothFileName);
% switch smoothFileType
%     case 0
%         msgfig = msgbox(sprintf("\\fontsize{%d}\\fontname{%s} No approximation saved"),"Warning","warn","non-modal");
%         uiwait(msgfig);
%     case 1
%         Comments = cell2mat(inputdlg( ...
%             'Enter Comment of the feed rate smoothing processing:', ...
%             'Saving Comments', ...
%             [5 60], ...
%             string(datestr(now))));
%         save(smoothName,"Comments","Fr","rTheta");
% end

%% spiral tool path generation with the smoothing result
% initialize the spiral path
spiralPath0 = zeros(3,accumPtNum(end)); % the spiral tool path & orientation
spiralNorm0 = zeros(3,accumPtNum(end));
spiralCut0 = zeros(3,accumPtNum(end));
spiralQuat0 = zeros(accumPtNum(end),4);
spiralContactU0 = zeros(1,accumPtNum(end));

interpR = fnval(Fr,1:SurfEach(end));
interpR(1:SurfEach(1)) = 0;
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
for kk = 2:numLoop % begin with the 1st loop
    angleN = toolPathAngle(accumPtNumlength(kk - 1) + 1:accumPtNumlength(kk));
    for indInterp = accumPtNum(kk - 1) + 1:accumPtNum(kk)
        %get the inner closest point and linearly interpolate them
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
if strcmp(angularIncrement,'Constant Arc')
    [spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
        spiralAngle0,spiralPath0,spiralContactU0,{spiralNorm0;spiralCut0},spiralQuat0,toolData,'interpType','linear');
%     [spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
%         spiralAngle0,spiralPath0,spiralContactU0,spiralQuat0,{spiralNorm0;spiralCut0},'interpType','linear');
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
figure;
plot3(spiralPath(1,:),spiralPath(2,:),spiralPath(3,:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

figure;
plot3(spiralPath0(1,:),spiralPath0(2,:),spiralPath0(3,:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

tSpiral = toc(tSpiral0);
fprintf('The time spent in the spiral toolpath generation process is %fs.\n',tSpiral);

% msgfig = questdlg({'Spiral tool path was generated successfully!', ...
%     'Ready to continue to simulate?'}, ...
%     'Spiral tool path Generation','OK & Continue','Cancel & quit',questOpt);
% if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
%     return;
% end

spiralFolderName = getlastfoldername(workspaceDir);
[spiralPathFileName,spiralPathDirName,spiralPathFileType] = uiputfile({ ...
    '*.mat','MAT-file(*.mat)'; ...
    '*.txt','text-file(.txt)';...
    '*.*','all file(*.*)';...
    }, ...
    'Select the directory and filename to save the surface tool path', ...
    fullfile(workspaceDir,[spiralFolderName,'-spiralPath-',approxMethod, ...
    '-',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
spiralPathName = fullfile(spiralPathDirName,spiralPathFileName);
save(spiralPathName,"spiralAngle","spiralPath","spiralQuat", ...
    "spiralNorm","spiralCut");

diary off
return;

%% Spiral Residual height calculation of the spiral tool path
spiralRes = 5*xIncrement*ones(2,spiralPtNum);
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
%             toolSp,toolData.toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
        tmpres1 = 5*xIncrement;
        ptres1 = zeros(3,1);
    else
        [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
            spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolSp,toolData.toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
    end

    % outer ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + 2*pi,1,'first');
    ind3 = find(spiralAngle < spiralAngle(ind1) + 2*pi,1,'last');
    if isempty(ind2) || isempty(ind3)
        tmpres2 = 5*xIncrement;
        ptres2 = zeros(3,1);
    else
        [tmpres2,ptres2,spiralULim(:,ind1)] = residual3D( ...
            spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolSp,toolData.toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
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