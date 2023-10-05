% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

isAPP = false;
if isAPP
    questOpt.Interpreter = 'tex';
    questOpt.Default = 'OK & Continue';
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
    param = '';
    surfFunc = app.surfFuncs;
    surfSym = sym(surfFunc);
    surfFx = app.surfFx;
    surfFy = app.surfFy;
    if size(app.surfDomain,1) == 1
        surfDomain = [app.surfDomain;app.surfDomain];
    else
        surfDomain = app.surfDomain;
    end
    surfMesh = app.surfMesh;
    rMax = app.rMax;

    % machining paramters
    cutDirection = app.cutDirection;
    startDirection = app.startDirection;
    angularIncrement = app.angularDiscrete;
    arcLength = app.arcLength;
    maxAngPtDist = app.maxAngPtDist;
    angularLength = app.angularLength;
    radialIncrement = app.radialIncrement; % 'Surface'
    aimRes = app.aimRes;
    rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
    maxIter = app.maxIter;
    spiralMethod = app.spiralMethod;
    zAllowance = app.zAllowance;
    [~,toolFileName,~] = fileparts(app.toolDataFile);
else
    %% function-used
    close all;
    clear;
    clc;
    syms x y;
    questOpt.Interpreter = 'tex';
    questOpt.Default = 'OK & Continue';
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
    unitList = {'m','mm','\mum','nm'};
    
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
    presUnit = find(strcmp(unitList,toolData.unit),1);
    aimUnit = find(strcmp(unitList,unit),1);
    toolData.center = 1000^(aimUnit - presUnit)*toolData.center;
    toolData.radius = 1000^(aimUnit - presUnit)*toolData.radius;
    toolData.toolBform.coefs = 1000^(aimUnit - presUnit)*toolData.toolBform.coefs;
    toolData.toolCpts = 1000^(aimUnit - presUnit)*toolData.toolCpts;
    toolData.toolEdgePt = 1000^(aimUnit - presUnit)*toolData.toolEdgePt;
    toolData.toolFit = 1000^(aimUnit - presUnit)*toolData.toolFit;
    
    diaryFile = fullfile(workspaceDir,['diary',datestr(now,'yyyymmddTHHMMSS'),'.log']);
    diary(diaryFile);
    diary on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % machining paramters
    cutDirection = 'Edge to Center'; % 'Center to Edge'
    startDirection = 'X Minus'; % 'X Minus'
    angularIncrement = 'Constant Arc'; % 'Constant Angle'
    arcLength = 20; % um
    maxAngPtDist = 1*pi/180;
    angularLength = 1*pi/180;
    radialIncrement = 'On-Axis'; % 'Surface'
    aimRes = 0.5; % um
    rStep = toolData.radius/2; % rStep can also be determined by the axial differentiate of the surface
    maxIter = 100;
    spiralMethod = 'Radius-Number'; % Radius-Angle
    frMethodDefault = 'Approximation'; % 'Approximation'
    frParamDefault = 1-1e-5;
    dist2Surf = false;

    % Explanation: 
    % - cutDirection is the direction along which the tool feeds, and it
    %       affects the order of the r-value of the tool path, i.e., rRange
    % - startDirection is the direction where the tool startd to feed. It
    %       determines whether the start point is positive or not. 
    %   Notice that the tool is fixed in the MCS, the startDirection also 
    %       determines the spindle rotation direction. E.g., If the r-value
    %       of the start point is positive, it means that the tool feeds 
    %       from the X+ direction. Therefore, the spindle should rotate in
    %       the counterclockwise direction since the tool rake face is
    %       fixed facing the top. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % concentric surface generation / import
    % A = tand(20)/(2*2000);
    c = 0.69/1000/(1000^(aimUnit - presUnit));
    param = sprintf('c = %f',c);
    syms C;
    surfSymDisp = C.*(x.^2 + y.^2)./(1 + sqrt(1 - C.^2.*(x.^2 + y.^2)));
    surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-500,500;-500,500];
    zAllowance = 1.2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% related parameters
isUIncrease = toolData.toolBform.coefs(end,1) - toolData.toolBform.coefs(1,1);
switch startDirection
    case 'X Plus' % plus both in this program and in moore
        rMax = max(zAllowance*surfDomain(1,2),zAllowance*surfDomain(2,2));
        rStep = -1*rStep;
        % cutDirect = [0;1;0];
    case 'X Minus' % minus both in this program and in moore
        rMax = min(zAllowance*surfDomain(1,1),zAllowance*surfDomain(2,1)); % reverse
        rStep = 1*rStep;
        % cutDirect = [0;-1;0];
end
cutDirect = [0;-1;0]; % aimed cut direction

% uDirection represents the direction of the parameter u which increase
if isUIncrease*rStep > 0
    % ([1,0] & X Plus) or ([0,1] & X Minus)
    uDirection = 'U Plus';
else
    % ([1,0] & X Minus) or ([0,1] & X Plus)
    uDirection = 'U Minus';
end

switch cutDirection
    case 'Edge to Center'
        rRange = [rMax,0];
    case 'Center to Edge'
%             rRange = [0,rMax];
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

symdisp(surfSymDisp);
msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s}', ...
    'Surface was generated successfully!\n'],textFontSize,textFontType), ...
    'The workspace directory name is: ', ...
    sprintf('%s\n',getlastfoldername(workspaceDir)), ...
    sprintf('The parameters are listed below:'), ...
    sprintf('1. Surface radius: %f%s',abs(rMax),unit), ...
    sprintf('2. Surface parameters: %s',param), ...
    sprintf('3. Tool file: %s (radius: %f%s)',toolFileName,toolData.radius,unit), ...
    '4. X increment: ', ...
    sprintf('\tX direction (in program): %s',startDirection), ...
    sprintf('\tAimed residual error: %f%s',aimRes,unit), ...
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
% the first toolpath pt
% [curvePathPt,curveQuat,curveContactU,curvePt] = ...
%     curvepos(curveFunc,curveFx,toolData,curvePathPt,[0;0;-1],[0;-1;0]);
[curvePathPt,curveQuat,curveContactU] = curvetippos(toolData,curvePt, ...
    curveNorm,[0;0;-1],cutDirect,'directionType','norm-cut');
curveNorm = quat2rotm(curveQuat)*toolData.toolEdgeNorm;
% rRange(1) = curvePt(1);

%     scatter(curvePathPt(1,1),curvePathPt(3,1),36,[0.4940,0.1840,0.5560]);
%     toolSp0 = toolData.toolBform;
%     toolSp0.coefs = quat2rotm(curveQuat(1,:))*toolSp0.coefs + curvePathPt(:,1);
%     toolPt0 = fnval(toolSp0,0:0.001:1);
%     toolContactPt0 = fnval(toolSp0,curveContactU(1));
%     plot(toolPt0(1,:),toolPt0(3,:),'Color',[0.7,.7,.7]);
%     scatter(toolContactPt0(1),toolContactPt0(3),18,[0.929,0.694,0.1250],"filled");

t1O = toc(t1);
fprintf('-----\nNo.1\t toolpath %f\t[r = %f] is calculated within %fs.\n-----\n',curvePathPt(1,1),curvePt(1,1),t1O);

% the rest
% opt = optimset('Display','iter','MaxIter',50,'TolFun',1e-6,'TolX',1e-6); % fzero options

% opt = optimoptions(@fsolve,'Display','iter','useParallel',true, ...
%     'MaxIterations',50,'FunctionTolerance',1e-6,'StepTolerance',1e-6); % fsolve options
% % 'PlotFcn',{@optimplotx,@optimplotfval},

opt.XTol = 1e-3;

% opt = optimoptions('ga','UseParallel',false,'Display','iter', ...
%     'MaxGenerations',10,'FitnessLimit',aimRes*0.5);

% opt = optimoptions('particleswarm','UseParallel',true,'Display','iter');

[curvePathPt,curveQuat,curveContactU,curvePt,curveRes,curvePeakPt, ...
    curveInterPt,curveULim] = iterfunc_curvepath_multi_solve(curveFunc,curveFx, ...
    toolData,curvePathPt,curveQuat,curveContactU,curvePt,rStep,aimRes,rRange, ...
    'uDirection',uDirection, ...
    'algorithm','search-bisection','directionType','norm-cut','optimopt',opt,'curveFunc',dist2Surf);

%% the last point
% it should be noticed that the last tooltippt can be a minus value.But if 
% thetoolpathpt is minus, the spiral path would be difficult to generate.
tic;
criterum = curvePt(1,end);
if criterum*curvePathPt(1,end) >= 0
    % which means in the last tooltip lies over the symetrical axis as well
    % as the toolpathpt ( both are <= 0, or both ar >= 0)
    ii = find(curvePt(1,:).*curvePathPt(1,:) < 0,1,'last');
    if ~isempty(ii)
        % spiral pitch is too small that all the curvePt(end)*curvePathPt(end) > 0
        curvePathPt = [curvePathPt,zeros(3,1)];
        % ensure whether to delete the rest point
        questStr = {sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
            'Whether to delete?'],textFontSize,textFontType), ...
            'The last curve path point exceeds the symetrical axis: ', ...
            sprintf('  curvePathPt  curvePt')};
        for tmpii = 1:size(curvePathPt,2) - ii
            questStr{tmpii + 3} = sprintf('  %f  %f',curvePathPt(1,tmpii),curvePt(1,tmpii));
        end
        questStr{tmpii + 4} = sprintf('Ensure to delete the last %d curve path point?', ...
            length(curvePathPt(1,:)) - ii);
        warningTone = load('gong');
        sound(warningTone.y,warningTone.Fs);
        msg = questdlg(questStr,'Whether to delete the exceeding point', ...
            'OK & Continue','Cancel',questOpt);
        waitfor(msg);
        switch msg
            case 'OK'
                curvePathPt(1,ii) = 0;
                % delete the rest point
                curvePathPt(:,(ii + 1):end) = [];
                curveQuat((ii + 1):end,:) = [];
                curveContactU((ii + 1):end) = [];
                curvePt(:,(ii + 1):end) = [];
                curveRes((ii + 1):end) = [];
                curvePeakPt(:,(ii + 1):end) = [];
                curveInterPt(:,(ii + 1):end) = [];
                curveULim((ii + 1):end) = [];
            case 'Cancel'
                curvePathPt(1,end) = 0;
        end
    else
        curvePathPt(1,end) = 0;
    end
else
    % which means in the last tooltip lies over the symetrical axis while
    % the toolpathpt within it
    curvePathPt = [curvePathPt,zeros(3,1)];
    curveULim{size(curvePathPt,2)} = [];
end
ind = size(curvePathPt,2);
[curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind),curvePt(:,ind)] = ...
    curvepos(curveFunc,curveFx,toolData,curvePathPt(:,ind),[0;0;-1],[0;-1;0]);
% calculate the residual height of the loop and the inner nearest loop
toolSp1 = toolSp;
toolSp1.coefs = quat2rotm(curveQuat(ind,:))*toolSp1.coefs + curvePathPt(:,ind);
toolSp2 = toolSp;
toolSp2.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp2.coefs + curvePathPt(:,ind - 1);
prevUlimInd = size(curveULim{ind},2); % get the curveUlim use of the last pathpt
if prevUlimInd > 1
    curveULim{ind - 1}(:,(end - prevUlimInd + 2):end) = [];
end
switch uDirection
    case 'U Plus'
        curveULim{ind - 1}(end) = 1;
        curvePlotSpar = 0.0001;
    case 'U Minus'
        curveULim{ind - 1}(end) = 0;
        curvePlotSpar = -0.0001;
end
if dist2Surf
    [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim{ind}, ...
        curveULim{ind - 1}] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
        curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1}, ...
        'uDirection',uDirection,'aimRes',aimRes,'curveFunc',curveFunc);
else
    [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim{ind}, ...
        curveULim{ind - 1}] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
        curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1}, ...
        'uDirection',uDirection,'aimRes',aimRes);
end
% curvePeakPt(5,ind) = curvePeakPt(5,ind) + length(curveContactU);
fprintf('-----\nNo.%d\t toolpath point at [r = 0] is calculated within %fs.\n-----\n',length(curveContactU),toc);

% get rid of the redundant part of the uLim of the last point
tmpU = 0:abs(curvePlotSpar):1;
toolSp1 = toolSp;
toolSp1.coefs = quat2rotm(curveQuat(end,:))*toolSp1.coefs;
toolSp1Pt = fnval(toolSp1,tmpU);
[~,tmpUInd] = min(abs(toolSp1Pt(2,:)));
innermostU = tmpU(tmpUInd);

if strcmp(startDirection,'X Minus')
    % uLim reverse
    curveULim{1}(1) = 1;
    curveULim{end}(end) = 0;

%     for ii = 1:length(curveULim)
%         tmp = curveULim{ii};
%         tmpSorted = sort(tmp(:),'ascend');
%         curveULim{ii} = reshape(tmpSorted,2,[]);
%     end

    curveULim{end}(end) = innermostU;
end
clear tmpU toolSp1 toolSp2;

fprintf('The toolpath concentric optimization process causes %f seconds.\n',toc(tRes0));

% diary off;

%% plot the result above
figure('Name','tool path optimization');
tiledlayout(2,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile(1,[1,1]);
plot(curvePathPt(1,:),curvePathPt(3,:),'.','Color',[0.4940,0.1840,0.5560]);
hold on;
plot(curvePt(1,:),curvePt(3,:),'.','Color',[0.9290,0.6940,0.1250]);
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))),'Color',[0.7,0.7,0.7],'LineWidth',0.3);
for jj = 1:size(curvePathPt,2)
    plot(curveInterPt{jj}(1,:),curveInterPt{jj}(3,:),'.','Color',[0.850,0.3250,0.0980]);
    toolSp1 = toolData.toolBform;
    toolSp1.coefs = quat2rotm(curveQuat(jj,:))*toolSp1.coefs + curvePathPt(:,jj);
%     toolSp1Pt = fnval(toolSp1,curveULim{jj}(1):curvePlotSpar:curveULim{jj}(2));
%     toolSp1Pt(3,end) = NaN;
%     patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
%         'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
%         'LineWidth',0.5,'LineStyle','-');
    for ii = 1:size(curveULim{jj},2)
        toolSp1Pt = fnval(toolSp1,curveULim{jj}(1,ii):curvePlotSpar:curveULim{jj}(2,ii));
        plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410], ...
            'LineWidth',0.5);
%         toolSp1Pt = fnval(toolSp1,curveULim{jj}(2,ii):0.0001:curveULim{jj}(1,ii + 1));
%         toolSp1Pt(3,end) = NaN;
%         patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
%             'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
%             'LineWidth',0.5,'LineStyle','-');
    end
    toolSp1Pt = fnval(toolSp1,curveULim{jj}(1,end):curvePlotSpar:curveULim{jj}(2,end));
    plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410], ...
        'LineWidth',0.5);
%     toolSp1Pt = fnval(toolSp1,curveULim{jj}(2,end):0.0001:1);
%     toolSp1Pt(3,end) = NaN;
%     patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
%         'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
%         'LineWidth',0.5,'LineStyle','-');
end
legend('tool path point','tool contact point','ideal surface','peak point','','actual surface');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);

nexttile(3);
scatter(curvePeakPt(1,2:end),curveRes(2:end));
xlabel(['r (',unit,')']);
ylabel(['residual error (',unit,')']);
drawnow;

%% concentric toolpath generation for each loop
toolPathAngle = []; % tool path angle
loopPtNum = []; % the tool path number of each loop
accumPtNum = 0; % the tool path number of the existing loops
toolNAccum = []; % the loop No. of each concentric CL points
toolREach = curvePathPt(1,:);
toolRAccum = []; % the loop radius of each concentric CC points
toolQuat = [];
toolNormDirect = [];
toolCutDirect = [];
toolContactU = []; % the B-spline parameter of the contact point
res = []; % the residual height of each tool point and its closest point on the previous loop 
peakPt = [];
peakPlace = [];
uLim = {zeros(2,0)};
interPt = {zeros(3,0)};

if strcmp(startDirection,'X Plus') % 'X Minus'
    conThetaBound = [0,-2*pi];
    uLimOrder2 = 'last';
    uLimOrder3 = 'first';
else
    conThetaBound = [0,2*pi];
    uLimOrder2 = 'first';
    uLimOrder3 = 'last';
end

for ii = 1:size(curvePathPt,2)
    conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
        ceil(2*pi/min(maxAngPtDist,arcLength/abs(curvePathPt(1,ii)))) + 1);
    conTheta0(end) = [];
    loopPtNum = [loopPtNum,length(conTheta0)];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)];
    toolNAccum = [toolNAccum,(ii - 1)*ones(1,loopPtNum(end))];
    toolPathAngle = [toolPathAngle,conTheta0 + ...
        sign(conThetaBound(2) - conThetaBound(1))*2*pi*(ii - 1)];
    toolRAccum = [toolRAccum,curvePathPt(1,ii)*ones(1,loopPtNum(end))];
    toolContactU = [toolContactU,curveContactU(ii)*ones(1,loopPtNum(end))];
    res = [res,curveRes(ii)*ones(1,loopPtNum(end))];
    peakPlace = [peakPlace,curvePeakPt(5,ii)*ones(1,loopPtNum(end))];
end
accumPtNum(1) = [];
% figure;
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.5,'LineStyle','none');
% hold on;
% colormap('summer');
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;
for ii = 1:length(toolPathAngle)
    R = rotz(toolPathAngle(ii)/pi*180);
    kk = find(ii <= accumPtNum,1);
    loopQuat(ii,:) = rotm2quat(R);
    toolPathPt(:,ii) = R*curvePathPt(:,kk); % the concentric tool path point
    toolQuat(ii,:) = quatmul(curveQuat(kk,:),loopQuat(ii,:)); 
    toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
    peakPt(1:3,ii) = R*curvePeakPt(1:3,kk); % the concentric peak point and its state
    peakPt(4,ii) = curvePeakPt(4,kk);
    peakPt(5,ii) = curvePeakPt(5,kk);
    uLim{ii} = curveULim{kk}; % the u limitation of the tool tip of each tool path point
    interPt{ii} = R*curveInterPt{kk}; % the intersection point of each tool path point
    surfPt(:,ii) = R*curvePt(:,kk); % the contact point of each tool path point
%     plot3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
%         '.','MarkerSize',6,'Color',[0.8500 0.3250 0.0980]);
%     q1 = quiver3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
%         toolCutDirect(1,ii),toolCutDirect(2,ii),toolCutDirect(3,ii),100, ...
%         'filled','AutoScale','on','Color',[0.9290 0.6940 0.1250]);
%     q2 = quiver3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
%         toolNormDirect(1,ii),toolNormDirect(2,ii),toolNormDirect(3,ii),100, ...
%         'filled','AutoScale','on','Color',[0.6350 0.0780 0.1840]);
%     drawnow;
%     clear q1 q2;
end


%% 
nexttile(2,[2,1]);
% figure;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
colormap('summer');
% cb = colorbar;
% cb.Label.String = ['Height (',unit,')'];
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
grid on;
axis equal;
% toolCoefs = toolSp.coefs;
% stepNum = abs(log10(abs(curvePlotSpar)));
% waitBar = waitbar(0,'Drawing ...','Name','Concentric Results Drawing', ...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% for ii = 1:accumPtNum(end)
%     % Check for clicked Cancel button
%     if getappdata(waitBar,'canceling')
%         break;
%     end
%     displayData = ii/accumPtNum(end); % Calculate percentage
%     waitbar(displayData,waitBar,['Figure Plotting ... ', ...
%         num2str(roundn(displayData*100,-2),'%.2f'),'%']); % Progress bar dynamic display
%     toolSp1 = toolSp;
%     toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolPathPt(:,ii);
%     for jj = 1:size(uLim{ii},2)
%         uLimRound = round(uLim{ii},stepNum);
%         Q = fnval(toolSp1,uLimRound(1,jj):curvePlotSpar:uLimRound(2,jj));
%         plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
%         drawnow;
%     end
% end
% delete(waitBar);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('designed surface','tool center point','Location','northeast');
drawnow;

warningTone = load('handel');
% sound(warningTone.y,warningTone.Fs);
s6_visualize_concentric_multi;

msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
    'Concentric tool path was generated successfully!'],textFontSize,textFontType), ...
    'Ready to continue?'}, ...
    'Concentric tool path Generation','OK & Continue','Cancel & quit',questOpt);
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% Feed rate smoothing
% to smooth the loopR to get the real tool path

% approximation
% the function between the numeric label of tool path and surf radius R
diffR = abs(diff(toolREach));
% Fr = csape(accumPtNum,toolREach,[1,1]);
% toolNoTheta = linspace(2*pi*1,2*pi*length(accumPtNum),length(accumPtNum));
% rTheta = csape(toolNoTheta,toolREach,[1,1]);

switch spiralMethod
    case 'Radius-Number'
        surfEach = accumPtNum;
        surfAccum = 1:surfEach(end);
    case 'Radius-Angle'
        surfEach = sign(conThetaBound(2) - conThetaBound(1))* ...
            linspace(2*pi,2*pi*(length(accumPtNum)),length(accumPtNum));
        surfAccum = toolPathAngle;
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
        rhoTheta = csape(surfEach,toolREach,[1 1]);
        thetaRho = csape(toolREach,surfEach,[1 1]);
    case 'pchip'
        rhoTheta = pchip(surfEach,toolREach);
        thetaRho = pchip(toolREach,surfEach);
end

    yyaxis left;
    scatter(surfEach,toolREach);
    hold on;
    FrPlot = fnval(rhoTheta,surfAccum);
    plot(surfAccum,FrPlot,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
    plot(surfAccum,toolRAccum,':','Color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
    ylim1 = [min(toolREach),max(toolREach)];
    if strcmp(startDirection,'X Minus')
        set(gca,'YDir','reverse');
        set(gca,'XDir','reverse');
        ylim1 = [abs(max(toolREach)),abs(min(toolREach))];
    end
    set(gca,'YLim',[2*ylim1(1) - ylim1(2),ylim1(2)],'YColor','k','YMinorGrid','on');
    ylabel(['Radius of the Loop (',unit,')']);
    yyaxis right;
    bar(surfEach(1:end - 1),diffR,'EdgeColor','none');
    hold on;
    xlim = get(gca,'XLim');
    line([xlim(1),xlim(2)],[mean(diffR),mean(diffR)],'Color',[0.5,0.5,0.5]);
    ylim2 = [0,max(diffR)];
    set(gca,'YLim',[ylim2(1),2*ylim2(2) - ylim2(1)],'YColor','k','YMinorGrid','on');
    ylabel(['Cutting width of each Loop (',unit,')']);
    % line([0,loopRcsape(end)/(2*pi/maxAngPtDist/rStep)],[0,loopRcsape(end)], ...
    %     'Color',[0.929,0.694,0.1250]);
    set(gca,'FontSize',textFontSize,'FontName',textFontType);
    grid on;
    xlabel('Concentric Angle');
    legend('No.-R scatters','Approx result','Concentric result','Pitch');

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

interpR = fnval(rhoTheta,surfAccum);
% interpR(surfAccum > surfEach(1)) = 0;
accumPtNumlength = [0,accumPtNum];

numLoop = length(accumPtNum);

% video capture & debug
% figure('Name','concentric-to-spiral debug & video');
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.3,'LineStyle','none'); hold on;
% colormap('summer');
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel('x'); ylabel('y'); zlabel('z');
% % view(0,90);
% axis equal;
% plot3(toolPathPt(1,:),toolPathPt(2,:),toolPathPt(3,:), ...
%     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);

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
        % get the inner closest point and linearly interpolate them
        angleDel = angleN - (toolPathAngle(indInterp) - sign(conThetaBound(2) - conThetaBound(1))*2*pi);
        [ind1,ind2] = getInnerLoopToolPathIndex(angleN,angleDel); % get the closest tol point in the inner loop
        ind1 = ind1 + accumPtNumlength(kk - 1); % get the index of the closest in the whole list
        ind2 = ind2 + accumPtNumlength(kk - 1);
        [spiralPath0(:,indInterp),spiralContactU0(indInterp),spiralQuat0(indInterp,:)] = toolInterp( ...
            interpR(indInterp),toolRAccum,indInterp,ind1,ind2,toolPathPt,toolContactU,toolQuat,toolCutDirect(:,indInterp));
        spiralNorm0(:,indInterp) = quat2rotm(spiralQuat0(indInterp,:))*toolData.toolEdgeNorm;
        spiralCut0(:,indInterp) = quat2rotm(spiralQuat0(indInterp,:))*toolData.cutDirect;
        % test & debug & video
        % plot3(spiralPath0(1,indInterp),spiralPath0(2,indInterp),spiralPath0(3,indInterp), ...
        %     '.','MarkerSize',6,'Color',[0.8500 0.3250 0.0980]);
        % q1 = quiver3(spiralPath0(1,indInterp),spiralPath0(2,indInterp),spiralPath0(3,indInterp), ...
        %     spiralCut0(1,indInterp),spiralCut0(2,indInterp),spiralCut0(3,indInterp),100, ...
        %     'filled','AutoScale','on','Color',[0.9290 0.6940 0.1250]);
        % q2 = quiver3(spiralPath0(1,indInterp),spiralPath0(2,indInterp),spiralPath0(3,indInterp), ...
        %     spiralNorm0(1,indInterp),spiralNorm0(2,indInterp),spiralNorm0(3,indInterp),100, ...
        %     'filled','AutoScale','on','Color',[0.6350 0.0780 0.1840]);
        % drawnow;
        % % toolSp1 = toolSp;
        % % toolSp1.coefs = quat2rotm(spiralQuat0(indInterp,:))*toolSp.coefs + toolVec0(:,indInterp);
        % % Q = fnval(toolSp1,uLim(1,indInterp):0.01:uLim(2,indInterp));
        % % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
        % clear q1 q2;
    end
end

% elliminate the first loop, which focuses on the z=0 point
spiralPath0(:,1:accumPtNum(1)) = [];
spiralNorm0(:,1:accumPtNum(1)) = [];
spiralCut0(:,1:accumPtNum(1)) = [];
spiralQuat0(1:accumPtNum(1),:) = [];
spiralContactU0(1:accumPtNum(1)) = [];

% change the spiral tool path to constant arc length one
spiralAngle0 = toolPathAngle; % + sign(conThetaBound(1) - conThetaBound(2))*2*pi*toolNAccum;
spiralAngle0(:,1:accumPtNum(1)) = [];
if strcmp(angularIncrement,'Constant Arc')
    [spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
        spiralAngle0,spiralPath0,spiralContactU0,toolData,spiralQuat0,{spiralNorm0;spiralCut0},'interpType','spline');
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
for ii = 1:length(accumPtNum)
    plot3(toolPathPt(1,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        toolPathPt(2,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        toolPathPt(3,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',0.1);
end
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');

figure;
plot3(spiralPath0(1,:),spiralPath0(2,:),spiralPath0(3,:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
for ii = 1:length(accumPtNum)
    plot3(toolPathPt(1,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        toolPathPt(2,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        toolPathPt(3,accumPtNumlength(ii) + 1:accumPtNumlength(ii + 1)), ...
        'Color',[0.8500 0.3250 0.0980],'LineStyle','--','LineWidth',0.1);
end
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');

tSpiral = toc(tSpiral0);
fprintf('The time spent in the spiral toolpath generation process is %fs.\n',tSpiral);

spiralFolderName = getlastfoldername(workspaceDir);
[spiralPathFileName,spiralPathDirName,spiralPathFileType] = uiputfile({ ...
    '*.mat','MAT-file(*.mat)'; ...
    '*.txt','text-file(.txt)';...
    '*.*','all file(*.*)';...
    }, ...
    'Select the directory and filename to save the surface tool path', ...
    fullfile(workspaceDir,append(spiralFolderName(1),'-spiralPath-',approxMethod, ...
    '-',datestr(now,'yyyymmddTHHMMSS'),'.mat')));
spiralPathName = fullfile(spiralPathDirName,spiralPathFileName);
if spiralPathFileName
    save(spiralPathName,"spiralAngle","spiralPath","spiralQuat", ...
        "spiralNorm","spiralCut");
end

msgfig = questdlg({'Spiral tool path was generated successfully!', ...
    'Ready to continue to simulate?'}, ...
    'Spiral tool path Generation','OK & Continue','Cancel & quit',questOpt);
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% Spiral Residual height calculation of the spiral tool path
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

% figure;
% surf( ...
%     surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');
% hold on;
% waitBar = waitbar(0,'Figure Plotting ...','Name','Residual Results Plot', ...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(waitBar,'canceling',0);
% if dist2Surf
%     curveFunc3D = curveFunc;
% else
%     curveFunc3D = [];
% end
curveFunc3D = [];
for ind1 = 1:spiralPtNum
    % find the outer side of ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) - conThetaBound(end),1,uLimOrder2);
    ind3 = find(spiralAngle < spiralAngle(ind1) - conThetaBound(end),1,uLimOrder3);
    if isempty(ind2) || isempty(ind3)
%         ind2 = find(spiralAngle >= spiralAngle(ind1) + pi,1,'first');
%         ind3 = find(spiralAngle < spiralAngle(ind1) + pi,1,'last');
%         [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
%             spiralPath,spiralNorm,spiralCut,spiralContactU, ...
%             toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
        spiralULim{ind1} = uLimIni;
        tmpRes1 = 5*aimRes;
        tmpPeak1 = zeros(5,1);
    else
        [tmpRes1,tmpPeak1,spiralInterPtIn{ind1},spiralULim{ind1}] = ...
            residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
            toolData,toolRadius,spiralULim{ind1},aimRes,uLimIni,curveFunc3D,ind1,ind2,ind3);
    end

    % find the inner side of ulim & residual height
    ind2 = find(spiralAngle >= spiralAngle(ind1) + conThetaBound(end),1,uLimOrder2);
    ind3 = find(spiralAngle < spiralAngle(ind1) + conThetaBound(end),1,uLimOrder3);
    if isempty(ind2) || isempty(ind3)
        tmpRes2 = 5*aimRes;
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
            toolData,toolRadius,spiralULim{ind1},aimRes,uLimIni,curveFunc3D,ind1,ind2,ind3);
    end
    
    spiralRes(:,ind1) = [tmpRes1;tmpRes2];
    spiralPeakPt(:,ind1) = [tmpPeak1;tmpPeak2];

    % debug
%     h = scatter3(spiralPath(1,ind1),spiralPath(2,ind1),spiralPath(3,ind1),12); %'MarkerEdgeColor',[0,0.4450,0.7410]);
%     toolSp1 = toolData.toolBform;
%     toolSp1.coefs = quat2rotm(spiralQuat(ind1,:))*toolSp1.coefs + spiralPath(:,ind1);
%     for ii = 1:size(spiralULim{ind1},2)
%         tmp = fnval(toolSp1,spiralULim{ind1}(1,ii):curvePlotSpar*100:spiralULim{ind1}(2,ii));
%         plot3(tmp(1,:),tmp(2,:),tmp(3,:),'.','Color',h.CData);
%     end
%     displayData = ind1/spiralPtNum; % Calculate percentage
%     waitbar(displayData,waitBar,['Figure Plotting ... ', ...
%         num2str(roundn(displayData*100,-2),'%.3f'),'%']); % Progress bar dynamic display
%     if getappdata(waitBar,'canceling'), break; end
disp(ind1);
end

% get rid of the redundant part of the uLim of the innermost circle
rdomain = abs(surfDomain(1,2)/zAllowance);
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
diary off;
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));