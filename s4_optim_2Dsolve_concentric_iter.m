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
    % save('input_data/surface/ellipsoidAray.mat', ...
    %    "surfMesh","surfNorm","surfCenter");
    
    surfType = '2D';
    Geometry2DCell = {'Rotating Paraboloid','Aspheric'};

    % machining paramters
    cutDirection = 'Center to Edge'; % 'Center to Edge'
    spindleDirection = 'Counterclockwise'; % 'Clockwise'
    angularDiscrete = 'Constant Arc'; % 'Constant Angle'
    aimRes = 0.5;
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

%% Tool path adjusting
t0 = tic;
tRes0 = tic;
t1 = tic;

% get the f(r) function of the surface
surfFuncr = matlabFunction(subs(surfSym,y,0),'Vars',x);
surfFxr = matlabFunction(subs(surfFy,y,0),'Vars',x);

switch cutDirection
    case 'Edge to Center'
        rRange = [rMax,0];
        conThetaBound = [2*pi,0];
        rStep = -1*rStep;
        % initialize the cutting pts: the outest loop
        surfPt = [rRange(1);0;surfFuncr(rRange(1))];
        surfNorm = [surfFxr(rRange(1));0;-1];
        surfNorm = surfNorm./norm(surfNorm);
        % the first toolpath pt
        [curvePathPt,curveQuat,curveContactU] = curvetippos(toolData,surfPt,surfNorm, ...
            [0;0;-1],[0;-1;0],'directionType','norm-cut');
        toolNormDirect = quat2rotm(curveQuat)*toolData.toolEdgeNorm;
        fprintf('No.1\t toolpath point is calculated.\n-----\n');
    case 'Center to Edge'
        rRange = [0,rMax];
        conThetaBound = [0,2*pi];
        % initialize the cutting pts: the rotation center is seen as the first loop
        curvePathPt = [0;rRange(1);surfFuncr(rRange(1));];
        
        % the first toolpath pt
        [curvePathPt(:,1),curveQuat,curveContactU,surfPt] = curvepos( ...
            surfFuncr,surfFxr,toolData,curvePathPt(:,1),[0;0;-1],[0;1;0], ...
            'useParallel',true,'finalDisplay','iter-detailed');
        toolNormDirect = quat2rotm(curveQuat)*toolData.toolEdgeNorm;
        fprintf('No.1\t toolpath point is calculated.\n-----\n');
    otherwise
        errordlg('Invalid cut-direction value!','Parameter Error');
end
        
% the rest
[curvePathPt,curveQuat,curveContactU,surfPt,curveRes,curvePeakPt,curveULim] = ...
    iterfunc_curvepath_solve(surfFuncr,surfFxr,toolData,curvePathPt, ...
    curveQuat,curveContactU,surfPt,rStep,aimRes,rRange,"useParallel",true);

fprintf('The toolpath concentric optimization process causes %f seconds.\n',toc(tRes0));

%% plot the result above
figure('Name','tool path optimization');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile;
plot(curvePathPt(1,:),curvePathPt(3,:),'.','Color',[0.8500,0.3250,0.0980]);
hold on;
plot(surfPt(1,:),surfPt(3,:),'.','Color',[0.9290,0.6940,0.1250]);
plot(peakPt(1,:),peakPt(3,:),'.','Color',[0.850,0.3250,0.0980]);
rSpar = linspace(0,rMax,spar);
plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))),'Color',[0.7,0.7,0.7],'LineWidth',0.3);
for jj = 1:size(curvePathPt,2)
    toolSp1 = toolData.toolBform;
    toolSp1.coefs = quat2rotm(curveQuat(jj,:))*toolSp1.coefs + curvePathPt(:,jj);
    toolSp1Pt = fnval(toolSp1,curveULim(1,jj):0.0001:curveULim(2,jj));
    plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410],'LineWidth',0.5);
end
legend('tool path point','tool contact point','peak point','ideal surface','actual surface');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['r (',unit,')']);
ylabel(['z (',unit,')']);

nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.6,'LineStyle','none');
hold on;
% concentric toolpath for each loop
spiralAngle0 = [];
loopPtNum = [];
accumPtNum = 0;
res = [];
uLim = [];
for ii = 1:size(curvePathPt,2)
    conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
        ceil(2*pi/min(maxAngPtDist,arcLength/curvePathPt(2,ii))) + 1);
    conTheta0(end) = [];
    spiralAngle0 = [spiralAngle0,conTheta0];
    loopPtNum = [loopPtNum,length(conTheta0)];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)];
    res = [res,curveRes(ii)*ones(1,loopPtNum(end))];
    uLim = [uLim,ndgrid(curveULim(:,ii),1:loopPtNum(end))];
end
accumPtNum(1) = [];
for ii = 1:length(spiralAngle0)
    R = rotz(spiralAngle0(ii));
    kk = find(ii <= accumPtNum,1);
    loopQuat(ii,:) = rotm2quat(R);
    toolPathPt(:,ii) = R*curvePathPt(:,kk);
    peakPt(ii,:) = R*curvePeakPt(:,kk);
    toolQuat(ii,:) = quatmul(curveQuat(kk,:),loopQuat(ii,:));
end
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

s4_visualize_concentric; % res for two lines

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
toolNoTheta = 0:2*pi:2*pi*size(toolPathPt,2);
rTheta = csape(toolNoTheta,toolPathPt(2,:),[1,1]);

figure('Name','Feed Rate Smoothing');
scatter(toolNoTheta,toolPathPt(2,:));
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

% switch spindleDirection
%     case 'Clockwise'
%         conThetaBound = [2*pi,0];
%     case 'Counterclockwise'
%         conThetaBound = [0,2*pi];
% end





% initialize the spiral path
spiralPath = zeros(3,accumPtNum(end)); % the spiral tool path
spiralNorm = zeros(3,accumPtNum(end));
spiralCutDir = zeros(3,accumPtNum(end));
spiralPath(:,1:accumPtNum(1)) = toolPathPt(:,1:accumPtNum(1));
spiralNorm(:,1:accumPtNum(1)) = toolNormDirect(:,1:accumPtNum(1));
spiralCutDir(:,1:accumPtNum(1)) = toolCutDirect(:,1:accumPtNum(1));

% the concentric angle of each tool path
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));

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

switch angularDiscrete
    case 'Constant Arc'
        [spiralAngle,spiralPath,spiralContactU,spiralQuat,spiralVec] = arclengthparam(arcLength,maxAngPtDist, ...
            spiralAngle0,spiralPath0,spiralContactU0,{spiralNorm0;spiralCut0},toolData,'interpType','linear'); % spiralQuat0,
        spiralNorm = spiralVec{1};
        spiralCut = spiralVec{2};
    case 'Constant Angle'
        spiralAngle = spiralAngle0;
        spiralPath = spiralPath0;
        spiralNorm = spiralNorm0;
        spiralCut = spiralCut0;
        spiralQuat = spiralQuat0;
        spiralContactU = spiralContactU0;
    otherwise
        error('Invalid angular discretization type.');
end







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


%%
% delete(parObj);
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));