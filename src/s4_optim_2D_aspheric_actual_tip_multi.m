% optimization of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path adjusting for the first loop
% Step three: tool path adjusting for the rest
% Step four: simulation of the machining surface
% Step Five: generate the actual toolpath

close all;
clear;
clc;
addpath(genpath('funcs'));
% global variables
% global textFontSize textFontType;
% workspaceDir = fullfile('..','workspace','\20220925-contrast\nagayama_concentric';
workspaceDir = fullfile('..','workspace','\20221020-tooltip\tooltip fitting result');
unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

msgOpts.Default = 'Cancel and quit';
msgOpts.Interpreter = 'tex';

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

% surface data import
A = tand(20)/(2*2000);
C = 0;
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

% machining paramters
cutDirection = 'Edge to Center'; % 'Center to Edge'
spindleDirection = 'Counterclockwise'; % 'Clockwise'
angularDiscrete = 'Constant Arc'; % 'Constant Angle'
aimRes = 0.2;
rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
maxIter = 100;
arcLength = 30;
maxAngPtDist = 6*pi/180;
angularLength = 6*pi/180;

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
nexttile;
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
%     legend('Original Points','Orthogonal direction','Location','northeast');
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

%% Tool path adjusting
t0 = tic;
tRes0 = tic;
t1 = tic;

% get the f(r) function of the surface
curveFunc = matlabFunction(subs(surfSym,y,0),'Vars',x);
curveFx = matlabFunction(subs(surfFx,y,0),'Vars',x);

conThetaBound = [2*pi,0];
rRange = [rMax,0];
rStep = -1*rStep;

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius
toolCoefs = toolSp.coefs;

if rRange(2) > rRange(1) % to ensure the rake face on top
    % if range(2) > rRange(1), then feed direction is center-to-edge, and
    % the toolDirect is [0;1;0]. So cut direction is [-1;0;0]
    cutDirect = [0;1;0];
else
    cutDirect = [0;-1;0];
end

% initialize the cutting pts: the outest loop
curvePt = [rRange(1);0;curveFunc(rRange(1))];
curveNorm = [curveFx(rRange(1));0;-1];
curveNorm = curveNorm./norm(curveNorm);
% the first toolpath pt
[curvePathPt,curveQuat,curveContactU] = curvetippos(toolData,curvePt,curveNorm, ...
    [0;0;-1],[0;-1;0],'directionType','norm-cut');
curveNorm = quat2rotm(curveQuat)*toolData.toolEdgeNorm;

%     scatter(curvePathPt(1,1),curvePathPt(3,1),36,[0.4940,0.1840,0.5560]);
%     toolSp0 = toolData.toolBform;
%     toolSp0.coefs = quat2rotm(curveQuat(1,:))*toolSp0.coefs + curvePathPt(:,1);
%     toolPt0 = fnval(toolSp0,0:0.001:1);
%     toolContactPt0 = fnval(toolSp0,curveContactU(1));
%     plot(toolPt0(1,:),toolPt0(3,:),'Color',[0.7,.7,.7]);
%     scatter(toolContactPt0(1),toolContactPt0(3),18,[0.929,0.694,0.1250],"filled");

fprintf('No.1\t toolpath point is calculated.\n-----\n');

% the rest
curveULim = {[0;1]}; % the interval of each toolpath
curvePeakPt = zeros(5,1);
curveInterPt = {zeros(3,1)};
curveRes = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes
r = rRange(1);
ind = 1;
while (r - rRange(2))*rStep < 0
    delta = 2*sqrt((2*toolData.radius*aimRes - aimRes^2))*cos(atan(curveFx(r)));
    % direction to iterate
    delta = sign(rStep)*delta;
    r = r + delta;
    ind = ind + 1;
    curveULim{ind} = [0;1];

    curvePt(:,ind) = [r;0;curveFunc(r)];
    surfNorm = [curveFx(r);0;-1];
    surfNorm = surfNorm./norm(surfNorm);
    % calculate the surfPt and toolpathPt from center to edge
    [curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind)] = ...
        curvetippos(toolData,curvePt(:,ind),surfNorm,[0;0;-1],cutDirect, ...
        "directionType",'norm-cut');
    % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

    % calculate the residual height of the loop and the inner nearest loop
    toolSp1 = toolSp;
    toolSp1.coefs = quat2rotm(curveQuat(ind,:))*toolSp1.coefs + curvePathPt(:,ind);
    toolContactPt1 = fnval(toolSp1,curveContactU(ind));
    toolSp2 = toolSp;
    toolSp2.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp2.coefs + curvePathPt(:,ind - 1);
    toolContactPt2 = fnval(toolSp2,curveContactU(ind - 1));

    [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim{ind}, ...
        curveULim{ind - 1}] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
        curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1});
    curvePeakPt(5,ind) = curvePeakPt(5,ind) + ind;

    if isinf(curveRes(ind))
        fprintf('No intersection of the current tool path.\n');
        return;
    end

%         scatter(curvePathPt(1,ind),curvePathPt(3,ind),36,[0.4940,0.1840,0.5560]);
%         scatter(curvePathPt(1,ind - 1),curvePathPt(3,ind - 1),36,[0.4940,0.1840,0.5560]);
%         toolPt1 = fnval(toolSp1,0:0.001:1);
%         plot(toolPt1(1,:),toolPt1(3,:),'Color',[0.7,.7,.70]);
%         toolPt2 = fnval(toolSp2,0:0.001:1);
%         plot(toolPt2(1,:),toolPt2(3,:),'Color',[0.7,.7,.70]);
%         scatter(toolContactPt1(1),toolContactPt1(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(toolContactPt2(1),toolContactPt2(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(curveInterPt{ind}(1,:),curveInterPt{ind}(3,:),18,[0.850,0.325,0.0980],"filled");

    fprintf('No.%d\t toolpath point is calculated.\n-----\n',ind);
end

fprintf('The toolpath concentric optimization process causes %f seconds.\n',toc(tRes0));


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

nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on;
% colormap('summer');
% cb = colorbar;
% cb.Label.String = ['Height (',unit,')'];

% concentric toolpath for each loop
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
peakPt = [];
uLim = {zeros(2,0)};
interPt = {zeros(3,0)};
for ii = 1:size(curvePathPt,2)
    conTheta0 = linspace(conThetaBound(1),conThetaBound(2), ...
        ceil(2*pi/min(maxAngPtDist,arcLength/abs(curvePathPt(1,ii)))) + 1);
    conTheta0(end) = [];
    toolPathAngle = [toolPathAngle,conTheta0];
    loopPtNum = [loopPtNum,length(conTheta0)];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNum(end)];
    toolNAccum = [toolNAccum,ii*ones(1,loopPtNum(end))];
    toolRAccum = [toolRAccum,curvePt(1,ii)*ones(1,loopPtNum(end))];
    toolContactU = [toolContactU,curveContactU(ii)*ones(1,loopPtNum(end))];
    res = [res,curveRes(ii)*ones(1,loopPtNum(end))];
end
accumPtNum(1) = [];
for ii = 1:length(toolPathAngle)
    R = rotz(toolPathAngle(ii));
    kk = find(ii <= accumPtNum,1);
    loopQuat(ii,:) = rotm2quat(R);
    toolPathPt(:,ii) = R*curvePathPt(:,kk);
    peakPt(1:3,ii) = R*curvePeakPt(1:3,kk);
    peakPt(4:5,ii) = curvePeakPt(4:5,kk);
    uLim{ii} = curveULim{kk};
    interPt{ii} = R*curveInterPt{kk};
    toolQuat(ii,:) = quatmul(curveQuat(kk,:),loopQuat(ii,:));
    surfPt(:,ii) = R*curvePt(:,kk);
    toolNormDirect(ii,:) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    toolCutDirect(ii,:) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
end
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
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
% legend('designed surface','tool center point','tool edge','Location','northeast');

s6_visualize_concentric_multi;

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
toolREach = curvePt(1,:);
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
spiralRes = 5*aimRes*ones(1,spiralPtNum);
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
        [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D_multi( ...
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