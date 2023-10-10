% simulation process for the (3-axes) STS diamond turning spiral tool path
%
% A cnc file is expected to load, and the machining results as well as the
% residual error distribution map will be generated.
% Notice that the post process is suitable for the Nanotech 650FG V2 only.

isAPP = false;
if isAPP
    return;
else
    %% function-used
%     close all;
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
    % profile on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % design machining paramters
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
    % load cnc file
    cncMode = 'C%f X%f Z%f';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % design concentric surface generation / import
    % A = tand(20)/(2*2000);
    c = 0.69/1000/(1000^(aimUnit - presUnit));
    param = sprintf('c = %f',c);
    syms C;
    surfSymDisp = C.*(x.^2 + y.^2)./(1 + sqrt(1 - C.^2.*(x.^2 + y.^2)));
    surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfNormFunc = matlabFunction([surfFx;surfFy;-1],'Vars',{'x','y'});
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

if isUIncrease*rStep > 0 % the direction of the parameter u while feeding
    % ([1,0] & X Plus) or ([0,1] & X Minus)
    uDirection = 'U Plus';
else
    % ([1,0] & X Minus) or ([0,1] & X Plus)
    uDirection = 'U Minus';
end

switch cutDirection
    case 'Edge to Center'
        rRange = [rMax,0];
        if strcmp(startDirection,'X Plus') % 'X Minus'
            conThetaBound = [0,-2*pi];
            uLimOrder2 = 'last';
            uLimOrder3 = 'first';
        else
            conThetaBound = [0,2*pi];
            uLimOrder2 = 'first';
            uLimOrder3 = 'last';
        end
    case 'Center to Edge'
        rRange = [0,rMax];
        if strcmp(startDirection,'X Plus') % 'X Minus'
            conThetaBound = [0,2*pi];
            uLimOrder2 = 'first';
            uLimOrder3 = 'last';
        else
            conThetaBound = [0,-2*pi];
            uLimOrder2 = 'last';
            uLimOrder3 = 'first';
        end
end

% load nc file
cncData = read_my_STS(workspaceDir,cncMode);
% convert the unit
cncData = unitconversion({'angle','length','length'},cncData, ...
    'mm','\mum','deg','rad');
% work out the curveQuat (rotate the tool from standard position to that
% perpenticular to the Z-axis
% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
curveRot = axesRot(toolData.toolEdgeNorm,toolData.cutDirect, ...
    [0;0;-1],cutDirect,'zx');
curveQuat = rotm2quat(curveRot);
% forward kinematics
[spiralAngle,spiralPath,spiralQuat,spiralNorm,spiralCut] = ...
    moore650kine('CXZ',cncData,curveQuat,conThetaBound, ...
    'shift',true,'toolData',toolData);

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
plot3(spiralPath(1,:),spiralPath(2,:),spiralPath(3,:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');

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

%% 2D simulation
t0 = tic;
if exist('curvePathPt','var')
    % 2D simulation
    if isRecal
        curveNum = size(curvePathPt,2);
        curveRes = 5*aimRes*ones(1,curveNum);
        curvePeakPt = zeros(5,curveNum);
        curveInterPt = cell(1,34);
        curveInterPt{1} = [0;0;0];
        curveULim = cell(1,34);
        switch uDirection % the interval of each toolpath
            case 'U Plus'
                curveULim{1} = [0;1];
            case 'U Minus'
                curveULim{1} = [1;0];
        end
        toolSp = toolData.toolBform;
        for ind = 2:curveNum
            tic
            % calculate the residual height of the loop and the inner nearest loop
            toolSp1 = toolSp;
            toolSp1.coefs = quat2rotm(curveQuat(ind,:))*toolSp1.coefs + curvePathPt(:,ind);
            % toolContactPt1 = fnval(toolSp1,curveContactU(ind));
            toolSp2 = toolSp;
            toolSp2.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp2.coefs + curvePathPt(:,ind - 1);
            % toolContactPt2 = fnval(toolSp2,curveContactU(ind - 1));
        
            [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim1, ...
                curveULim2] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
                curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1}, ...
                'uDirection',uDirection,'aimRes',aimRes);

            curveULim{ind} = curveULim1;
            curveULim{ind - 1} = curveULim2;
            fprintf('No.%d\t toolpath %f\t is calculated within %fs.\n-----\n',ind,curvePathPt(1,ind),toc);
        end
    end
    figure('Name','tool path optimization');
    tiledlayout(1,2);
    pos = get(gcf,'position');
    set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
    nexttile;
    plot(curvePathPt(1,:),curvePathPt(3,:),'.','Color',[0.4940,0.1840,0.5560]); % curve path points
    hold on;
    plot(curvePt(1,:),curvePt(3,:),'.','Color',[0.9290,0.6940,0.1250]); % curve contact points
    rSpar = linspace(0,rMax,spar);
    plot(rSpar,surfFunc(rSpar,zeros(1,length(rSpar))),'Color',[0.7,0.7,0.7],'LineWidth',0.3); % curve
    for jj = 1:size(curvePathPt,2)
        plot(curveInterPt{jj}(1,:),curveInterPt{jj}(3,:),'.','Color',[0.850,0.3250,0.0980]); % intersection point
        toolSp1 = toolData.toolBform;
        toolSp1.coefs = quat2rotm(curveQuat(jj,:))*toolSp1.coefs + curvePathPt(:,jj);
        for ii = 1:size(curveULim{jj},2)
            toolSp1Pt = fnval(toolSp1,curveULim{jj}(1,ii):curvePlotSpar:curveULim{jj}(2,ii));
            plot(toolSp1Pt(1,:),toolSp1Pt(3,:),'Color',[0,0.4470,0.7410], ...
                'LineWidth',0.5);
%             toolSp1Pt = fnval(toolSp1,curveULim{jj}(2,ii):curvePlotSpar:curveULim{jj}(1,ii + 1));
%             toolSp1Pt(3,end) = NaN;
%             patch('XData',toolSp1Pt(1,:),'YData',toolSp1Pt(3,:), ...
%                 'EdgeColor',[0,0.4470,0.7410],'EdgeAlpha',0.1, ...
%                 'LineWidth',0.5,'LineStyle','-');
        end
    end
    legend('tool path point','tool contact point','ideal surface','peak point','','actual surface');
    set(gca,'FontSize',textFontSize,'FontName',textFontType);
    xlabel(['r (',unit,')']);
    ylabel(['z (',unit,')']);
    nexttile;
    plot(curveRes(2:end),'.','Color',[0,0.4470,0.7410],'MarkerSize',18); hold on;
    line([0;curveNum + 1],[mean(curveRes(2:end));mean(curveRes(2:end))], ...
        'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    xlabel('No. of CL points');
    ylabel(['Residual error (',unit,')']);
    grid on;
end

%% 2-axes convertion to cartesian
% spiralPtNum = size(cncData,2);
% spiralPath = zeros(3,spiralPtNum);
% spiralQuat = zeros(spiralPtNum,4);
% 
% spiralAngle = cncData(1,:)*pi/180;
% angAdd = find(diff(spiralAngle) < 0);
% for ii = 1:length(angAdd)
%     spiralAngle(angAdd(ii) + 1:end) = spiralAngle(angAdd(ii) + 1:end) + 2*pi;
% end
% 
% figure;
% xxx = linspace(0,surfDomain(1,2),1000);
% zzz = A*xxx.^2 + C;
% plot(xxx,zzz);
% hold on;
% curvePath = cncData(:,angAdd);
% plot(curvePath(2,:),curvePath(3,:));
% 
% spiralPath(3,:) = cncData(3,:);
% 
% spiralPath(1,1) = cncData(2,1);
% spiralQuat(1,:) = [1,0,0,0];
% for ii = 2:spiralPtNum
%     spiralQuat(ii,:) = rotm2quat(rotz(cncData(1,ii)));
%     spiralPath(1,ii) = cncData(2,ii)*cos(spiralAngle(ii));
%     spiralPath(2,ii) = cncData(2,ii)*sin(spiralAngle(ii));
% end
% % spiralPath = 1000*spiralPath;
% 
% figure;
% surf(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
% hold on;
% plot3(spiralPath(1,:),spiralPath(2,:),spiralPath(3,:),'-','Color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% 
% 
% % for ii = 1:spiralPtNum
% %     toolSp1 = toolData.toolBform;
% %     toolSp1.coefs = toolSp1.coefs + spiralPath(:,ii);
% %     toolPt1 = fnval(toolSp1,0:0.01:1);
% %     plot3(toolPt1(1,:),toolPt1(2,:),toolPt1(3,:),'Color',[0,0.4450,0.7410]);
% % end

%% spiral tool path simulation and residual height calculation
spiralPtNum = size(cncData,2);
spiralPathZ = zeros(3,spiralPtNum);
spiralContactU = zeros(1,spiralPtNum);
surfPt = zeros(3,spiralPtNum);

% method 1: directly solve the z-coordinate of CL points before residual calculation
for ind1 = 1:spiralPtNum
    [spiralPathZ(:,ind1),~,spiralContactU(ind1),surfPt(:,ind1)] = toolpathpos( ...
        surfFunc,surfNormFunc,toolData,spiralPath(:,ind1),spiralNorm(:,ind1), ...
        spiralCut(:,ind1),'directionType','norm-cut');
    if abs(spiralPathZ(3,ind1) - spiralPath(3,ind1)) > 1
        fprintf("spiralPathZ = %f, spiralPath = %f\n",spiralPathZ(3,ind1),spiralPath(3,ind1));
        pause;
    else
        disp(ind1);
    end
end
figure;
plot3(spiralPath(1,:),spiralPath(2,:),spiralPathZ(:), ...
    'Color',[0,0.4470,0.7410],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
hold on;
plot3(surfPt(1,:),surfPt(2,:),surfPt(:), ...
    'Color',[0.8500 0.3250 0.0980],'LineStyle',':','LineWidth',0.1, ...
    'Marker','.','MarkerSize',6);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.2,'LineStyle','none');
colormap('summer');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
% set(gca,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('designed surface','tool center point','Location','northeast');
drawnow;


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

%%
% delete(parObj);
diary off;
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));