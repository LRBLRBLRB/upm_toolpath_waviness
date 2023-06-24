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

% workspaceDir = fullfile('..','workspace','\20220925-contrast\nagayama_concentric';
% workspaceDir = fullfile('..','workspace','\20221020-tooltip\tooltip fitting result';
workspaceDir = fullfile('..','workspace','20230510');
diaryFile = fullfile('..','workspace','\diary',['diary',datestr(now,'yyyymmddTHHMMSS')]);
% diary diaryFile;
% diary on;

%% load tool file
[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)'; ...
    '*,*','all files(*.*)'}, ...
    'Select one tool edge data file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');
toolName = fullfile(dirName,fileName);
% toolName = 'output_data\tool\toolTheo_3D.mat';
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

%% load tool data file
[fileName,dirName,fileInd] = uigetfile({ ...
    '*.mat','text-files(*.mat)'; ...
    '*.nc;.pgm','CNC-files(*.nc,*pgm)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one cnc file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');

switch fileInd
    case 1
        cncName = fullfile(dirName,fileName);
        load(cncName);
    case 2
        cncName = fullfile(dirName,fileName);
        cncData = load(cncName,'-ascii');
        cncData = cncData.';
        % cncFid = fopen(cncName,'r');
        % % load the X Z C data
        % cncData = zeros(2,0);
        % while ~feof(cncFid)
        %     tmpLine = fgetl(cncFid);
        %     % if the line begins with %d%d or -%d, then break
        %     if strcmp(tmpLine,'(linking block)')
        %         break;
        %     end
        %     cncData = [cncData,sscanf(tmpLine,'C%fX%fZ%f')];
        % end
        % [ncData,nData] = fscanf(cncFid,'C%fX%fZ%f\n');
        % fclose(cncFid);
    otherwise
        return;
end

%% surface import

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 0.69/1000/(1000^(aimUnit - presUnit));
syms x y;
surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
surfFunc = matlabFunction(surfSym);
surfFx = diff(surfFunc,x);
surfFy = diff(surfFunc,y);
surfDomain = [-500,500;-500,500];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surfDomain = 1.1*surfDomain;
rMax = max(surfDomain(1,2),surfDomain(2,2));
% sampling density
spar = 501;
conR = linspace(0,rMax,spar); % concentric radius vector
conTheta = linspace(0,2*pi,spar);
[rMesh,thetaMesh] = meshgrid(conR,conTheta);
surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));

%% parameter settings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% machining paramters
cutDirection = 'Edge to Center'; % 'Center to Edge'
startDirection = 'X Minus'; % 'X Minus'
angularIncrement = 'Constant Arc'; % 'Constant Angle'
arcLength = 20; % um
maxAngPtDist = 1*pi/180;
angularLength = 1*pi/180;
radialIncrement = 'On-Axis'; % 'Surface'
aimRes = 1; % um
rStep = toolData.radius/2; % 每步步长可通过曲面轴向偏导数确定
maxIter = 100;
spiralMethod = 'Radius-Number'; % Radius-Angle
frMethodDefault = 'Approximation'; % 'Approximation'
frParamDefault = 1-1e-5;

% simulation parameters
isRecal = false; % whether to recalculate the residual error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isUIncrease = toolData.toolBform.coefs(end,1) - toolData.toolBform.coefs(1,1);
switch startDirection
    case 'X Plus' % plus in this program, but minus in moore
        rMax = max(surfDomain(1,2),surfDomain(2,2));
        rStep = -1*rStep;
    case 'X Minus' % minus in this program, but plus in moore
        rMax = min(surfDomain(1,1),surfDomain(2,1)); % reverse
        rStep = 1*rStep;
end

if isUIncrease*rStep < 0
    % ([1,0] & X minus) or ([0,1] & X plus)
    uDirection = 'U Minus';
    curvePlotSpar = -0.0001;
else
    % ([1,0] & X plus) or ([0,1] & X minus)
    uDirection = 'U Plus';
    curvePlotSpar = 0.0001;
end

switch cutDirection
    case 'Edge to Center'
        rRange = [rMax,0];
    case 'Center to Edge'
%             rRange = [0,rMax];
end

%% 2D simulation
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
% spiralRes = nan(2,spiralPtNum);
% spiralPeakPt = zeros(10,spiralPtNum);
% spiralInterPtIn = cell(1,spiralPtNum);
% spiralInterPtOut = cell(1,spiralPtNum);
% spiralULim = cell(1,spiralPtNum);
% tSpiralRes0 = tic;
% 
% parfor ind1 = 1:spiralPtNum
%     % inner ulim & residual height
%     ind2 = find(spiralAngle >= spiralAngle(ind1) + conThetaBound(end),1,'first');
%     ind3 = find(spiralAngle < spiralAngle(ind1) + conThetaBound(end),1,'last');
%     if isempty(ind2) || isempty(ind3)
% %         ind2 = find(spiralAngle >= spiralAngle(ind1) + pi,1,'first');
% %         ind3 = find(spiralAngle < spiralAngle(ind1) + pi,1,'last');
% %         [tmpres1,ptres1,spiralULim(:,ind1)] = residual3D( ...
% %             spiralPath,spiralNorm,spiralCut,spiralContactU, ...
% %             toolSp,toolRadius,spiralULim(:,ind1),ind1,ind2,ind3);
%         spiralULim{ind1} = [0;1];
%         tmpRes1 = nan;
%         tmpPeak1 = zeros(5,1);
%     else
%         [tmpRes1,tmpPeak1,spiralInterPtIn{ind1},spiralULim{ind1}] = ...
%             residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
%             toolData,toolRadius,spiralULim{ind1},ind1,ind2,ind3);
%     end
% 
%     % outer ulim & residual height
%     ind2 = find(spiralAngle >= spiralAngle(ind1) + 2*pi,1,'first');
%     ind3 = find(spiralAngle < spiralAngle(ind1) + 2*pi,1,'last');
%     if isempty(ind2) || isempty(ind3)
%         tmpRes2 = nan;
%         tmpPeak2 = zeros(5,1);
%     else
% 
% %         scatter3(spiralPath(1,ind1),spiralPath(2,ind1),spiralPath(3,ind1));
% %         quiver3(spiralPath(1,ind1),spiralPath(2,ind1),spiralPath(3,ind1), ...
% %             spiralCut(1,ind1),spiralCut(2,ind1),spiralCut(3,ind1));
% %         toolSp1 = toolData.toolBform;
% %         toolSp1.coefs = quat2rotm(spiralQuat(ind1,:))*toolSp1.coefs + spiralPath(:,ind1);
% %         fnplt(toolSp1);
% %         scatter3(spiralPath(1,ind2),spiralPath(2,ind2),spiralPath(3,ind2));
% %         quiver3(spiralPath(1,ind2),spiralPath(2,ind2),spiralPath(3,ind2), ...
% %             spiralCut(1,ind2),spiralCut(2,ind2),spiralCut(3,ind2),'AutoScale','on');
% %         scatter3(spiralPath(1,ind3),spiralPath(2,ind3),spiralPath(3,ind3));
% %         quiver3(spiralPath(1,ind3),spiralPath(2,ind3),spiralPath(3,ind3), ...
% %             spiralCut(1,ind3),spiralCut(2,ind3),spiralCut(3,ind3));
% 
%         [tmpRes2,tmpPeak2,spiralInterPtOut{ind1},spiralULim{ind1}] = ...
%             residual3D_multi(spiralPath,spiralNorm,spiralCut,spiralContactU, ...
%             toolData,toolRadius,spiralULim{ind1},ind1,ind2,ind3);
%     end
%     spiralRes(:,ind1) = [tmpRes1;tmpRes2];
%     spiralPeakPt(:,ind1) = [tmpPeak1;tmpPeak2];
% 
% %     if spiralRes(:,ind1) > 5
% %         1;
% %     end
% 
%     % debug
%     % plot3(toolPathPtRes(1,ii),toolPathPtRes(2,ii),toolPathPtRes(3,ii), ...
%     %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
%     % toolSp1 = toolSp;
%     % R1 = axesRot([0;0;1],[1;0;0],toolNormDirectRes(:,ii),toolCutDirectRes(:,ii),'zx');
%     % toolSp1.coefs = R1*toolSp.coefs + toolPathPtRes(:,ii);
%     % Q = fnval(toolSp1,uLimTmp(1,ii):0.01:uLimTmp(2,ii));
%     % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1);
% end
% 
% tSpiralRes = toc(tSpiralRes0);
% fprintf('The time spent in the residual height calculation for spiral toolpath process is %fs.\n',tSpiralRes);


