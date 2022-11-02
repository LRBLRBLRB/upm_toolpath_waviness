% to deal with the tool tip measurement data
% and to get the 3D point cloud of the tool tip

if true
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20221020-tooltip';
    fitOpts.arcFitMethod = 'Levenberg-Marquardt';
    paramMethod = 'centripetal';
    unit = 'mm';
    textFontSize = 12;
    textFontType = 'Times New Roman';
end

%% 3D curve results
isTest = false;
isSelf = false;
if isTest
    cx0 = 0*1000; % unit: mu m
    cy0 = 0*1000; % unit: mu m
    r0 = 0.1*1000; % unit: mu m
    noise = 0.03;
    theta = linspace(0,2*pi/3,200);
    r = r0*(1 - noise + 2*noise*rand(1,length(theta)));
    toolOri(1,:) = r.*cos(theta) + cx0;
    toolOri(2,:) = r.*sin(theta) + cy0;
    toolOri(3,:) = zeros(1,length(theta));
    rmse0 = sqrt(...
                sum(...
                    ((toolOri(1,:) - cx0).^2 + (toolOri(2,:) - cy0).^2 - r0^2).^2) ...
            /length(theta));
    clear theta r;
else
    if isSelf
        optsInput.Resize = 'on';
        optsInput.WindowStyle = 'normal';
        optsInput.Interpreter = 'tex';
        toolInput = inputdlg({'Tool fitting method','Tool parameterization method', ...
            'Unit','Font type in figure','Font size in figure'}, ...
            'Tool Processing Input', ...
            [1 50; 1 50; 1 50; 1 50; 1 50], ...
            {'Levenberg-Marquardt','centripetal','\mu m','Times New Roman','14'},optsInput);
        fitOpts.arcFitMethod = toolInput{1};
        paramMethod = toolInput{2};
        unit = toolInput{3};
        textFontType = toolInput{4};
        textFontSize = str2double(toolInput{5});
        % tool measurement file loading
        [fileName,dirName] = uigetfile({ ...
            '*.csv','Comma-Separated Values-files(*.csv)'; ...
            '*.mat','MAT-files(*.mat)'; ...
            '*.txt','text-files(*.txt)'; ...
            '*.*','all files(*.*)'...
            }, ...
            'Select one tool tip measurement data', ...
            workspaceDir, ...
            'MultiSelect','off');
        pathName = fullfile(dirName,fileName);
    else
        pathName = fullfile(workspaceDir,"tooltip result/20221019-strategy-2+40-5.csv");
    end
    % get rid of the header of the csv file
    numHeader = 0;
    tooltipFile = fopen(pathName);
    while ~feof(tooltipFile)
        tmpLine = fgets(tooltipFile);
        if ~isnan(str2double(tmpLine(1:2)))
            break;
        end
        numHeader = numHeader + 1;
    end
    fclose(tooltipFile);
    toolOri = importdata(pathName,',',numHeader);
    toolOri = toolOri.data;
    toolOri(:,3) = [];
    toolOri = toolOri';
end


fig1 = figure('Name','Original tool data');
ax1 = plot(toolOri(1,:),toolOri(2,:),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);

% outliners removed
% rmoutliers(oriPts,2,"median");
% 如果多边形的直线段和圆弧段重合部分的数据差别比较大，说明目前的测量模式有问题。


%% tool tip arc fitting and extraction by RANSAC

% % ransac
% sampleSz = 3; % number of points to sample per trial
% maxDist = 0.003; % max allowable distance for inliers
% 
% fitLineFcn = @(pts) arcFit2D(pts','displayType','off');  % fit function
% evalLineFcn = ...   % distance evaluation function
%   @(mdl, pts) abs(vecnorm(pts - (mdl{1})',2,2) - mdl{2});
% 
% % test whetger the functions above is true
% isTest = false;
% if isTest
%     cx0 = 0*1000; % unit: mu m
%     cy0 = 0*1000; % unit: mu m
%     r0 = 0.1*1000; % unit: mu m
%     noise = 0.03;
%     theta = (linspace(0,2*pi/3,200))';
%     r = r0*(1 - noise + 2*noise*rand(length(theta),1));
%     oriPts(1,:) = r.*cos(theta) + cx0;
%     oriPts(2,:) = r.*sin(theta) + cy0;
%     rmse0 = sqrt(...
%                 sum(...
%                     ((oriPts(1,:) - cx0).^2 + (oriPts(2,:) - cy0).^2 - r0^2).^2) ...
%             /length(theta));
%     fitCirc = fitLineFcn(oriPts');
%     figure('Name','Function Testification');
%     plot(oriPts(1,:),oriPts(2,:),'.','Color',[0,0.45,0.74]); hold on;
%     % plot the fitting center of the circle
%     scatter(fitCirc{1}(1),fitCirc{1}(2),36,[0.6350,0.0780,0.1840],'filled');
%     quiver(fitCirc{1}(1),fitCirc{1}(2), ...
%         0.6*fitCirc{2}*fitCirc{4}(1),0.6*fitCirc{2}*fitCirc{4}(2), ...
%         'filled','Color',[0.6350,0.0780,0.1840]);
%     % plot the fitting circle
%     scaThe = linspace(0,2*pi);
%     scat(1,:) = fitCirc{2}*cos(scaThe) + fitCirc{1}(1);
%     scat(2,:) = fitCirc{2}*sin(scaThe) + fitCirc{1}(2);
%     plot(scat(1,:),scat(2,:),'k--','LineWidth',1);
%     % plot the fitting arc
%     R = rotz(atan2(fitCirc{4}(2),fitCirc{4}(1)) - 0.5*fitCirc{3});
%     scaThe = linspace(0,fitCirc{3});
%     scat(1,:) = fitCirc{2}*cos(scaThe);
%     scat(2,:) = fitCirc{2}*sin(scaThe);
%     circFit = R(1:2,1:2)*scat + fitCirc{1};
%     plot(circFit(1,:),circFit(2,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',3);
%     hold off;
%     grid on;
%     axis equal
%     evalPer = evalLineFcn(fitCirc,oriPts');
% %     evalSum = sum(evalLineFcn(fitCirc,toolOri));
% end
% 
% [modelRANSAC,inlierInd] = ransac(toolOri',fitLineFcn,evalLineFcn, ...
%   sampleSz,maxDist);
% 
% % Least Square Fitting Based on the Inliers
% fitOpts.arcFitdisplayType = 'iter-detailed';
% [circ2D,toolFit,rmseLsc] = toolFit2D(toolOri(:,inlierInd), ...
%     'arcFitMethod',fitOpts.arcFitMethod,'displayType',fitOpts.arcFitdisplayType);
% radius = circ2D{2};
% openAngle = circ2D{3};

%% tool tip arc fitting and extraction by self-selection

% line fitting based on ransac
figure(fig1);
% get the start and end point of the line
fprintf(['Please select the points of the two edegs of the tool tip: \n' ...
    'The first two are the start and end points of the left one, ' ...
    'while the last two are those of the right.\n'])

% we have three methods to get the points
% [lineX,lineY] = ginput(2);
% [lineX,lineY] = getpts(fig1);
% ax1.Children

lineMid = ginput(2);
leftPts = toolOri(:,toolOri(1,:) < lineMid(1,1));
rightPts = toolOri(:,toolOri(1,:) > lineMid(2,1));

% ransac line fitting
sampleSz = 2; % number of points to sample per trial
maxDist = 0.00001; % max allowable distance for inliers
fitLineFcn = @(pts) polyfit(pts(:,1),pts(:,2),1); % fit function using polyfit
evalLineFcn = ...   % distance evaluation function
  @(mdl, pts) sum((pts(:, 2) - polyval(mdl, pts(:,1))).^2,2);
[leftLine,leftInlierIdx] = ransac(leftPts',fitLineFcn,evalLineFcn, ...
  sampleSz,maxDist);
[rightLine,rightInlierIdx] = ransac(rightPts',fitLineFcn,evalLineFcn, ...
  sampleSz,maxDist);

fig2 = figure('Name','Remove Outliers of the Tool Tip Arc');

grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
plot(leftPts(1,leftInlierIdx),leftLine(1)*leftPts(1,leftInlierIdx) + leftLine(2), ...
    '.','Color',[0.9290    0.6940    0.1250],'MarkerSize',8);
hold on;
plot(rightPts(1,rightInlierIdx),rightLine(1)*rightPts(1,rightInlierIdx) + rightLine(2), ...
    '.','Color',[0.9290    0.6940    0.1250],'MarkerSize',8);
plot(toolOri(1,:),toolOri(2,:),'LineWidth',0.5,'Color',[0    0.4470    0.7410]);

% Least Square Fitting Based on the Inliers
circIdx(1) = find(leftInlierIdx,1,"last");
circIdx(2) = size(toolOri,2) - find(flipud(rightInlierIdx),1,"last");
fitOpts.arcFitdisplayType = 'iter-detailed';
[circ2D,toolFit,rmseLsc] = toolFit2D(toolOri(:,circIdx(1):circIdx(2)), ...
    'arcFitMethod',fitOpts.arcFitMethod,'displayType',fitOpts.arcFitdisplayType);
radius = circ2D{2};
openAngle = circ2D{3};

%% plot the fitting results
f2 = figure('Name','Tool Sharpness Fitting Result');
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
hold on;
xLim = 1.1*max(toolFit(1,:));
quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
hold on;
text(0.9*xLim,-.05*radius,'x');
quiver(0,-0.2*radius,0,1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
text(0.05*xLim,1.05*radius,'y');
theta = (pi/2 - openAngle/2):0.01:(pi/2 + openAngle/2);
xtmp = radius*cos(theta);
ytmp = radius*sin(theta);
plot(xtmp,ytmp,'Color',[0.85,0.33,0.10],'LineWidth',1,'LineStyle','--'); % tool edge circle
scatter(0,0,'MarkerFaceColor',[0.85,0.33,0.10],'MarkerEdgeColor',[0.85,0.33,0.10]); % tool edge center
quiver(0,0,0,0.5*radius,'AutoScale','off','Color',[0.93,0.69,0.13], ...
    'LineWidth',2.5,'MaxHeadSize',0.3); % tool edge normal
line([0,xtmp(1)],[0,ytmp(1)],'LineStyle','--','Color',[0.85,0.33,0.10]);
line([0,xtmp(end)],[0,ytmp(end)],'LineStyle','--','Color',[0.85,0.33,0.10]);
% xlim([-1.1*xLim,1.1*xLim]);
axis equal;
% set(gca,'TickLabelInterpreter','tex');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
title('Tool fitting result');
legend('','','tool edge','tool fitting arc','tool center', ...
    'tool normal vector','Location','northeast');


%% tool modelling
% s1_toolModel

% rmpath(genpath('.')