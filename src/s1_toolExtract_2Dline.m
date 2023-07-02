% to deal with the tool tip measurement data
% and to get the 3D point cloud of the tool tip

%% 2D curve results
isAPP = false;
if isAPP
    workspaceDir = app.workspaceDir;
    toolOri = app.toolOri;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;
    fitOpts.toolFitType = app.toolFitType;
    fitOpts.arcRansacMaxDist = app.arcRansacMaxDist;
    fitOpts.arcFitMethod = app.arcFitMethod;
    fitOpts.lineFitMaxDist = app.lineFitMaxDist;
    fitOpts.lineFitMethod = app.lineFitMethod;
else
%     close all;
    clear; clc;
    isAPP = false;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = fullfile('..','workspace','');
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';

    toolUnit = 'mm';
    
    unitList = {'m','mm','\mum','nm'};
    aimUnit = find(strcmp(unitList,unit),1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tool process parameters
    fitOpts.toolFitType = 'lineArc';
    fitOpts.arcRansacMaxDist = 1e-2;
    fitOpts.arcFitMethod = 'levenberg-marquardt';
    fitOpts.lineFitMaxDist = 1000^(aimUnit - 3);
    fitOpts.lineFitMethod = 'ransac';
    paramMethod = 'centripetal';

    radius0 = 192; % unit: um
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % tool measurement file loading
    [fileName,dirName] = uigetfile({ ...
        '*.csv','Comma-Separated Values-files(*.csv)'; ...
        '*.txt','text-files(*.txt)'; ...
        '*.*','all files(*.*)'...
        }, ...
        'Select one tool tip measurement data', ...
        workspaceDir, ...
        'MultiSelect','off');
    if ~fileName
        return;
    end
    pathName = fullfile(dirName,fileName);

    [~,~,fileExt] = fileparts(pathName);
    switch fileExt
        case {'.csv'}
            % get rid of the header of the csv file
            numHeader = 0;
            tooltipFile = fopen(pathName);
            while ~feof(tooltipFile)
                tmpLine = fgets(tooltipFile);
                % if the line begins with %d%d or -%d, then break
                if ~isnan(str2double(tmpLine(1:2)))
                    break;
                end
                numHeader = numHeader + 1;
            end
            fclose(tooltipFile);
            toolOri1 = importdata(pathName,',',numHeader);
            if size(toolOri1,2) ~= 3 && size(toolOri1,2) ~= 2
                toolOri1 = toolOri1.data;
            end
            toolOri1(:,3) = [];
%             toolOri2 = toolini(toolOri1);
            toolOri3 = sortrows(toolOri1,1,'ascend');
            toolOri = toolOri3';
        case {'.txt'}
            % tool data file that has been processed in mmt software
            toolOri1 = importdata(pathName,' ',0);
%             toolOri2 = toolini(toolOri1);
%             toolOri3 = sortrows(toolOri2,1,'ascend');
            toolOri = toolOri1';
    end

    % unit convertion
    presUnit = find(strcmp(unitList,toolUnit),1);
    toolOri = 1000^(aimUnit - presUnit)*toolOri;
end

%% plot the importing result
fig1 = figure('Name','Original tool data');
ax1 = plot(toolOri(1,:),toolOri(2,:),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);

%% outliners removed
% rmoutliers(toolOri,2,"median");
% 如果多边形的直线段和圆弧段重合部分的数据差别比较大，说明目前的测量模式有问题。


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
fitOpts.arcFitdisplayType = 'iter-detailed';
[circ2D,toolFitUnsorted,RMSE,lineFitMaxDist] = toolFit2D(toolOri,fitOpts.arcRansacMaxDist,fitOpts.lineFitMaxDist, ...
    'toolFitType',fitOpts.toolFitType,'lineFitMethod',fitOpts.lineFitMethod, ...
    'arcFitMethod',fitOpts.arcFitMethod,'radius0',radius0, ...
    'arcFitdisplayType',fitOpts.arcFitdisplayType);
radius = circ2D.radius;
openAngle = circ2D.openAng;

% tool data resort & averaging (to avoid loops)
toolAngle = atan2(toolFitUnsorted(2,:),toolFitUnsorted(1,:)); % polar angle of the tool point
[~,sortInd] = sort(toolAngle,'descend');
toolFit = toolFitUnsorted(:,sortInd);
% toolAngle1 = atan2(toolFit(2,:),toolFit(1,:)); % polar angle of the tool point
% toolDiff = toolAngle1(2:end) - toolAngle1(1:end - 1);
% find([0,toolDiff] > 0)

if isAPP
    app.lineFitMaxDist = lineFitMaxDist;
end

%% plot the fitting results
f2 = figure('Name','Tool Sharpness Fitting Result');
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
hold on;
xLim = 1.1*max(toolFit(1,:));
quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
hold on;
text(0.9*xLim,-.05*radius,'x');
quiver(0,0.2*radius,0,-1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
text(0.05*xLim,-1.05*radius,'y');
theta = (-pi/2 - openAngle/2):0.01:(-pi/2 + openAngle/2);
xtmp = radius*cos(theta);
ytmp = radius*sin(theta);
plot(xtmp,ytmp,'Color',[0.85,0.33,0.10],'LineWidth',1,'LineStyle','--'); % tool edge circle
scatter(0,0,'MarkerFaceColor',[0.85,0.33,0.10],'MarkerEdgeColor',[0.85,0.33,0.10]); % tool edge center
quiver(0,0,0,-0.5*radius,'AutoScale','off','Color',[0.93,0.69,0.13], ...
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

% clear xtmp ytmp theta xLim; % 删除画图的临时变量

nCPts = size(toolFit,2);
toolFit = [zeros(1,nCPts);toolFit];

%% post-processing
if isAPP
    app.toolFit = toolFit;
    app.openAngle = openAngle;
    app.radius = radius;
else
    s1_toolModel;
end

% rmpath(genpath('.')