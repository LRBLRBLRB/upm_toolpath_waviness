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
    close all;
    clear; clc;
    isAPP = false;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20221020-tooltip';
    fitOpts.toolFitType = 'lineArc';
    fitOpts.arcRansacMaxDist = 1e-3;
    fitOpts.arcFitMethod = 'levenberg-marquardt';
    fitOpts.lineFitMaxDist = 1e-3;
    fitOpts.lineFitMethod = 'polyfit';
    paramMethod = 'centripetal';
    unit = 'mm';
    textFontSize = 12;
    textFontType = 'Times New Roman';

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
                {'lineArc','centripetal','\mu m','Times New Roman','14'},optsInput);
            fitOpts.toolFitType = toolInput{1};
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
            % if the line begins with %d%d or -%d, then break
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
    'arcFitMethod',fitOpts.arcFitMethod, ...
    'arcFitdisplayType',fitOpts.arcFitdisplayType);
radius = circ2D.radius;
openAngle = circ2D.openAng;

% tool data resort & averaging
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