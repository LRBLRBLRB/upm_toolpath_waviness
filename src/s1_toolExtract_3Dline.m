% to deal with the tool tip measurement data
% and to get the 3D point cloud of the tool tip

%% 3D curve results
isAPP = true;
if isAPP
    workspaceDir = app.workspaceDir;
    toolOri = app.toolOri;
    fitOpts.toolFitType = app.toolFitType;
    paramMethod = app.paramMethod;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;
else
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = '../workspace/20221020-tooltip';
    fitOpts.fitMethod = 'Levenberg-Marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';

    isTest = false;
    isSelf = false;
    if isTest
        cx0 = 1*1000; % unit:mu m
        cy0 = 2*1000; % unit:mu m
        cz0 = 3*1000; % unit:mu m
        r0 = 0.1*1000; % unit:mu m
        openAng = pi/3; % unit: rad
        edgePV = 200; % low-frequency error
        k = -edgePV/openAng;
        noise = r0*5e-3; % mid-frequency error
        zNoise = r0*0.05; % data pre-processing error
        theta = linspace(0,openAng,300);
        r = r0 + edgePV/2 + k*theta + (noise*rand(1,length(theta)) - 0.5*noise);
        toolOri(1,:) = cx0 + r.*cos(theta);
        toolOri(2,:) = cy0 + r.*sin(theta);
        toolOri(3,:) = cz0 + (zNoise*rand(1,length(theta)) - 0.5*zNoise);
        rmse0 = sqrt( ...
            sum((toolOri - ndgrid([cx0;cy0;cz0],1:length(theta)).^2),2) ...
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
            if ~isnan(str2double(tmpLine(1:2)))
                break;
            end
            numHeader = numHeader + 1;
        end
        fclose(tooltipFile);
        toolOri = importdata(pathName,',',numHeader);
        toolOri = toolOri.data;
        toolOri = toolOri';
    end
end

figure('Name','Original tool data');
plot3(toolOri(1,:),toolOri(2,:),toolOri(3,:),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

%% outliners removed
% rmoutliers(toolOri,2,"median");
% 如果多边形的直线段和圆弧段重合部分的数据差别比较大，说明目前的测量模式有问题。


%% ransac line fitting
% it seems to be useless at present, since the tool tip arc can be
% recognized by directly ransac circle fitting.

% sampleSz = 2; % number of points to sample per trial
% maxDist = 2; % max allowable distance for inliers
% 
% fitLineFcn = @(pts) polyfit(pts(:,1),pts(:,2),1); % fit function using polyfit
% evalLineFcn = ...   % distance evaluation function
%   @(mdl, pts) sum((pts(:, 2) - polyval(mdl, pts(:,1))).^2,2);
% 
% [modelRANSAC,inlierIdx] = ransac(oriPts,fitLineFcn,evalLineFcn, ...
%   sampleSz,maxDist);
% 
% inlierPts = oriPts(inlierIdx,:);
% x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
% y = modelRANSAC(1)*x + modelRANSAC(2);
% plot(x, y, 'm-');

%% tool tip arc fitting and extraction

% ransac
sampleSz = 3; % number of points to sample per trial
maxDist = 100; % max allowable distance for inliers

fitLineFcn = @(pts) arcFit3D(pts','displayType','off');  % fit function
evalLineFcn = ...   % distance evaluation function
  @(mdl, pts) sum(abs(vecnorm(pts - (mdl{1})',2,2) - mdl{2}^2));


[modelRANSAC,inlierIdx] = ransac(toolOri',fitLineFcn,evalLineFcn, ...
  sampleSz,maxDist);



%% plot the fitting results
% f2 = figure('Name','Tool Sharpness Fitting Result');
% xLim = 1.1*max(toolFit(1,:));
% quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
% hold on;
% text(0.9*xLim,-.05*radius,'x');
% quiver(0,-0.2*radius,0,1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
% text(0.05*xLim,1.05*radius,'y');
% plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
% theta = (pi/2 - openAngle/2):0.01:(pi/2 + openAngle/2);
% xtmp = radius*cos(theta);
% ytmp = radius*sin(theta);
% plot(xtmp,ytmp,'Color',[0.85,0.33,0.10],'LineWidth',1,'LineStyle','--'); % tool edge circle
% scatter(0,0,'MarkerFaceColor',[0.85,0.33,0.10],'MarkerEdgeColor',[0.85,0.33,0.10]); % tool edge center
% quiver(0,0,0,0.5*radius,'AutoScale','off','Color',[0.93,0.69,0.13], ...
%     'LineWidth',2.5,'MaxHeadSize',0.3); % tool edge normal
% line([0,xtmp(1)],[0,ytmp(1)],'LineStyle','--','Color',[0.85,0.33,0.10]);
% line([0,xtmp(end)],[0,ytmp(end)],'LineStyle','--','Color',[0.85,0.33,0.10]);
% % xlim([-1.1*xLim,1.1*xLim]);
% axis equal;
% % set(gca,'TickLabelInterpreter','tex');
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel(['x(',unit,')']);
% ylabel(['y(',unit,')']);
% title('Tool fitting result');
% legend('','','tool edge','tool fitting arc','tool center', ...
%     'tool normal vector','Location','northeast');


%% tool modelling
% s1_toolModel

% rmpath(genpath('.')