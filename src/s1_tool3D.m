% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高

%% simulation initialization
isAPP = false;
if isAPP
    addpath(genpath('funcs'));
    fitOpts.toolFitType = app.toolFitType;
    paramMethod = app.paramMethod;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;
else
    close all;
    clear; clc;
    isAPP = false;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = fullfile('..','workspace','/20220925-contrast/nagayama_concentric');
    fitOpts.arcFitMethod = 'levenberg-marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 14;
    textFontType = 'Times New Roman';
    
    debug = 3;
    switch debug
        case 2 % 2D tool profile simulation
            cx0 = 0*1000; % unit: mu m
            cy0 = 0*1000; % unit: mu m
            r0 = 0.1*1000; % unit: mu m
            noise = 0.03;
            zNoise = 0.01;
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
        case 3 % 3D tool profile simulation
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
        otherwise
            % tool measurement file loading
            [fileName,dirName] = uigetfile({ ...
                '*.mat','MAT-files(*.mat)'; ...
                '*.txt','text-files(*.txt)'; ...
                '*.*','all files(*.*)'...
                }, ...
                'Select one tool tip measurement data', ...
                workspaceDir, ...
                'MultiSelect','off');
            pathName = fullfile(dirName,fileName);
            dataOri = load(pathName); % load from file, n*2 matrix
            % data pre-processing
            toolOri = dataOri;
    end
end

%% plot the importing result
f1 = figure('Name','Original Scatters of the Tool');
% set(0,"Units","centimeters");
% winSize = get(0,"ScreenSize");
% set(gcf,"Units","centimeters");
% figSize = get(gcf,"Position");
% set(gcf,"Position",[(winSize(3) - figSize(3))/2,(winSize(4) - figSize(4))/2, ...
%     8.3,figSize(4)/figSize(3)*8.3]);
scatter3(toolOri(1,:),toolOri(2,:),toolOri(3,:));
hold on; grid on;
% axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
title('Original scatters of the tool');
view(-45,60);
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
fitOpts.arcFitdisplayType = 'iter-detailed';
[circ3D,toolFit] = toolFit3D(toolOri, ...
    'arcFitMethod',fitOpts.arcFitMethod,'arcFitdisplayType',fitOpts.arcFitdisplayType);
radius = circ3D.radius;
openAngle = circ3D.openAng;
pause(1);

% plot the fitting results
f2 = figure('Name','Tool Sharpness Fitting Result');
xLim = 1.1*max(toolFit(1,:));
quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
hold on;
text(0.9*xLim,-.05*radius,'x');
quiver(0,0.2*radius,0,-1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
text(0.05*xLim,-1.05*radius,'y');
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
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
clear theta xtmp ytmp

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

% rmpath(genpath('funcs'));