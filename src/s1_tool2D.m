% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高

%% simulation initialization
isAPP = false;
if isAPP
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
    workspaceDir = fullfile('..','workspace','/20221020-tooltip';
    fitOpts.arcFitMethod = 'levenberg-marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';

    debug = 2;
    switch debug
        case 2 % 2D tool profile simulation
            cx0 = 0*1000; % unit: mu m
            cy0 = 0*1000; % unit: mu m
            r0 = 0.1*1000; % unit: mu m
            openAng = pi/3; % unit: rad
            edgePV = 200; % low-frequency error
            k = -edgePV/openAng;
            noise = r0*5e-3; % mid-frequency error
            theta = linspace(0,openAng,200);
            r = r0 + edgePV/2 + k*theta + (noise*rand(1,length(theta)) - 0.5*noise);
            toolOri(1,:) = r.*cos(theta) + cx0;
            toolOri(2,:) = r.*sin(theta) + cy0;
            rmse0 = sqrt(...
                        sum(...
                            ((toolOri(1,:) - cx0).^2 + (toolOri(2,:) - cy0).^2 - r0^2).^2) ...
                    /length(theta));
            clear theta r;
        otherwise
            toolInput = inputdlg({'Tool fitting method','Tool parameterization method', ...
                'Unit','Font type in figure','Font size in figure'}, ...
                'Tool Processing Input', ...
                [1 20; 1 20; 1 20; 1 20; 1 20], ...
                {'Levenberg-Marquardt','centripetal','\mu m','Times New Roman','14'}, ...
                'WindowStyle');
            fitOpts.arcFitMethod = toolInput{1};
            paramMethod = toolInput{2};
            unit = toolInput{3};
            textFontType = toolInput{4};
            textFontSize = str2double(toolInput{5});
            % tool measurement file loading
            [fileName,dirName] = uigetfile({ ...
                '*.mat','MAT-files(*.mat)'; ...
                '*.txt','text-files(*.txt)'; ...
                '*.*','all files(*.*)'...
                }, ...
                'Select one tool tip measurement data', ...
                'Tool\', ...
                'MultiSelect','off');
            pathName = fullfile(dirName,fileName);
            dataOri = load(pathName); % load from file, n*2 matrix
            % data pre-processing
            toolOri = dataOri;
    end
end

%% plot the importing result
% dim = size(toolOri,1);
f1 = figure('Name','original xy scatters of the tool');
plot(toolOri(1,:),toolOri(2,:),'.','MarkerSize',2);
hold on;

% theta1 = atan2(toolOri(2,1),toolOri(1,1));
% theta2 = atan2(toolOri(2,end),toolOri(1,end));
% theta = transpose(linspace(theta1,theta2,100));
% plot(cx0 + r0*cos(theta),cy0 + r0*sin(theta),'blue');
% scatter(cx0,cy0,'blue');
axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
fitOpts.arcFitdisplayType = 'iter-detailed';
[circ2D,toolFit,rmseLsc] = toolFit2D(toolOri, ...
    'arcFitMethod',fitOpts.arcFitMethod,'arcFitdisplayType',fitOpts.arcFitdisplayType);
radius = circ2D.radius;
openAngle = circ2D.openAng;

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
clear xtmp ytmp theta; % 删除画图的临时变量

nCPts = size(toolFit,2);
toolFit = [zeros(1,nCPts);toolFit];

%% tool modelling
s1_toolModel

%% 相邻刀位点的残高
% x1 = [0,0,0]; vec1 = [0,1];
% x2 = [r0,r0/30,0]; vec2 = [cosd(90+8),sind(90+8)];
% R = rotz(5,'deg'); % rotation matrix about z axis
% tool1 = toolPt + x1;
% tool2 = toolPt*transpose(R) + x2;
% 
% [res,interPt] = residualHigh(x1(:,1:2),vec1,tool1(:,1:2),x2(:,1:2),vec2,tool2(:,1:2));
% 
% figure('Name','Residual of the adjacent tool');
% plot(tool1(:,1),tool1(:,2),'Color',[0,0.45,0.74]); hold on;
% plot(tool2(:,1),tool2(:,2),'Color',[0.85,0.33,0.10]);
% plot(interPt(1),interPt(2),'*','Color',[0.49,0.18,0.56]);
% plot(x1(1),x1(2),'o','Color',[0,0.45,0.74],'MarkerFaceColor',[0,0.45,0.74]);
% line([x1(1),x1(1)+radius*vec1(1)],[x1(2),x1(2)+radius*vec1(2)],'LineWidth',2,'Color',[0,0.45,0.74]);
% plot(x2(1),x2(2),'o','Color',[0.85,0.33,0.10],'MarkerFaceColor',[0.85,0.33,0.10]);
% line([x2(1),x2(1)+radius*vec2(1)],[x2(2),x2(2)+radius*vec2(2)],'LineWidth',2,'Color',[0.85,0.33,0.10]);
% axis equal

%%
% rmpath(genpath('funcs'));