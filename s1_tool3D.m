% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高
if true
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
    FitMethod = 'Levenberg-Marquardt';
    ParamMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
end

%% simulation initialization
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
        toolInput = inputdlg({'Tool fitting method','Tool parameterization method', ...
            'Unit','Font type in figure','Font size in figure'}, ...
            'Tool Processing Input', ...
            [1 20; 1 20; 1 20; 1 20; 1 20], ...
            {'Levenberg-Marquardt','centripetal','\mu m','Times New Roman','14'}, ...
            'WindowStyle');
        FitMethod = toolInput{1};
        ParamMethod = toolInput{2};
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
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
zlabel(['z(',unit,')']);
title('Original scatters of the tool');
view(-45,60);
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
nCPts = size(toolOri,2);
[radius,openAngle,toolFit] = toolFit3D(toolOri,FitMethod);

% plot the fitting results
f2 = figure('Name','Tool Sharpness Fitting Result');
xLim = 1.1*max(toolFit(1,:));
quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
hold on;
text(0.9*xLim,-.05*radius,'x');
quiver(0,-0.2*radius,0,1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
text(0.05*xLim,1.05*radius,'y');
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
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


%% 分离轮廓误差和波纹度误差
% Cartesian coordinate to cylindrical coordinate
toolTheta = atan2(toolFit(2,:),toolFit(1,:));
toolR = vecnorm(toolFit,2,1);

% plot the geometric error
figure('name','Tool Geometric Error');
t = tiledlayout(2,1);
nexttile;
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74],'LineWidth',0.5); % tool edge scatters
hold on;
plot(xtmp,ytmp,'Color',[0.85,0.33,0.10],'LineWidth',1,'LineStyle','--'); % tool edge circle
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
title('Tool contour');
legend('tool edge','tool fitting arc','Location','northeast');

nexttile;
plot(toolTheta*180/pi - 90,toolR - radius); hold on;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlim([-openAngle*180/pi/2,openAngle*180/pi/2]);
title('Tool geometric error');
ylabel({'polar diameter error',['(',unit,')']});
xlabel('central angle \theta(°)');
grid on;

title(t,'Tool Geometric Error');

% fft to filter different geometric error 



% plot the sharness and waviness, respectively
% nexttile;
% plot(toolTheta*180/pi,toolR); hold on;
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% title('Tool contour error');
% xlabel('central angle \theta(°)');
% ylabel({'polar diameter',['(',unit,')']});

% nexttile;
% plot(toolTheta*180/pi - 90,toolR - radius); hold on;
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlim([-openAngle*180/pi/2,openAngle*180/pi/2]);
% title('Tool waviness error');
% ylabel({'polar diameter error',['(',unit,')']});
% xlabel('central angle \theta(°)');
% grid on;

clear xtmp ytmp theta xLim; % 删除画图的临时变量

%% 车刀轮廓插值 two methods to interpolate
k = 3; % degree of the B-spline
u = 0:0.0002:1;
nPts = length(u);

toolFit = [zeros(1,nCPts);toolFit];
% B-spline interpolate in the polar coordinate
[toolEdgePt,toolBform] = bsplinePts_spapi(toolFit,k,u, ...
    'ParamMethod',ParamMethod,'CoordinateType','Cartesian'); 
toolCpts = toolBform.coefs;

% [toolCpts,U] = bSplineCpts(toolFit',k,'chord');
% toolDim = size(toolCpts,2);
% toolPt2 = zeros(nPts,toolDim);
% for i = 1:nPts
%     toolPt2(i,:) = bSplinePt(toolCpts,k,u(i),U);
% end
% toolPt2 = toolPt2';
% 
% plot the difference between the two
% figure('Name','Comparison between the self function and the Mathworks ons');
% tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
% nexttile;
% plot(u',toolPt2(1,:) - toolPt(1,:)); hold on;
% ylabel('\Delta x');
% nexttile;
% plot(u',toolPt2(2,:) - toolPt(2,:)); hold on;
% ylabel('\Delta y');
% xlabel('u');

% plot the interpolation results
figure('Name','Tool Interpolation Results');
plot(toolFit(2,:),toolFit(3,:),'--.', ...
    'MarkerSize',8,'Color',[0,0.447,0.741]); hold on;
plot(toolCpts(2,:),toolCpts(3,:),'x','Color',[0.32,0.55,0.19],'MarkerSize',5);
plot(toolEdgePt(2,:),toolEdgePt(3,:),'Color',[0.635,0.078,0.184]);
axis equal
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['y(',unit,')']);
ylabel(['z(',unit,')']);
title('Tool interpolation results')
legend('Measured Pts','Control Pts','Fitting Pts','Location','best');

% h4 = axes('Position',[0.3,0,0.4,0.4]);
% plot(toolFit(2,:),toolFit(3,:),'--.', ...
%     'MarkerSize',8,'Color',[0,0.447,0.741]); hold on;
% plot(toolCpts(2,:),toolCpts(3,:),'x','Color',[0.32,0.55,0.19],'MarkerSize',5);
% plot(toolEdgePt(2,:),toolEdgePt(3,:),'Color',[0.635,0.078,0.184]);
% set(h4,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[],'XColor','w','YColor','w');
% annotation('rectangle',[0.51+(l1MaxPrec-0.00025)/2/l1Prec(end),(l2MaxPrec-0.00025)/l2Prec(end),...
%     0.006/l1Prec(end),0.01/l2Prec(end)],'LineStyle','-','Color','w','LineWidth',1);


% plot the interpolation error
% figure('Name','Tool interpolation Error');

%% save the tool interpolation results
center = [0;0;0];
toolEdgeNorm = [0;0;1];
cutDirect = [1;0;0];
toolDirect = [0;1;0];
[toolFileName,toolDirName,toolFileType] = uiputfile({ ...
        '*.mat','MAT-file(*.mat)'; ...
        '*.txt','text-file(.txt)';...
        '*.*','all file(*.*)';...
        }, ...
        'Select the directory and filename to save the tool model', ...
        fullfile(workspaceDir,['toolTheo',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
toolFile = fullfile(toolDirName,toolFileName);
switch toolFileType
    case 0 % no saving files
        msgfig = msgbox("No tool model saved","Warning","warn","modal");
        uiwait(msgfig);
    case 1 % *.mat
        Comments = cell2mat(inputdlg( ...
            'Enter Comment of the tool model:', ...
            'Saving Comments', ...
            [5 60], ...
            string(datestr(now))));
        save(toolFile,"Comments","unit","FitMethod","ParamMethod", ... % comments and notes
            "center","radius","openAngle", ... % tool fitting results
            "toolEdgeNorm","toolDirect","cutDirect","toolBform", ... % tool interpolation results
            "toolEdgePt","toolFit"); % auxiliary data
        % toolEdgePt, toolCpts, toolFit are useless in the following process at present
    otherwise
        msgfig = msgbox("File type error","Error","error","modal");
        uiwait(msgfig);
end

%%
% rmpath(genpath('funcs'));