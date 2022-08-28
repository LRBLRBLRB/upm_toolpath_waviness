% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高
close all;
clear; clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType unit fitMethod paramMethod;
FitMethod = 'Levenberg-Marquardt';
ParamMethod = 'centripetal';
unit = 'mm';
textFontSize = 14;
textFontType = 'Times New Roman';

%% simulation initialization
debug = 3;
switch debug
    case 2 % 2D tool profile simulation
        cx0 = 0; % unit: mm
        cy0 = 0; % unit: mm
        r0 = 0.1; % unit: mm
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
        cx0 = 1; % unit:mm
        cy0 = 2; % unit:mm
        cz0 = 3; % unit:mm
        r0 = 0.1; % unit:mm
        noise = 0.05;
        zNoise = 0.01;
        theta = linspace(0,2*pi/3,300);
        r = r0*(1 - noise + 2*noise*rand(1,length(theta)));
        toolOri(1,:) = cx0 + r.*cos(theta);
        toolOri(2,:) = cy0 + r.*sin(theta);
        toolOri(3,:) = cz0 + r0*(1 - zNoise + 2*zNoise*rand(1,length(theta)));
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

f1 = figure('Name','original scatters of the tool');
scatter3(toolOri(1,:),toolOri(2,:),toolOri(3,:));
hold on; grid on;
% axis equal;
xlabel(['x(',unit,'m)'],'FontSize',textFontSize,'FontName',textFontType);
ylabel(['y(',unit,'m)'],'FontSize',textFontSize,'FontName',textFontType);
zlabel(['z(',unit,'m)'],'FontSize',textFontSize,'FontName',textFontType);
view(-45,60);
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
nCPts = size(toolOri,2);
[~,radius,includedAngle,toolFit,rmseLsc] = toolFit3D(toolOri,FitMethod);
% plot the fitting results
f2 = figure('Name','tool sharpness fitting result');
xLim = 1.1*max(toolFit(1,:));
quiver(-xLim,0,2*xLim,0,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % X axis
hold on;
text(0.9*xLim,-.05*radius,'x','FontSize',textFontSize,'FontName',textFontType);
quiver(0,-0.2*radius,0,1.3*radius,'AutoScale','off','Color',[0,0,0],'MaxHeadSize',0.1); % Y axis
text(0.05*xLim,1.2*radius,'y','FontSize',textFontSize,'FontName',textFontType);
plot(toolFit(1,:),toolFit(2,:),'Color',[0,0.45,0.74]); % tool edge scatters
theta = (pi/2 - includedAngle/2):0.01:(pi/2 + includedAngle/2);
xtmp = radius*cos(theta);
ytmp = radius*sin(theta);
plot(xtmp,ytmp,'Color',[0.85,0.33,0.10]); % tool edge circle
scatter(0,0,'MarkerFaceColor',[0.85,0.33,0.10],'MarkerEdgeColor',[0.85,0.33,0.10]); % tool edge center
quiver(0,0,0,0.5*radius,'AutoScale','off','Color',[0.93,0.69,0.13], ...
    'LineWidth',2.5,'MaxHeadSize',0.3); % tool edge normal
line([0,xtmp(1)],[0,ytmp(1)],'LineStyle','--','Color',[0.85,0.33,0.10]);
line([0,xtmp(end)],[0,ytmp(end)],'LineStyle','--','Color',[0.85,0.33,0.10]);
% xlim([-1.1*xLim,1.1*xLim]);
axis equal;
% set(gca,'TickLabelInterpreter','tex');
xlabel(['x(',unit,'m)'],'FontSize',textFontSize,'FontName',textFontType);
ylabel(['y(',unit,'m)'],'FontSize',textFontSize,'FontName',textFontType);
legend('','','tool edge','tool fitting arc','tool center', ...
    'tool normal vector','Location','northeast');
clear xtmp ytmp theta xLim; % 删除画图的临时变量

% Cartesian coordinate to cylindrical coordinate
toolTheta = atan2(toolFit(2,:),toolFit(1,:));
toolR = vecnorm(toolFit,2,1);
figure('name','tool curve');
tiledlayout(2,1);
nexttile;
plot(toolTheta*180/pi,toolR); hold on;
title('tool curve','FontSize',textFontSize,'FontName',textFontType);
ylabel({'polar diameter',['(',unit,'m)']},'FontSize',textFontSize,'FontName',textFontType);
nexttile;
plot(toolTheta*180/pi,toolR - radius); hold on;
title('tool waviness','FontSize',textFontSize,'FontName',textFontType);
ylabel({'polar diameter error',['(',unit,'m)']},'FontSize',textFontSize,'FontName',textFontType);
xlabel('central angle\theta(°)','FontSize',textFontSize,'FontName',textFontType);

%% 车刀轮廓插值 two methods to interpolate
k = 3; % degree of the B-spline
u = 0:0.0002:1;
nPts = length(u);

% B-spline interpolate in the polar coordinate
[toolPt,toolBform] = bsplinePts_spapi(toolFit,k,u, ...
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
plot(toolFit(1,:),toolFit(2,:),'--.', ...
    'MarkerSize',8,'Color',[0,0.447,0.741]); hold on;
plot(toolCpts(1,:),toolCpts(2,:),'x','Color',[0.32,0.55,0.19],'MarkerSize',5);
plot(toolPt(1,:),toolPt(2,:),'Color',[0.635,0.078,0.184]);
axis equal
legend('Measured Pts','Control Pts','Fitting Pts','Location','best');

% plot the interpolation error
% figure('Name','Tool interpolation Error');

%% save the tool interpolation results
center = [0;0;0];
toolFit = [toolFit(1,:);zeros(1,nCPts);toolFit(2,:)];
toolBform.coefs = [toolBform.coefs(1,:);zeros(1,nCPts);toolBform.coefs(2,:)];
toolPt = [toolPt(1,:);zeros(1,nPts);toolPt(2,:)]; % get
toolEdgeNorm = [0;0;1];
cutDirect = [0;-1;0]; % Caution!!!
toolDirect = [1;0;0];
[toolFileName,toolDirName,toolFileType] = uiputfile({ ...
        '*.mat','MAT-file(*.mat)'; ...
        '*.txt','text-file(.txt)';...
        '*.*','all file(*.*)';...
        }, ...
        'Select the directory and filename to save the tool model', ...
        ['output_data/tool/toolTheo',datestr(now,'yyyymmddTHHMMSS'),'.mat']);
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
            "center","radius","includedAngle", ... % tool fitting results
            "toolEdgeNorm","toolDirect","cutDirect","toolBform", ... % tool interpolation results
            "toolPt","toolFit"); % auxiliary data
        % toolPt, toolCpts, toolFit are useless in the following process at present
    otherwise
        msgfig = msgbox("File type error","Error","error","modal");
        uiwait(msgfig);
end

%%
% rmpath(genpath('funcs'));