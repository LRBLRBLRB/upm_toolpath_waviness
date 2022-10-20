% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高
if true
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
    fitMethod = 'Levenberg-Marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
end

%% simulation initialization
debug = 2;
switch debug
    case 2
        cx0 = 0; % unit: nm
        cy0 = 0;
        r0 = 0.1;
        noise = 0.05;
        theta = transpose(linspace(0,2*pi/3,50));
        r = r0*(1 - noise + 2*noise*rand(length(theta),1));
        toolOri(:,1) = r.*cos(theta) + cx0;
        toolOri(:,2) = r.*sin(theta) + cy0;
        % toolOri(:,3) = 0;
        rmse0 = sqrt(...
                    sum(...
                        ((toolOri(:,1) - cx0).^2 + (toolOri(:,2) - cy0).^2 - r0^2).^2) ...
                /length(toolOri));
        clear theta r;
    otherwise
        toolInput = inputdlg({'Tool fitting method','Tool parameterization method', ...
            'Unit','Font type in figure','Font size in figure'}, ...
            'Tool Processing Input', ...
            [1 20; 1 20; 1 20; 1 20; 1 20], ...
            {'Levenberg-Marquardt','centripetal','\mu m','Times New Roman','14'}, ...
            'WindowStyle');
        fitMethod = toolInput{1};
        paramMethod = toolInput{2};
        unit = toolInput{3};
        textFontType = toolInput{4};
        textFontSize = str2double(toolInput{5});
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

% dim = size(toolOri,1);
f1 = figure('Name','original xy scatters of the tool');
plot(toolOri(:,1),toolOri(:,2));
hold on;

theta1 = atan2(toolOri(1,2),toolOri(1,1));
theta2 = atan2(toolOri(end,2),toolOri(end,1));
theta = transpose(linspace(theta1,theta2,100));
plot(cx0 + r0*cos(theta),cy0 + r0*sin(theta),'blue');
scatter(cx0,cy0,'blue');
axis equal;
xlabel('x(nm)');
ylabel('y(nm)');
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
[center,radius,includedAngle,toolXY,rmseLsc] = toolFit2D(toolOri,'Levenberg-Marquardt');
f2 = figure('Name','tool sharpness fitting result');
plot(toolXY(:,1),toolXY(:,2));
hold on;
theta = (pi/2 - includedAngle/2):0.01:(pi/2 + includedAngle/2);
xtmp = radius*cos(theta);
ytmp = radius*sin(theta);
plot(xtmp,ytmp,'black');
scatter(0,0,'black');
axis equal;
xlabel('x(nm)');
ylabel('y(nm)');
clear xtmp ytmp theta; % 删除画图的临时变量

% 直角坐标系转极坐标系
toolTheta = atan2(toolXY(:,2),toolXY(:,1));
toolR = vecnorm(toolXY,2,2);

%% 车刀轮廓插值

% 轮廓插值的测试输入
% center = [0,0]; % um
% radius = 4; % um
% toolTheta = (pi/6:0.01:5*pi/6)';
% toolR = -0.2 + 0.1 * rand(length(toolTheta),1) + radius; % 实际上是由报告获得，单位：um

figure('name','tool curve');
subplot(2,1,1); plot(toolTheta*180/pi,toolR); title('tool curve');
subplot(2,1,2); plot(toolTheta*180/pi,toolR - radius); title('tool waviness');

[toolCpts,U] = bSplineCpts(toolXY,3,'chord'); 
u = 0:0.002:1;
n = length(toolCpts);
toolPt = zeros(length(u),2);
for i = 1:length(u)
    toolPt(i,:) = bSplinePt(toolCpts,3,u(i),U);
end

% plot the interpolation results
figure('Name','Tool Interpolation Results');
plot(toolXY(:,1),toolXY(:,2),'Color',[0,0.45,0.74]); hold on;
plot(toolCpts(:,1),toolCpts(:,2),'x','Color',[0.85,0.33,0.10]);
axis equal
legend('Measured Pts','Interpolation Cpts','Location','best');

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
        save(toolFile,"Comments","unit","fitMethod","paramMethod", ... % comments and notes
            "center","radius","openAngle", ... % tool fitting results
            "toolEdgeNorm","toolDirect","cutDirect","toolBform", ... % tool interpolation results
            "toolEdgePt","toolFit"); % auxiliary data
        % toolEdgePt, toolCpts, toolFit are useless in the following process at present
    otherwise
        msgfig = msgbox("File type error","Error","error","modal");
        uiwait(msgfig);
end

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