% 给定刀具波纹度拟合曲线和已加工曲面拟合曲线，给定残高，求各个刀位点的切宽
% 方案：先求切宽和残高的定量关系；然后考虑如何移动各个刀位点，微调来获得所要求的残高
close all;
clear; clc;
addpath(genpath('funcs'));

%% simulation initialization
debug = 3;
switch debug
    case 2 % 2D tool profile simulation
        cx0 = 0; % unit: nm
        cy0 = 0;
        r0 = 0.1;
        noise = 0.05;
        zNoise = 0.01;
        theta = transpose(linspace(0,2*pi/3,100)); % need debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        r = r0*(1 - noise + 2*noise*rand(length(theta),1));
        toolOri(:,1) = r.*cos(theta) + cx0;
        toolOri(:,2) = r.*sin(theta) + cy0;
        toolOri(:,3) = zeros(length(theta),1);
        rmse0 = sqrt(...
                    sum(...
                        ((toolOri(:,1) - cx0).^2 + (toolOri(:,2) - cy0).^2 - r0^2).^2) ...
                /length(theta));
        clear theta r;
    case 3 % 3D tool profile simulation
        cx0 = 1;
        cy0 = 2;
        cz0 = 3;
        r0 = 0.1; % unit:um
        noise = 0.05;
        zNoise = 0.01;
        theta = transpose(linspace(0,2*pi/3,100));
        r = r0*(1 - noise + 2*noise*rand(length(theta),1));
        toolOri(:,1) = cx0 + r.*cos(theta);
        toolOri(:,2) = cy0 + r.*sin(theta);
        toolOri(:,3) = cz0 + r0*(1 - zNoise + 2*zNoise*rand(length(theta),1));
        rmse0 = sqrt( ...
            sum((toolOri - meshgrid([cx0,cy0,cz0],1:length(theta)).^2),2) ...
            /length(theta));
        clear theta r;
    otherwise
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
plot3(toolOri(:,1),toolOri(:,2),toolOri(:,3));
hold on;
% axis equal;
xlabel('x(nm)');
ylabel('y(nm)');
clear theta theta1 theta2;

%% 由(x,y)图获得刃口圆心半径以及波纹度函数
[center,radius,includedAngle,toolFit,rmseLsc] = toolFit3D(toolOri,'Levenberg-Marquardt');
f2 = figure('Name','tool sharpness fitting result');
plot(toolFit(:,1),toolFit(:,2));
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

% Cartesian coordinate to cylindrical coordinate？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
toolTheta = atan2(toolFit(:,2),toolFit(:,1));
toolR = vecnorm(toolFit,2,2);

%% 车刀轮廓插值
figure('name','tool curve');
subplot(2,1,1); plot(toolTheta*180/pi,toolR); title('tool curve');
subplot(2,1,2); plot(toolTheta*180/pi,toolR - radius); title('tool waviness');

k = 3; % order of the B-spline
[toolCpts,U] = bSplineCpts(toolFit,k,'chord'); 
u = 0:0.002:1;
nPts = length(u);
toolDim = size(toolCpts,2);
toolPt = zeros(nPts,toolDim);
for i = 1:nPts
    toolPt(i,:) = bSplinePt(toolCpts,k,u(i),U);
end

% plot the interpolation results
figure('Name','Tool Interpolation Results');
plot(toolFit(:,1),toolFit(:,2),'--.', ...
    'MarkerSize',8,'Color',[0.32,0.55,0.19]); hold on;
plot(toolCpts(:,1),toolCpts(:,2),'x','Color',[0.85,0.33,0.10]);
plot(toolPt(:,1),toolPt(:,2),'Color',[0,0.45,0.74]);
axis equal
legend('Measured Pts','Control Pts','Fitting Pts','Location','best');

%% save the tool interpolation results
toolEdgeNorm = [0,1,0];
center(1,3) = 0;
[toolFileName,toolDirName,toolFileType] = uiputfile({ ...
        '*.mat','MAT-file(*.mat)'; ...
        '*.txt','text-file(.txt)';...
        '*.*','all file(*.*)';...
        }, ...
        'Select the directory and filename to save the tool model', ...
        ['Tool/toolTheo',datestr(now,'yyyymmddTHHMMSS'),'.mat']);
toolFile = fullfile(toolDirName,toolFileName);
switch toolFileType
    case 0
        msgfig = msgbox("No tool model saved","Warning","warn","modal");
        uiwait(msgfig);
    case 1
        Comments = cell2mat(inputdlg('Enter Comment of the tool model:', ...
            'Input Saving Comments',[5 60],string(datestr(now))));
        save(toolFile,"center","radius","Comments","includedAngle", ...
            "toolPt","toolEdgeNorm");
    otherwise
        msgfig = msgbox("File type error","Error","error","modal");
        uiwait(msgfig);
end

%% 相邻刀位点的残高
x1 = [0,0]; vec1 = [0,1];
x2 = [r0,r0/30]; vec2 = [cosd(90+8),sind(90+8)];
R = rotz(5,'deg'); % rotation matrix about z axis
tool1 = toolPt + x1;
tool2 = toolPt*(R(1:2,1:2)') + x2;

[res,interPt] = residualHigh(x1(:,1:2),vec1,tool1(:,1:2),x2(:,1:2),vec2,tool2(:,1:2));

figure('Name','Residual of the adjacent tool');
plot(tool1(:,1),tool1(:,2),'Color',[0,0.45,0.74]); hold on;
plot(tool2(:,1),tool2(:,2),'Color',[0.85,0.33,0.10]);
plot(interPt(1),interPt(2),'*','Color',[0.49,0.18,0.56]);
plot(x1(1),x1(2),'o','Color',[0,0.45,0.74],'MarkerFaceColor',[0,0.45,0.74]);
line([x1(1),x1(1)+radius*vec1(1)],[x1(2),x1(2)+radius*vec1(2)],'LineWidth',2,'Color',[0,0.45,0.74]);
plot(x2(1),x2(2),'o','Color',[0.85,0.33,0.10],'MarkerFaceColor',[0.85,0.33,0.10]);
line([x2(1),x2(1)+radius*vec2(1)],[x2(2),x2(2)+radius*vec2(2)],'LineWidth',2,'Color',[0.85,0.33,0.10]);
axis equal

rmpath(genpath('funcs'));