% to continue the tool modelling process, this file aims to interpolate the
% tool tip arc based on the data that is extracted and analyzed in the file
% "s1_tool3D.m" or "s1_toolExtract_line.m"

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
theta = (pi/2 - openAngle/2):0.01:(pi/2 + openAngle/2);
xtmp = radius*cos(theta);
ytmp = radius*sin(theta);
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

nCPts = size(toolFit,2);
toolFit = [zeros(1,nCPts);toolFit];
% B-spline interpolate in the polar coordinate
[toolEdgePt,toolBform] = bsplinePts_spapi(toolFit,k,u, ...
    'paramMethod',paramMethod,'cptsType','Cartesian'); 
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
switch toolFileType
    case 0 % no saving files
        toolDataFile = '';
        msgfig = msgbox("No tool model saved","Warning","warn","modal");
        uiwait(msgfig);
    case 1 % *.mat
        toolDataFile = fullfile(toolDirName,toolFileName);
        Comments = cell2mat(inputdlg( ...
            'Enter Comment of the tool model:', ...
            'Saving Comments', ...
            [5 60], ...
            string(datestr(now))));
        save(toolDataFile,"Comments","unit","fitOpts","paramMethod", ... % comments and notes
            "center","radius","openAngle", ... % tool fitting results
            "toolEdgeNorm","toolDirect","cutDirect","toolBform", ... % tool interpolation results
            "toolEdgePt","toolFit"); % auxiliary data
        % toolEdgePt, toolCpts, toolFit are useless in the following process at present
    otherwise
        toolDataFile = fullfile(toolDirName,toolFileName);
        msgfig = msgbox("File type error","Error","error","modal");
        uiwait(msgfig);
end

%%
% rmpath(genpath('funcs'));