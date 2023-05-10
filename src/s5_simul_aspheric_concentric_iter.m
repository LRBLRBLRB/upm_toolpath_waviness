% simulation of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path calculation
% Step three: residual height calculation and correction
% Step four: simulation of the machining surface

isAPP = true;
if isAPP
    workspaceDir = app.workspaceDir;
    unit = app.unit;
    textFontSize = app.fontSize;
    textFontType = app.fontName;

    % machining paramters
    aimRes = app.aimRes;
    toolData = app.toolData;
    rStep = toolData.radius; % 每步步长可通过曲面轴向偏导数确定
    rMax = app.rMax;
    arcLength = app.arcLength;
    maxAngPtDist = app.maxAngPtDist;
    angularLength = app.angularLength;

    surfFunc = app.surfFuncs;
    surfFx = app.surfFx;
    surfFy = app.surfFy;
    surfDomain = app.surfDomain;
    surfMesh = app.surfMesh;
    
    msgOpts.Default = 'Cancel and quit';
    msgOpts.Interpreter = 'tex';
    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);
else
    close all;
    clear;
    clc;
    addpath(genpath('funcs'));
    t0 = tic;
    
    workspaceDir = fullfile('..','workspace','/20220925-contrast/nagayama_concentric';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    
    msgOpts.Default = 'Cancel and quit';
    msgOpts.Interpreter = 'tex';
    % msgOpts.modal = 'non-modal';
    profile on
    tPar0 = tic;
    parObj = gcp;
    tPar = toc(tPar0);
    fprintf('The time spent in the parallel computing activating process is %fs.\n',tPar);

    %% concentric surface generation / import
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select one tool edge data file', ...
        fullfile(workspaceDir,'tooltheo.mat'), ...
        'MultiSelect','off');
    toolName = fullfile(dirName,fileName);
    % toolName = 'output_data\tool\toolTheo_3D.mat';
    toolData = load(toolName);
    
    default = false;
    if default
        [fileName,dirName] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','all files(*.*)'}, ...
            'Select the surface edge data file', ...
            fullfile('..','workspace','\input_data\surface\ellipsoidAray.mat', ...
            'MultiSelect','off');
        surfName = fullfile(dirName,fileName);
        load(surfName);
    % % % % % % % % % % % % % % % % % % % % % % % % %     sparTheta
    % % % % % % % % % % % % % % % % % % % % % % % % %     surfMesh
    % % % % % % % % % % % % % % % % % % % % % % % % %     surfNorm
    else % ellipsoid
        R = 10/2*1000;
        A = 3.5/2;
        B = 4/2;
        C = 5/2;
        % machining surface
        syms x y;
        surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
        surfFunc = matlabFunction(surfSym);
        surfFx = matlabFunction(diff(surfFunc,x));
        surfFy = matlabFunction(diff(surfFunc,y));
        rMax = R/2;
        % sampling density
        r = [0,R/4]; % concentric radius range
        sparTheta = 101;
        surfCenter = [0,0,sqrt(C^2*(R.^2-r(2).^2))]; % concentric circle center
        conR = (toolData.radius/8):(toolData.radius/2):r(2); % concentric radius vector
        densR = length(conR);
        conTheta = linspace(0,2*pi,sparTheta);
        % conTheta(end) = [];
    
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = sqrt(C^2*(R.^2-rMesh.^2));
    %  y = reshape(surfMesh(:,:,2),[],1);
    % z = reshape(surfMesh(:,:,3),[],1);
    % surfXYZ = [x,y,z];
    
    % calculate the normal vector of the analytic surface
    % syms X Y z (X,Y);
    % z (X,Y) = sqrt(C^2*(R^2 - X.^2/A^2 - Y.^2/B^2));
    % ZDX = diff(Z,X);
    % ZDY = diff(Z,Y);
    % surfNorm(:,1) = eval(subs(ZDX,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,2) = eval(subs(ZDY,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,3) = -ones(densTheta*densR,1);
    [surfNorm(:,:,1),surfNorm(:,:,2),surfNorm(:,:,3)] = surfnorm( ...
        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));
%     save(fullfile('..','workspace','/input_data/surface/ellipsoidAray.mat', ...
%         "surfMesh","surfNorm","surfCenter");
    end
end

sparMode = 'angle';
if strcmp(sparMode,'angle')
    surfPt = transpose(reshape(surfMesh,[],3));
    surfNorm = transpose(reshape(surfNorm,[],3));
    ptNum = size(surfPt,2);
    surfDirect = cutdirection(surfPt,[0;0;0],'method','concentric');
    % surfDirect = cutdirection(surfPt,'method','vertical');
    % surfDirect = zeros(3,ptNum);
    % for ii = 1:densR
    %     surfDirect(:,(ii - 1)*sparTheta + 1:ii*sparTheta) = ...
    %         cutDirection(surfPt(:,(ii - 1)*sparTheta + 1:ii*sparTheta));
    % end
else % still need debugging!
    syms x y;
    surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
    surfFunc = matlabFunction(surfSym);
    surfFx = matlabFunction(diff(surfFunc,x));
    surfFy = matlabFunction(diff(surfFunc,y));
    load("..\workspace\20220925-contrast\nagayama_concentric\loopR.mat");
    % r = loopR;
    r = (toolData.radius/2):(toolData.radius/2):R/2;
    arcLength = 30;
    maxAngPtDist = 6*pi/180;

    surfX = [];
    surfY = [];
    surfNorm = [];
    for ii = 1:length(r)
        conTheta = linspace(0,2*pi, ...
            ceil(2*pi/min(maxAngPtDist,arcLength/r(ii))) + 1);
        conTheta(end) = [];
        sparTheta(ii) = length(conTheta);
        surfX = [surfX,r(ii)*cos(conTheta)];
        surfY = [surfY,r(ii)*sin(conTheta)];
    end
    ptNum = length(surfX);
    surfPt(1:2,:) = [surfX;surfY];
    surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:));
    surfNorm(1,:) = surfFx(surfPt(1,:),surfPt(2,:));
    surfNorm(2,:) = surfFy(surfPt(1,:),surfPt(2,:));
    surfNorm(3,:) = -1*ones(1,ptNum);
    surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1));
%     surfDirect = cutdirection(surfPt,[0;0;0],'method','concentric');
    surfDirect = cutdirection(surfPt,'method','vertical');
end

%% plot the aspheric surface
figure('Name','original xyz scatters of the surface (sparsely)');
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
tiledlayout(1,2);
nexttile;
plot3(surfPt(1,:),surfPt(2,:),surfPt(3,:),'.','Color',[0,0.45,0.74]);
hold on;
quiver3(surfPt(1,:),surfPt(2,:),surfPt(3,:), ...
    surfNorm(1,:),surfNorm(2,:),surfNorm(3,:), ...
    'AutoScale','on','Color',[0.85,0.33,0.10],'DisplayName','Normal Vectors');
quiver3(surfPt(1,:),surfPt(2,:),surfPt(3,:), ...
    surfDirect(1,:),surfDirect(2,:),surfDirect(3,:), ...
    'AutoScale','on','Color',[0.4660,0.6740,0.1880],'DisplayName','Normal Vectors');
legend('Original Points','Orthogonal direction','Cutting direction','Location','northeast');
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
axis equal; grid on;
% title({'Radially & circunferentially sparse','by 50 and 2 times, respectively'}, ...
%     'FontSize',textFontSize,'FontName',textFontType);
nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on; axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% cb2 = colorbar;
msgfig = questdlg({'Surface was generated successfully!', ...
    'Ready to continue?'}, ...
    'Surface Generation','OK & continue','Cancel & quit','OK & continue');
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    msgbox('Exit for the program','Exit','help','modal');
    uiwait(msgbox);
    profile off;
    tTol = toc(t0);
    fprintf("The time spent in the whole process is %fs.\n",tTol);
    return;
end

%% Calculation of Tool Path & Spindle Drection
toolQuat = zeros(ptNum,4);
toolVec = zeros(3,ptNum);
toolPathPt = zeros(3,ptNum);
toolCutDirect = zeros(3,ptNum);
toolNormDirect = zeros(3,ptNum);
toolContactU = zeros(1,ptNum);
isCollision = false(1,ptNum);
tToolpath0 = tic;
% figure;
parfor ii = 1:ptNum
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),[0;0;1],surfDirect(:,ii));
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
%     plot3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii)); hold on;
%     quiver3(toolPathPt(1,ii),toolPathPt(2,ii),toolPathPt(3,ii), ...
%         toolCutDirect(1,ii),toolCutDirect(2,ii),toolCutDirect(3,ii), ...
%         'AutoScale','on');
%     toolSp1 = toolData.toolBform;
%     toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolData.toolBform.coefs + toolVec(:,ii);
%     Q = fnval(toolSp1,0:0.01:1);
%     plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
%     xlabel('x'); ylabel('y'); axis equal; hold on; 
end
tToolpath = toc(tToolpath0);
fprintf('The time spent in the tool path generating process is %fs.\n',tToolpath);

% still need procedures to deal with the invalid points, which will cause
% interference between the tool edge and the designed surface. 

figure('Name','tool center position & tool normal vector');
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','Color',[0.6350,0.0780,0.1840]);
hold on;
quiver3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    toolCutDirect(1,1:plotSpar:end), ...
    toolCutDirect(2,1:plotSpar:end), ...
    toolCutDirect(3,1:plotSpar:end), ...
    'AutoScale','on','Color',[0,0.4470,0.7410]);
quiver3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    toolNormDirect(1,1:plotSpar:end), ...
    toolNormDirect(2,1:plotSpar:end), ...
    toolNormDirect(3,1:plotSpar:end), ...
    'AutoScale','on','Color',[0.85,0.33,0.10]);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.6,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = 'Height (mm)';
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('tool center point','tool cutting direction', ...
    'tool spindle direction','','Location','northeast');
msgfig = questdlg({'Tool Path was calculated successfully!', ...
    'Ready to continue?'}, ...
    'Tool Path Simulation','OK & continue','Cancel & quit','OK & continue');
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    msgfig = msgbox('Exit for the program','Exit','help','modal');
    uiwait(msgfig);
    profile off;
    tTol = toc(t0);
    fprintf("The time spent in the whole process is %fs.\n",tTol);
    return;
end



%% Calculation of Residual Height & Cutting Surface
toolSp = toolData.toolBform;
toolRadius = toolData.radius;
resNum = ptNum - sparTheta;
res = zeros(1,2*resNum);
peakPt = zeros(3,2*resNum);
uLim = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of each toolpath
% uLimTmp = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of the projective toolpath
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));

tTip0 = tic;
% outer side of each point in the tool path
parfor ii = 1:resNum
    % 如果是沿同一个极径的，就可以直接不用投影；否则还是需要这样子找
    nLoop = floor((ii - 1)/sparTheta) + 1;
    angleN = angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1));
    angleDel = angleN - angle(ii);
    % ind2(ii) remains the index of angle nearest to angle(ii) within 
    % those which is larger than the angle(ii) and in angleN
    if isempty(angleN(angleDel >= 0))
        % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
        angleDel = angleDel + 2*pi;
    end
    ind2 = sparTheta*nLoop + find(angleN == min(angleN(angleDel >= 0)));
    % ind3(ii) remains the index of angle nearest to angle(ii) within 
    % those which is smaller than the angle(ii) and in angleN
    if isempty(angleN(angleDel < 0))
        angleDel = angleDel - 2*pi;
    end
    ind3 = sparTheta*nLoop + find(angleN == max(angleN(angleDel < 0)));
    [res(ii),peakPt(:,ii),uLim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp,toolRadius, ...
        uLim(:,ii),ii,ind2,ind3);
end
% inner side of each point on the tool path
parfor ii = (sparTheta + 1):ptNum
    % 如果是沿同一个极径的，就可以直接不用投影；否则还是需要这样子找
    nLoop = floor((ii - 1)/sparTheta) - 1;
    angleN = angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1));
    angleDel = angleN - angle(ii);
    % ind2(ii) remains the index of angle nearest to angle(ii) within 
    % those which is larger than the angle(ii) and in angleN
    if isempty(angleN(angleDel >= 0))
        % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
        angleDel = angleDel + 2*pi;
    end
    ind2 = sparTheta*nLoop + find(angleN == min(angleN(angleDel >= 0)));
    % ind3(ii) remains the index of angle nearest to angle(ii) within 
    % those which is smaller than the angle(ii) and in angleN
    if isempty(angleN(angleDel < 0))
        angleDel = angleDel - 2*pi;
    end
    ind3 = sparTheta*nLoop + find(angleN == max(angleN(angleDel < 0)));
    [res(resNum + ii),peakPt(:,resNum + ii),uLim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp,toolRadius, ...
        uLim(:,ii),ii,ind2,ind3);
end
tTip = toc(tTip0);
fprintf('The time spent in the tool tip calculation process is %fs.\n',tTip);

res(resNum + 1:ptNum) = [];
peakPt(:,resNum + 1:ptNum) = [];

% post-processing of the residual height data
% ( this process is in the visualization part now)

% for ii = 1:ptNum
%     uLim(1,ii) = max([uLim(1,ii),uLimTmp(1,ind2 == ii)]);
%     uLim(2,ii) = min([uLim(2,ii),uLimTmp(2,ind2 == ii)]);
% end
clear angle ind2 ind3 uLimTmp;

figure('Name','residual height calculation');
plotSpar = 1;
tPlot0 = tic;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
hold on;
% quiver3(toolPathPt(1,1:plotSpar:end), ...
%     toolPathPt(2,1:plotSpar:end), ...
%     toolPathPt(3,1:plotSpar:end), ...
%     toolCutDirect(1,1:plotSpar:end), ...
%     toolCutDirect(2,1:plotSpar:end), ...
%     toolCutDirect(3,1:plotSpar:end), ...
%     'AutoScale','on','Color',[0.6350,0.0780,0.1840]);
% quiver3(toolPathPt(1,1:plotSpar:end), ...
%     toolPathPt(2,1:plotSpar:end), ...
%     toolPathPt(3,1:plotSpar:end), ...
%     toolNormDirect(1,1:plotSpar:end), ...
%     toolNormDirect(2,1:plotSpar:end), ...
%     toolNormDirect(3,1:plotSpar:end), ...
%     'AutoScale','on','Color',[0.85,0.33,0.10]);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
parfor ii = 1:ptNum
    toolSp = toolData.toolBform;
    toolCoefs = toolData.toolBform.coefs;
    toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp,uLim(1,ii):0.01:uLim(2,ii));
    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
end
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
legend('tool center point','','tool edge','Location','northeast');
tPlot = toc(tPlot0);
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);


%% Visualization & Simulation
nLoop = ceil(ptNum/sparTheta);
s2_visualize_process;

% delete(parObj);
profile off
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));