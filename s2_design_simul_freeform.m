% simulation of the designed tool path of a freeform surface
% Step one: tool and surface data import
% Step two: tool path calculation
% Step three: residual height calculation and correction
% Step four: simulation of the machining surface
close all;
clear;
clc;
addpath(genpath('funcs'));
t0 = tic;

% global variables
% global textFontSize textFontType;
workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
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

%% to load the tool data
[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)'; ...
    '*,*','all files(*.*)'}, ...
    'Select one tool edge data file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');
toolName = fullfile(dirName,fileName);
% toolName = 'output_data\tool\toolTheo_3D.mat';
toolData = load(toolName);

%% to load and reconstruct the freeform surface
% aim: to get the function of the surface ,and its partial derivative
xlim = 1000*input('x range (unit: mm): \n');
if ~isequal(size(xlim),[1,2])
    error('invalid x range: the dimension is wrong');
end
ylim = 1000*input('y range (unit: mm): \n');
if ~isequal(size(ylim),[1,2])
    error('invalid y range: the dimension is wrong');
end
surfType = 'function';
if strcmp(surfType,'function')
    syms x y;
    surfSym = input("Type the function expression of the surface," + ...
        "with the independent variables of \'x\' and \'y\': \n");
    if isempty(surfSym) % default: biconic surface
        rx = 1000*100; % radius in x direction
        cx = 1/rx; % curvature in x direction
        ry = 1000*100; % radius in y direction
        cy = 1/ry; % curvature in y direction
        kx = -1; % conic in x direction
        ky = -1; % conic in y direction
        Ar = 0; % polynomials (same below)
        Ap = 0;
        Br = 0;
        Bp = 0;
        Cr = 0;
        Cp = 0;
        Dr = 0;
        Dp = 0;
        surfSym = (cx*x^2 + cy*y^2) ...
            /(1 + sqrt(1 - (1 + kx)*cx^2*x^2 - (1 + ky)*cy^2*y^2)) ...
            + Ar*((1 - Ap)*x^2 + (1 + Ap)*y^2)^2 ...
            + Br*((1 - Bp)*x^2 + (1 + Bp)*y^2)^3 ...
            + Cr*((1 - Cp)*x^2 + (1 + Cp)*y^2)^4 ...
            + Dr*((1 - Dp)*x^2 + (1 + Dp)*y^2)^5;
    elseif ~(contains(surfSym,'x') && contains(surfSym,'y'))
        error('invalid function expression: independent variables must be x & y.');
    end

    surfFunc = matlabFunction(surfSym);
    surfFx = matlabFunction(diff(surfFunc,x));
    surfFy = matlabFunction(diff(surfFunc,y));
        
    % discretization of the freeform surface (constant arc)
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
    surfDirect = cutDirection(surfPt,[0;0;0]);
    clear surfX surfY;

      
else % point cloud input
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select the surface edge data file', ...
        'workspace\input_data\surface\ellipsoidAray.mat', ...
        'MultiSelect','off');
    surfName = fullfile(dirName,fileName);
    surfData = load(surfName);
    dim = size(surfData.xyz,1);

    f1 = figure('Name','original xyz scatters of the surface');
    plot3(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3),'.');
    hold on; axis equal; grid on;
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % B-spline interpolation of the freeform surface in polar coordinate
    xPlot = linspace(xlim(1),xlim(2));
    yPlot = linspace(ylim(1),ylim(2));
    [surfMesh(:,:,1),surfMesh(:,:,2)] = meshgrid(); % 这里的meshgrid，x和y的顺序是不是我想要的？
    griddata(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3))
    % [surfCpts,U,V] = bSplineSurfCpts(surfData.mesh,3,3);
    u = 0:0.002:1;
    v = 0:0.002:1;
    % surfPts = zeros(length(u),length(v),dim);
    % surfPts = bSplineSurfPts(surfCpts,3,u,U,3,v,V);
% 这里的Q，是用极坐标的还是直角坐标的？
    [surfPt,surfBform] = bsplineSurfPts_spapi(Q,3,3,u,v,'paramMethod','centripetal');
    
    % calculate the derivative and cutting direction of the surface
    
    figSurfInt = figure('Name','Freeform surface interpolation results');
    pose = get(gcf,'Position');
    set(gcf,"Position",[pose(1)-pose(3)/2,pose(2),2*pose(3),pose(4)]);
    tilSurfInt = tiledlayout(1,2);
    nexttile;
    plot3(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3), ...
        '.','Color',[0,0.45,0.74]); 
    hold on;
    plot3(reshape(surfCpts(:,:,1),[],1), ...
        reshape(surfCpts(:,:,2),[],1), ...
        reshape(surfCpts(:,:,3),[],1), ...
        'x','Color',[0.85,0.33,0.10]);
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    grid on
    legend('Measured Pts','Interpolation Cpts', ...
        'Orientation','vertical','Location','best');
    nexttile;
    surf(surfPt(:,:,1),surfPt(:,:,2),surfPt(:,:,3), ...
        'FaceColor','interp','EdgeColor','none');
    cb = colorbar('Location','eastoutside');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
end




%% draw the freeform surface
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
parfor ii = 1:ptNum
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),[0;0;1],surfDirect(:,ii));
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
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
    msgbox('Exit for the program','Exit','help','modal');
    uiwait(msgbox);
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
for ii = 1:resNum
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
for ii = 1:ptNum
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
s5_visualize_process;
