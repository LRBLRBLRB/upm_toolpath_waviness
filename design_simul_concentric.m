close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = 'mm';
textFontSize = 14;
textFontType = 'Times New Roman';
msgOpts.Default = 'Cancel';
msgOpts.Interpreter = 'tex';
% msgOpts.modal = 'non-modal';
profile on
parObj = gcp;

%% concentric circles
% [fileName,dirName] = uigetfile({ ...
%     '*.mat','MAT-files(*.mat)'; ...
%     '*,*','all files(*.*)'}, ...
%     'Select one tool edge data file', ...
%     'output_data\tool\tooltheo.mat', ...
%     'MultiSelect','off');
% toolName = fullfile(dirName,fileName);
toolName = 'output_data\tool\toolTheo_3D.mat';
toolData = load(toolName);

default = false;
if default
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select the surface edge data file', ...
        'input_data\surface\ellipsoidAray.mat', ...
        'MultiSelect','off');
    surfName = fullfile(dirName,fileName);
    load(surfName);
% % % % % % % % % % % % % % % % % % % % % % % % %     sparTheta
% % % % % % % % % % % % % % % % % % % % % % % % %     surfMesh
% % % % % % % % % % % % % % % % % % % % % % % % %     surfNorm
else % ellipsoid
    R = 10/2;
    A = 3.5/2;
    B = 4/2;
    C = 5/2;
    % sampling density
    r = [0,R/4]; % concentric radius range
    sparTheta = 101;
    surfCenter = [0,0,sqrt(C^2*(R.^2-r(2).^2))]; % concentric circle center
    conR = (toolData.radius/2):(toolData.radius/2):r(2); % concentric radius vector
    densR = length(conR);
    conTheta = linspace(0,2*pi,sparTheta);
    
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = sqrt(C^2*(R.^2-rMesh.^2));
    %  y = reshape(surfMesh(:,:,2),[],1);
    % z = reshape(surfMesh(:,:,3),[],1);
    % surfXYZ = [x,y,z];
    
    % calculate the normal vector of the analytic surface
    % syms X Y Z(X,Y);
    % Z(X,Y) = sqrt(C^2*(R^2 - X.^2/A^2 - Y.^2/B^2));
    % ZDX = diff(Z,X);
    % ZDY = diff(Z,Y);
    % surfNorm(:,1) = eval(subs(ZDX,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,2) = eval(subs(ZDY,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,3) = -ones(densTheta*densR,1);
    [surfNorm(:,:,1),surfNorm(:,:,2),surfNorm(:,:,3)] = surfnorm( ...
        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));
%     save('input_data/surface/ellipsoidAray.mat', ...
%         "surfMesh","surfNorm","surfCenter");
end

%% plot the freeform surface
figure('Name','original xyz scatters of the surface (sparsely)');
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
tiledlayout(1,2);
nexttile;
rSpar = 1;
theSpar = 1;
plot3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    '.','Color',[0,0.45,0.74]);
hold on;
legend('Original Points','Location','northeast');
quiver3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    surfNorm(1:theSpar:end,1:rSpar:end,1), ...
    surfNorm(1:theSpar:end,1:rSpar:end,2), ...
    surfNorm(1:theSpar:end,1:rSpar:end,3), ...
    'AutoScale','on','Color',[0.85,0.33,0.10],'DisplayName','Normal Vectors');
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
zlabel(['z(',unit,')']);
axis equal; grid on;
% title({'Radially & circunferentially sparse','by 50 and 2 times, respectively'}, ...
%     'FontSize',textFontSize,'FontName',textFontType);
nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on; axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
zlabel(['z(',unit,')']);
% cb2 = colorbar;
msgfig = msgbox('Surface was generated successfully!', ...
    'Surface Generation','warn','non-modal');
uiwait(msgfig);
surfPt = transpose(reshape(surfMesh,[],3));
surfNorm = transpose(reshape(surfNorm,[],3));
surfDirect = cutDirection(surfPt,surfCenter);

%% Calculation of Tool Path & Spindle Drection
ptNum = size(surfPt,2);
toolQuat = zeros(ptNum,4);
toolVec = zeros(3,ptNum);
toolPathPt = zeros(3,ptNum);
toolCutDirect = zeros(3,ptNum);
toolNormDirect = zeros(3,ptNum);
toolContactU = zeros(1,ptNum);
isCollision = false(1,ptNum);
tic
parfor ii = 1:ptNum
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),surfDirect(:,ii),[0;0;1]);
    if isCollision(ii) == true
        pause('on');
    else
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
end
tToolpath = toc;
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
% plot3(surfPt(1,1:100:end), ...
%     surfPt(2,1:100:end), ...
%     surfPt(3,1:100:end), ...
%     '.','Color',[0,0.45,0.74]);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.6,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = 'Height (mm)';
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
zlabel(['z(',unit,')']);
legend('tool center point','tool cutting direction', ...
    'tool spindle direction','','Location','northeast');
msgfig = msgbox('Tool Path was calculated successfully!', ...
    'Tool Path Simulation','warn','non-modal');
uiwait(msgfig);

%% Calculation of Residual Height & Cutting Surface
toolSp = toolData.toolBform;
res = zeros(1,ptNum);
peakPt = zeros(3,ptNum);
uLim = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of each toolpath
tmpULim = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of the projective toolpath
ind2 = linspace(1,ptNum,ptNum);
ind3 = linspace(1,ptNum,ptNum);
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));
resNum = ptNum - sparTheta;
tic
parfor ii = 1:resNum
    nLoop = ceil(ii/sparTheta);
    closest = find(abs(angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1)) - angle(ii)),3);
    ind2(ii) = sparTheta*nLoop + closest(1);
    ind3(ii) = sparTheta*nLoop + closest(2);
    [res(ii),peakPt(:,ii),uLim(:,ii),tmpULim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp, ...
        uLim(:,ii),tmpULim(:,ii),ii,ind2(ii),ind3(ii));
end
tRes = toc;
fprintf('The time spent in the residual height calculation process is %fs.\n',tRes);

for ii = 1:ptNum
    uLim(1,ii) = max([uLim(1,ii),tmpULim(1,ind2(ii))]);
    uLim(2,ii) = min([uLim(2,ii),tmpULim(2,ind2(ii))]);
end
clear angle ind2 ind3 tmpULim;

figure;
plotSpar = 1;
tic
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
toolCoefs = toolData.toolBform.coefs;
for ii = 1:ptNum
    toolSp = toolData.toolBform;
    toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp,uLim(1,ii):0.01:uLim(2,ii));
    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
end
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.6,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = 'Height (mm)';
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x(',unit,')']);
ylabel(['y(',unit,')']);
zlabel(['z(',unit,')']);
legend('tool center point','tool cutting direction', ...
    'tool spindle direction','','','Location','northeast');
tPlot = toc;
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);
msgfig = questdlg({'Residual Height was calculated successfully!', ...
    'Ready for machining simulation?'}, ...
    'Tool Path Simulation','OK','Cancel',msgOpts);
if strcmp(msgfig,'Cancel')
    % delete(parObj);
    profile off
    % profsave(profile("info"),"profile_data");
    return; 
end

%% Visualization & Simulation
figure;
stepLength = 0.01;
nLoop = ceil(ptNum/sparTheta);
uLimRound = round(uLim,2);
toolPathMesh = [];
tic
for ii = 1:nLoop % each loop
    Q = cell(sparTheta,1);
    for jj = 1:sparTheta
        toolSp = toolData.toolBform;
        toolSp.coefs = quat2rotm(toolQuat((ii-1)*nLoop+jj,:))*toolCoefs + toolVec(:,(ii-1)*nLoop+jj);
        Q{jj} = fnval(toolSp,uLimRound(1,(ii-1)*nLoop+jj):stepLength:uLimRound(2,(ii-1)*nLoop+jj));
%         Q{jj} = tmp;
    end
    for u = 0:stepLength:1
        for jj = 1:sparTheta
            if u >= uLimRound(1,(ii-1)*nLoop+jj) && u <= uLimRound(2,(ii-1)*nLoop+jj)
                tmp = Q{jj}(:,round((u - uLimRound(1,(ii-1)*nLoop+jj))/stepLength + 1));
                toolPathMesh = [toolPathMesh,tmp];
            end
        end
    end
end
tSimul = toc;
fprintf('The time spent in the simulation calculation process is %fs.\n',tSimul);
plot3(toolPathMesh(1,:),toolPathMesh(2,:),toolPathMesh(3,:),'.','Color',[0,0.4450,0.7410]);
hold on;
grid on;

%% 
% delete(parObj);
profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));