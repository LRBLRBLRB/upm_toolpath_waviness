close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = 'mm';
textFontSize = 14;
textFontType = 'Times New Roman';

%% concentric circles
% [fileName,dirName] = uigetfile({ ...
%     '*.mat','MAT-files(*.mat)'; ...
%     '*,*','all files(*.*)'}, ...
%     'Select one tool edge data file', ...
%     'output_data\tool\tooltheo.mat', ...
%     'MultiSelect','off');
% toolName = fullfile(dirName,fileName);
toolName = 'output_data\tool\toolTheo.mat';
toolData = load(toolName);

default = false;
if default
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select the surface edge data file', ...
        'input_data\surface\ellipsoidAnay.mat', ...
        'MultiSelect','off');
    surfName = fullfile(dirName,fileName);
    load(surfName);
% % % % % % % % % % % % % % % % % % % % % % % % %     sparTheta
% % % % % % % % % % % % % % % % % % % % % % % % %     surfMesh
% % % % % % % % % % % % % % % % % % % % % % % % %     surfNorm
else % ellipsoid
    R = 10*2;
    A = 3.5*2;
    B = 4*2;
    C = 5*2;
    % sampling density
    r = [0,3*R/4]; % concentric radius range
    sparTheta = 101;
    surfCenter = [0,0,sqrt(C^2*(R.^2-r(2).^2))]; % concentric circle center
    conR = (toolData.radius/2):(toolData.radius):r(2); % concentric radius vector
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
%     save('Simulation/ellipsoidAnay.mat', ...
%         "conCenter","conR","surfMesh","surfNorm","surfCenter");
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
legend('Original Points','Location','best');
quiver3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    surfNorm(1:theSpar:end,1:rSpar:end,1), ...
    surfNorm(1:theSpar:end,1:rSpar:end,2), ...
    surfNorm(1:theSpar:end,1:rSpar:end,3), ...
    'AutoScale','on','Color',[0.85,0.33,0.10],'DisplayName','Normal Vectors');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x(',unit,'m)']);
ylabel(['y(',unit,'m)']);
zlabel(['z(',unit,'m)']);
axis equal; grid on;
% title({'Radially & circunferentially sparse','by 50 and 2 times, respectively'}, ...
%     'FontSize',textFontSize,'FontName',textFontType);
nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on; axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x(',unit,'m)']);
ylabel(['y(',unit,'m)']);
zlabel(['z(',unit,'m)']);
% cb2 = colorbar;
msgfig = msgbox('Surface was generated successfully!', ...
    'Surface Generation','warn','non-modal');
uiwait(msgfig);


%% Calculation of Tool Path & Spindle Drection
surfPt = transpose(reshape(surfMesh,[],3));
surfNorm = transpose(reshape(surfNorm,[],3));
surfDirect = cutDirection(surfPt,surfCenter);
ptNum = size(surfPt,2);
toolQuat = zeros(ptNum,4);
toolVec = zeros(3,ptNum);
toolPathPt = zeros(3,ptNum);
toolCutDirect = zeros(3,ptNum);
toolNormDirect = zeros(3,ptNum);
isCollision = false(1,ptNum);
% parpool(6);
tic
parfor ii = 1:ptNum
    [toolQuat(ii,:),toolVec(:,ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),surfDirect(:,ii),[0;0;1]);
    if isCollision(ii) == true
        pause('on');
    else
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
end
toc;

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
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x(',unit,'m)']);
ylabel(['y(',unit,'m)']);
zlabel(['z(',unit,'m)']);
legend('tool center point','tool cutting direction', ...
    'tool spindle direction','','Location','best');

%% Calculation of Residual Height & Cutting Surface
toolSp = toolData.toolBform;
res = zeros(1,ptNum);
tmpInd = linspace(1,ptNum);
uLim = [zeros(1,ptNum);ones(1,ptNum)];
tmpULim = [zeros(1,ptNum);ones(1,ptNum)];
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));
parfor ii = 1:ptNum
    nLoop = ceil(ii/sparTheta);
    closest = find(abs(angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1)) - angle(ii)),3);
    ind2 = sparTheta*nLoop + closest(2);
    ind3 = sparTheta*nLoop + closest(2);
    [res(ii),tmpInd(ii),uLim(:,ii),tmpULim(:,ii)] = residual3DPar( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolSp,ii,ind2,ind3);
end
clear angle;




%% Visualization & Simulation





%% 
% rmpath(genpath('funcs'));