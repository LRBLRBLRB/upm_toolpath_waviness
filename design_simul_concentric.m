close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
textFontSize = 14;
textFontType = 'Times New Roman';

%% concentric circles
% [fileName,dirName] = uigetfile({ ...
%     '*.mat','MAT-files(*.mat)'; ...
%     '*,*','all files(*.*)'}, ...
%     'Select one tool edge data file', ...
%     'Tool\tooltheo.mat', ...
%     'MultiSelect','off');
% toolName = fullfile(dirName,fileName);
toolName = 'output_data\tool\tooltheo.mat';
toolData = load(toolName);

default = false;
if default
    load('output_data/simulation/ellipsoidAnay.mat');
else % ellipsoid
    R = 1;
    A = .3;
    B = .4;
    C = .5;
    
    % sampling density
    r = [0,3*R/4]; % concentric radius range
    densTheta = 101; % 问题：r和theta的离散度不同会不会有影响？
    conCenter = [0,0,sqrt(C^2*(R.^2-r(2).^2))]; % concentric circle center
    conR = (toolData.center/2):(toolData.center):r(2); % concentric radius vector
    densR = length(conR);
    conTheta = linspace(0,2*pi,densTheta);
    
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = sqrt(C^2*(R.^2-rMesh.^2));
    %  y = reshape(surfMesh(:,:,2),[],1);
    % z = reshape(surfMesh(:,:,3),[],1);
    % surfXYZ = [x,y,z];
    
    % syms X Y Z(X,Y);
    % Z(X,Y) = sqrt(C^2*(R^2 - X.^2/A^2 - Y.^2/B^2));
    % ZDX = diff(Z,X);
    % ZDY = diff(Z,Y);
    % surfNorm(:,1) = eval(subs(ZDX,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,2) = eval(subs(ZDY,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,3) = -ones(densTheta*densR,1);
    [surfNorm(:,:,1),surfNorm(:,:,2),surfNorm(:,:,3)] = surfnorm( ...
        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));

    % save('Simulation/ellipsoidAnay.mat',"conCenter","conR","surfMesh","surfNorm");
end

%% plot the freeform surface
figure('Name','original xyz scatters of the surface (sparsely)');
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
tiledlayout(1,2);
nexttile;
rSpar = 50;
theSpar = 2;
plot3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    '.','Color',[0,0.45,0.74]);
hold on; axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');
quiver3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    surfNorm(1:theSpar:end,1:rSpar:end,1), ...
    surfNorm(1:theSpar:end,1:rSpar:end,2), ...
    surfNorm(1:theSpar:end,1:rSpar:end,3), ...
    'AutoScale','on','Color',[0.85,0.33,0.10]);
legend('Original Points','Normal Vectors','Location','best');
title({'Radially & circunferentially sparse','by 50 and 2 times, respectively'});
nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
cb = colorbar;
msgfig = msgbox('Surface was generated successfully!', ...
    'Surface Generation','warn','modal');
uiwait(msgfig);


%% Tool Edge Simulation
surfPt = reshape(surfMesh,[],3);
surfNorm = reshape(surfNorm,[],3);
ptNum = size(surfPt,1);
toolCenPt = zeros(ptNum,3);
for ii = 1:ptNum
    toolCenPt(ii,:) = toolPos(toolData,surfPt(ii,:),surfNorm(ii,:));
end

figure('Name','tool center position & tool normal vector');
plot3(toolCenPt(1:100:end,1), ...
    toolCenPt(1:100:end,2), ...
    toolCenPt(1:100:end,3), ...
    '.','Color',[0.85,0.33,0.10]);
hold on; axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');

plot3(surfPt(1:100:end,1), ...
    surfPt(1:100:end,2), ...
    surfPt(1:100:end,3), ...
    '.','Color',[0,0.45,0.74]);




%% 
rmpath(genpath('funcs'));