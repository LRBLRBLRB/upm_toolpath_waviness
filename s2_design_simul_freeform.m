close all;
clear;
clc;
addpath(genpath('funcs'));

%% load mesh data points of a designed freeform surface
[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)' ...
    }, ...
    'Select one surface point cloud file', ...
    'Surface\sphereHoleXYZ.mat', ...
    'MultiSelect','off');
pathName = fullfile(dirName,fileName);
surfData = load(pathName);
dim = size(surfData.xyz,1);

f1 = figure('Name','original xyz scatters of the surface');
plot3(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3),'.');
hold on; axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('z');

%% B-spline interpolation of the freeform surface
[surfCpts,U,V] = bSplineSurfCpts(surfData.mesh,3,3);
u = 0:0.002:1;
v = 0:0.002:1;
% surfPts = zeros(length(u),length(v),dim);
surfPts = bSplineSurfPts(surfCpts,3,u,U,3,v,V);


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
surf(surfPts(:,:,1),surfPts(:,:,2),surfPts(:,:,3), ...
    'FaceColor','interp','EdgeColor','none');
cb = colorbar('Location','eastoutside');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal

%% Simulaton with Tool Model

% tool
% 
% spiral
% 
% proj
% 
% tangent


%% 
rmpath(genpath('funcs'));