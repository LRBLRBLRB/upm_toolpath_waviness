% to recur the work by He Shuai's STPGM-NURBS
close all;
clear;
clc;
addpath(genpath('../funcs'));
t0 = tic;

% global variables
% global textFontSize textFontType;
workspaceDir = '../workspace/20220925-contrast/nagayama_concentric';
unit = '\mum';
textFontSize = 14;
textFontType = 'Times New Roman';
% profile on
parObj = gcp('nocreate');

% machining paramters, unit: mu m
toolRadius = 200;
aimRes = 10;
arcLength = 30;
maxAngPtDist = 6*pi/180;

% machining surface
R = 10/2*1000;
A = 3.5/2;
B = 4/2;
C = 5/2;
syms x y;
surfSym = C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2);
surfFunc = matlabFunction(surfSym);
surfFx = matlabFunction(diff(surfFunc,x));
surfFy = matlabFunction(diff(surfFunc,y));
rLin = linspace(0,R/4,101); % concentric radius vector
thetaLin = linspace(0,2*pi,101);
[rMesh,thetaMesh] = meshgrid(rLin,thetaLin);
xMesh = rMesh.*cos(thetaMesh);
yMesh = rMesh.*sin(thetaMesh);
zMesh = surfFunc(xMesh,yMesh);
surfOri = cat(3,rMesh,thetaMesh,zMesh);

%% NURBS interpolation
u = 0:0.001:1;
v = 0:0.001:1;
% actually he uses the NURBS to interpolate the surface
[surfPtsPolar,surfBform] = bsplineSurfPts_spapi(surfOri,3,3,u,v, ...
    'paramMethod','concentric','cptsType','Polar');
surfPtsPolar = permute(surfPtsPolar,[2,3,1]);
surfPts(:,:,1) = surfPtsPolar(:,:,1).*cos(surfPtsPolar(:,:,2));
surfPts(:,:,2) = surfPtsPolar(:,:,1).*sin(surfPtsPolar(:,:,2));
surfPts(:,:,3) = surfPtsPolar(:,:,3);
figure;
scatter3(xMesh,yMesh,zMesh,6,'MarkerEdgeColor',[0.635,0.078,0.1840],'MarkerEdgeAlpha',0.5);
hold on;
surf(surfPts(:,:,1),surfPts(:,:,2),surfPts(:,:,3), ...
    'FaceColor','interp','EdgeColor','none');

%% to generate cutting points









%% Tool radius and clearance angle compensation






%%
tTol = toc(t0);
% rmpath(genpath('../funcs'));