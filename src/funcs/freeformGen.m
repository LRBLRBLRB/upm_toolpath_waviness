% to generate a freeform surface by multiple methods
clear;
clc;
addpath(genpath('funcs'));

%% surface mesh generation
simul = 'ellipsoid';

switch simul
    case 'sphereHoleCartesian' %% Spheric surface with a pit in XYOZ
        x = -5:0.01:5;
        y = -5:0.01:5;
        R = 20;
        [xMesh,yMesh] = meshgrid(x,y);
        zMesh = sqrt(R.^2-xMesh.^2-yMesh.^2);
    case 'sphereHoleRadial' %% Spheric surface with a pit in R\thetaOZ
        R = 20;
        r = linspace(0,3*R/4,101);
        theta = linspace(0,2*pi,101);
        [rMesh,thetaMesh] = meshgrid(r,theta);
        xMesh = rMesh.*cos(thetaMesh);
        yMesh = rMesh.*sin(thetaMesh);
        zMesh = sqrt(R.^2-rMesh.^2);
    case 'ellipsoid' %% Ellipsoid Generation
        dens = 101;
        R = 100;
        A = 30;
        B = 40;
        C = 50;
        r = linspace(0,3*R/4,dens);
        theta = linspace(0,2*pi,dens);
        [rMesh,thetaMesh] = meshgrid(r,theta);
        xMesh = A*rMesh.*cos(thetaMesh);
        yMesh = B*rMesh.*sin(thetaMesh);
        zMesh = sqrt(C^2*(R.^2-rMesh.^2));
        center = [0,0,30];
end

figure;
mesh(xMesh,yMesh,zMesh);
colorbar
xlabel('x');
ylabel('y');
zlabel('z');
axis equal

mesh(:,:,1) = xMesh;
mesh(:,:,2) = yMesh;
mesh(:,:,3) = zMesh;
x = reshape(xMesh,[],1);
y = reshape(yMesh,[],1);
z = reshape(zMesh,[],1);
xyz = [x,y,z];
fileName = fullfile('Surface',[simul,'XYZ',num2str(dens),'.mat']);
save(fileName,'center','xyz','mesh');


%% Surface from .txt file
% SurfPC = readmatrix("Surface\sphereHole3000.txt",'NumHeaderLines',1);
% xyz = SurfPC(:,1:3);
% del = xyz(:,3)-2.5 == 0;
% xyz(del,:) = [];
% 
% figure;
% plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','Color',[0,0.45,0.74]);
% % [0,0.45,0.74] [0.49,0.18,0.56]
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis equal
% 
% save Surface\surfaceXYZ.mat xyz;

rmpath(genpath('funcs'));