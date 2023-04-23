% function pt = bSplineCurve(Q,k,u)
% use MATLAB functions to interpolate the curve, and compare it with the
% function that I write
% Inputs:
%   Q (n,3) matrix of data points
%   k (1,1) order of B spline（次数）

% arguments
%     Q {mustBeFinite}
%     k (1,1) {mustBeInteger} = 3
% end

clear;
clc;
addpath(genpath('..\'));
Q = [2 1.4;1 .5; 2 -.4; 5 1.4; 6 .5; 5 -.4];
k = 2;

%% my spline interpolation
% [cpts,U] = bSplineCpts(Q,k,'chord'); 
cpts = [2 1.4;1 .5; 2 -.4; 5 1.4; 6 .5; 5 -.4];
U = linspace(0,1,6);
UU = augknt(U,k+1);
u = 0:0.02:1;
nPts = length(u);
toolDim = size(cpts,2);
pt0 = zeros(nPts,toolDim);
% for i = 1:nPts
%     pt0(i,:) = bSplinePt(cpts,3,u(i),U);
% end

%% spline interpolation by Curve fitting toolbox
k = k + 1; % get the order of B spline (the same as the curvefitter toolbox)

% b-spline genrerating
sp = spmak(U,cpts');
pt = transpose(fnval(sp,u));
figure;
fnplt(sp,U([3 7]))

rmpath(genpath('..\'));
% end