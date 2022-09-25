function [r,ang,scatterDst,varargout] = toolFit3D(scatterOri,fitMethod)
% usage: [c,r,ang,scatterDst,RMSE] = toolFit3D(scatterOri,options)
%
% solve the edge sharpness of a arc turning tool: 
% calculate the radius of the tool tip arc, 
% as well as transforming the tool coordinate
%
% Inputs: 
%   scatterMat: original tool edge point from measuring (3,n)
%   options: 
%       fitMethod: the method of curve fitting
% Outputs: 
%   r: radius of the arc (1,1)
%   ang: the open angle of the tool (1,1)
%   scatterDst: tool edge points with pose adjustment (3,n)
%   c: center of the arc (3,1)
%   RMSE: the square-mean-root error of curve fitting
%
% method:
%   least square method by normal equation solving

arguments
    scatterOri (3,:) {mustBeFinite}
    fitMethod {mustBeMember(fitMethod, ...
        ['Gradient-decent','Normal-equation','Levenberg-Marquardt'])} ...
        = 'Levenberg-Marquardt'
end

n = size(scatterOri,2);

%% Method one: least-square fitting
% project the scatters to the least square fitting plane
% B = -1*ones(n,1);
% flat = lsqminnorm(scatterOri,B);
% % 投影再坐标转换和直接转换结果一样？
% k = (flat(1)*scatterOri(:,1) + flat(2)*scatterOri(:,2) + flat(3)*scatterOri(:,3)) ...
%     ./(flat(1)*flat(1) + flat(2)*flat(2) + flat(3)*flat(3));
% x = -flat(1)*k + scatterOri(:,1);
% y = -flat(2)*k + scatterOri(:,2);
% z = -flat(3)*k + scatterOri(:,3);
% 
% % coordinate transform to the plane
% 
% 
% 
% 
% % 2D circle fitting
% [cc,rr,ang,scatterPlaneDst,RMSE] = circleFit2D([xx,yy],options);
%
% reverse coordinate transform to 3D

%% Method two: SVD & least square fitting
% SVD to fit the plane
scatterCentered = scatterOri' - meshgrid(mean(scatterOri,2),1:n);
% svd: singular values are nonnegative and returned in decreasing order.
[~,~,V] = svd(scatterCentered,'vector'); 
planeNorm = V(:,3); % normVec remains the norm vector of th e fitting plane

% rotNorm = transpose(cross(planeNorm,[0;0;1]));
% rotNorm = rotNorm/norm(rotNorm);
% planeAng = vecAng(planeNorm,[0;0;1],1);
% scatterPlane = Rodrigues(scatterOri,rotNorm,planeAng);
scatterR = vecRot(planeNorm,[0;0;1]);
scatterPlane = scatterR*scatterOri;

% 2D circle fitting
% if any(scatterPlane(3,:))
%     error('Error: uncorrect projection from plane to xOz.');
% end
[c,r,ang,RMSE,startV,endV] = circleFit2D(scatterPlane(1:2,:),fitMethod);

% rigid transform: to get the standardized tool profile
mid = 0.5*(startV + endV);
rotAng = pi/2 - atan2(mid(2),mid(1));
rotMat = rotz(rotAng);
rotMat = rotMat(1:2,1:2);
scatterDst = rotMat*(scatterPlane(1:2,:) - ndgrid(c,1:n));

% optional output
varargout{1} = RMSE;
end