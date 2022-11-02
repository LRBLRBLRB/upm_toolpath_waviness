function [circ3D,RMSE,varargout] = arcFit3D(scatterOri,options)
% usage:
%   [c,r,ang,RMSE] = arcFit3D(scatterOri,fitMethod)
%   [c3D,r,ang,RMSE,scatterPlane,c2D,startV,endV] = arcFit3D(scatterOri,fitMethod)
%
% fit the circle from a cluster of 3D points
%
% Inputs: 
%   scatterMat: original 3D points (3,n)
%   fitMethod: the method of curve fitting
% Outputs: 
%   center3D: center of the arc (3,1)
%   radius: radius of the arc (1,1)
%   ang: the open angle of the tool (1,1)
%   RMSE: the square-mean-root error of curve fitting
%   startV: start edge vector of the arc
%   endV: end edge vector of the arc
%
% method:
%   least square method by normal equation solving


arguments
    scatterOri (3,:) {mustBeFinite}
    options.arcFitMethod {mustBeMember(options.arcFitMethod, ...
        ['gradient-decent','normal-equation','levenberg-marquardt'])} ...
        = 'levenberg-marquardt'
    options.displayType {mustBeMember(options.displayType, ...
        ['off','none','iter','iter-detailed','final','final-detailed'])} = 'final'
end

n = size(scatterOri,2);

%% Method One: Least Square Fitting
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
% [cc,rr,ang,scatterPlaneDst,RMSE] = arcFit2D([xx,yy],options);
%
% reverse coordinate transform to 3D

%% Method Two: SVD & Least Square Fitting
% SVD to fit the plane
scatterMean = mean(scatterOri,2);
scatterCentered = scatterOri' - meshgrid(scatterMean,1:n);
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
[circ2D,RMSE] = arcFit2D(scatterPlane(1:2,:), ...
    'arcFitMethod',options.arcFitMethod,'displayType',options.displayType);

% Inverse rigid transform to find the fitting center of the original data
scatterPlaneZ = mean(scatterPlane(3,:));
center3D = [circ2D{1};scatterPlaneZ];
center3D = scatterR'*center3D;
arcVec = [circ2D{4};0];
arcVec = scatterR'*arcVec;

circ3D = cell(5,1);
circ3D{1} = center3D;
circ3D{2} = circ2D{2};
circ3D{3} = circ2D{3};
circ3D{4} = arcVec/norm(arcVec);
circ3D{5} = planeNorm/norm(planeNorm);

%% Method Three: Projection & Least Square Fitting
% Method From 化春键, 熊雪梅, 陈莹. 基于拉格朗日乘子法的空间圆弧拟合优化方法[J]. 
% 工程设计学报, 2018, 25(06): 661-667.

% sysMat = [scatterOri';zeros(n,1)];
% planeParam = 
% 
% circ3D{1} = center3D;
% circ3D{2} = radius;
% circ3D{3} = ang;
% circ3D{4} = arcVec/norm(arcVec);
% circ3D{5} = planeNorm/norm(planeNorm);

%% optional output
varargout{1} = scatterPlane;
varargout{2} = circ2D;

end