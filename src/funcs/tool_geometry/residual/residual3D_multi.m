function [res,peakPt,interPt,uLim] = residual3D_multi(toolPt,toolNorm, ...
    toolCutDir,toolContactU,toolData,toolRadius,uLim,aimRes,uLimIni,curveFunc,varargin)
% to calculate the residual height among the adjacent tool points.
%
% Usage:
%
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDir,
%   toolContactU,sp,toolRadius,uLim1,uLim2,toolPt2,toolNorm2,toolPt3,
%   toolContactU2,toolNorm3,,toolContactU3)
% Inputs:
%   toolPt          (3,1) double    the current point on the tool path
%   toolNorm        (3,1) double    the spindle direction of the current point
%   toolCutDir      (3,1) double    the cutting direction of the current point
%   toolContactU    (1,1) double    
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolRadius      (1,1) double    the tool radius
%   uLim1           (1,2) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,2) double    the interval of the parameter for the
%                                   cloest tool point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolPt2         (3,1) double    the closest point on the tool path
%   toolNorm2       (3,1) double    the spindle direction of "toolPt2"
%   toolPt3         (3,1) double    the 2nd closest point on the tool path
%   toolNorm3       (3,1) double    the spindle direction of "toolPt3"
% Outputs:
%   res             (1,1) double    the residual height of the two point
%   peakPt          (3,1) double    the intersection point
%   uLim1           (1,2) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,2) double    the interval of the parameter for the
%                                   cloest tool point
%
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,
%   toolContactU,sp,toolRadius,uLim1,uLim2,ind1,ind2,ind3)
% Inputs:
%   toolPt          (3,:) double    the list of tool path points
%   toolNorm        (3,:) double    the list of the spindle direction
%   toolCutDir      (3,:) double    the list of the cutting direction
%   toolContactU    (1,:) double    
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolRadius      (1,1) double    the tool radius
%   uLim1           (1,1) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,1) double    the interval of the parameter for the
%                                   cloest tool point
%   ind1            (1,1) double    the index of the current tool path point
%   ind2            (1,1) double    the index of the 1st tool path point
%   ind3            (1,1) double    the index of the 2nd tool path point
% Outputs: same as the above

switch nargin
    case 13
        toolPt1 = toolPt(:,varargin{1});
        toolPt2 = toolPt(:,varargin{2});
        toolPt3 = toolPt(:,varargin{3});
        toolNorm1 = toolNorm(:,varargin{1});
        toolNorm2 = toolNorm(:,varargin{2});
        toolNorm3 = toolNorm(:,varargin{3});
        toolCutDir1 = toolCutDir(:,varargin{1});
%         toolCutDir2 = toolCutDir(:,varargin{2});
        toolContactU1 = toolContactU(varargin{1});
        toolContactUProj = 0.5*( ...
            toolContactU(varargin{2}) + toolContactU(varargin{3}));
    case 16
        toolPt1 = toolPt;
        toolPt2 = varargin{1};
        toolPt3 = varargin{4};
        toolNorm1 = toolNorm;
        toolNorm2 = varargin{2};
        toolNorm3 = varargin{5};
        toolCutDir1 = toolCutDir;
%         toolCutDir2 = varargin{3};
        toolContactU1 = toolContactU;
        toolContactUProj = 0.5*(toolContactU2 + toolContactU3);
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

toolSp = toolData.toolBform;

%% tool point and orientation projection
% method 1: to project the closest point and its tool edge to the current one
% toolPtProj = ptOnPlane(toolPt1,toolPt2,toolCutDir1);

% method 2: to simulate the runin process to get the projection tool edge
% the tool path point projection on the plane of tool path 1
if all(abs(toolPt2(1:2) - toolPt3(1:2)) < 1e-4*toolRadius) % toolPt2 and toolPt3 remains the same position
    toolPtProj = toolPt2;
    if norm(toolNorm2 - toolNorm1) < norm(toolNorm3 - toolNorm1)
        toolNormProj = toolNorm2;
    else
        toolNormProj = toolNorm3;
    end
    toolCutDirProj = toolCutDir1;
else
%     toolPtProj = curveIntersectLineOnPlane(toolPt2,toolPt3,toolPt1);
    toolPtProj = lineIntersectPlane(toolPt2,toolPt3,toolPt1,toolCutDir1);
    t = norm(toolPtProj - toolPt2)/norm(toolPt3 - toolPt2);
    % theta = t*vecAng(toolNorm2,toolNorm3,1);
    % axis = cross(toolNorm2,toolNorm3);
    % RInterp = expm(theta*axis);
    q2 = vecQuat(toolNorm2,toolNorm3);
    qInterp = slerp([1,0,0,0],q2,t);
    toolNormInterp = quat2rotm(qInterp)*toolNorm2;
    % toolCutDirProj = quat2rotm(qInterp)*toolCutDir2;
    toolNormProj = vecOnPlane(toolNormInterp,toolPtProj,toolCutDir1);
    toolNormProj = toolNormProj./norm(toolNormProj);
    % 讲道理两个应该是相等的但是并不等，是否意味着我的刃口并不在对应平面内？
    % toolCutDirProj1 = vecRot(toolNorm2,toolNormProj)*quat2rotm(qInterp)*toolCutDir2; 
    toolCutDirProj = toolCutDir1;
end

%% rigid transform of tool edge from the standard place to the corresponding
R1 = axesRot(toolData.toolEdgeNorm,toolData.cutDirect,toolNorm1,toolCutDir1,'zx');
toolSp1 = toolSp;
toolSp1.coefs = R1*toolSp.coefs + toolPt1;
RProj = axesRot(toolData.toolEdgeNorm,toolData.cutDirect,toolNormProj,toolCutDirProj,'zx');
toolSpProj = toolSp;
toolSpProj.coefs = RProj*toolSp.coefs + toolPtProj;

%% to solve the residual height between toolSp1 and toolSpProj
toolContactPt1 = fnval(toolSp1,toolContactU1);
toolContactPtProj = fnval(toolSpProj,toolContactUProj);

if uLimIni(1) > uLimIni(2)
    uDirection = 'U Minus';
else
    uDirection = 'U Plus';
end

if isempty(uLim)
    uLim2 = uLimIni;
    % the current point has not been calculated
    [res,peakPt,interPt,uLim,~] = residual2D_multi(toolSp1,toolSpProj, ...
        1e-5,toolContactPt1,toolContactPtProj,uLim2,'aimRes',aimRes,'uDirection',uDirection,'curveFunc',curveFunc);
else
    % the current point has been calculated once
    [res,peakPt,interPt,~,uLim] = residual2D_multi(toolSpProj,toolSp1, ...
        1e-5,toolContactPtProj,toolContactPt1,uLim,'aimRes',aimRes,'uDirection',uDirection,'curveFunc',curveFunc);
end

%% to update the valid U range of the two toolpath
% if norm(toolPt1(1:2)) < norm(toolPtProj(1:2))
%     % the projection point stays outside the current point
%     uLim
% else
%     % the projection point stays inside the current point
%     uLim
% end

end