function [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDir,toolContactU,toolSp,uLim1,uLim2,varargin)
% to calculate the residual height among the adjacent tool points.
%
% Usage:
%
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,
%   toolContactU,sp,uLim1,uLim2,toolPt2,toolNorm2,toolCutDir2,toolPt3,
%   toolContactU2,toolNorm3,toolCut3,toolContactU3)
% Inputs:
%   toolPt          (3,1) double    the current point on the tool path
%   toolNorm        (3,1) double    the spindle direction of the current point
%   toolCutDir      (3,1) double    the cutting direction of the current point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolPt2         (3,1) double    the closest point on the tool path
%   toolNorm2       (3,1) double    the spindle direction of "toolPt2"
%   toolCutDir2     (3,1) double    the cutting direction of "toolPt2"
%   toolPt3         (3,1) double    the 2nd closest point on the tool path
%   toolNorm3       (3,1) double    the spindle direction of "toolPt3"
%   toolCutDir3     (3,1) double    the cutting direction of "toolPt3"
% Outputs:
%   res             (1,1) double    the residual height of the two point
%   peakPt          (3,1) double    the intersection point
%   uLim1           (1,2) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,2) double    the interval of the parameter for the
%                                   cloest tool point
%
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,
%   toolContactU,sp,uLim1,uLim2,ind1,ind2,ind3)
% Inputs:
%   toolPt          (3,:) double    the current point on the tool path
%   toolNorm        (3,:) double    the spindle direction of the current point
%   toolCutDir      (3,:) double    the cutting direction of the current point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   ind1
%   ind2
%   ind3
% Outputs: same as the above

switch nargin
    case 10
        toolPt1 = toolPt(:,varargin{1});
        toolPt2 = toolPt(:,varargin{2});
        toolPt3 = toolPt(:,varargin{3});
        toolNorm1 = toolNorm(:,varargin{1});
        toolNorm2 = toolNorm(:,varargin{2});
        toolNorm3 = toolNorm(:,varargin{3});
        toolCutDir1 = toolCutDir(:,varargin{1});
        toolCutDir2 = toolCutDir(:,varargin{2});
        toolCutDir3 = toolCutDir(:,varargin{3});
        toolContactU1 = toolContactU(varargin{1});
        toolContactUProj = 0.5*( ...
            toolContactU(varargin{2}) + toolContactU(varargin{3}));
    case 15
        toolPt1 = toolPt;
        toolPt2 = varargin{1};
        toolPt3 = varargin{4};
        toolNorm1 = toolNorm;
        toolNorm2 = varargin{2};
        toolNorm3 = varargin{5};
        toolCutDir1 = toolCutDir;
        toolCutDir2 = varargin{3};
        toolCutDir3 = varargin{6};
        toolContactU1 = toolContactU;
        toolContactUProj = 0.5*(toolContactU2 + toolContactU3);
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

% to project the closest point and its tool edge to the current one
toolPtProj = ptOnPlane(toolPt1,toolPt2,toolCutDir1);
% Van der Waals: 主流的做法是对位置插补 然后刀轴是跟随运动；位置已经运动的长度/总长度=已经转过的角度/总角度
% 是否可以直接用向量插值、或者四元数来代替上述公式里的角度？
% 这里的问题：对于螺旋线，这种投影方式是否失效？
t = norm(toolPtProj - toolPt2)/norm(toolPt3 - toolPt2);
toolNormProj = (1 - t)*toolNorm2 + t*toolNorm3;
toolNormProj = vecOnPlane(toolNormProj,toolPtProj,toolCutDir1);
toolNormProj = toolNormProj./norm(toolNormProj);
% toolCutDirProj = (1 - t)*toolCutDir2 + t*toolCutDir3;
% toolCutDirProj = toolCutDirProj + toolCutDir1*(dot(toolCutDirProj,toolCutDir1) - 1);
toolCutDirProj1 = toolCutDir1; % 讲道理应该是同向的
R = vecRot(toolNorm2,toolNormProj);
toolCutDirProj = R*toolCutDir2; % 讲道理两个应该是相等的但是并不等，是否意味着我的刃口并不在对应平面内？

% rigid transform of tool edge from the standard place to the corresponding
R1 = axesRot(toolNorm1,toolCutDir1,[0;0;1],[1;0;0],'zx');
toolSp1 = toolSp;
toolSp1.coefs = R1*toolSp.coefs + toolPt1;
RProj = axesRot(toolNormProj,toolCutDirProj,[0;0;1],[1;0;0],'zx');
toolSpProj = toolSp;
toolSpProj.coefs = RProj*toolSp.coefs + toolPtProj;

% to solve the residual height between toolSp1 and toolSpProj
eps = 1e-3;
u = 0:eps:1;
Pt1 = fnval(toolSp1,u);
PtProj = fnval(toolSpProj,u);
[peakPt,ind1,indProj,~] = pcCrossDN(Pt1,PtProj);
toolContactPt1 = fnval(toolSp1,toolContactU1);
toolContactPtProj = fnval(toolSpProj,toolContactUProj);
res = norm(cross(toolContactPt1 - peakPt,toolContactPtProj - peakPt)) ...
    /norm(toolContactPtProj - toolContactPt1);
if res >= norm(Pt1(:,1)-Pt1(:,2)) % no intersection
    res = nan;
end

% to update the valid U range of the two toolpath
u1 = u(ind1);
uProj = u(indProj);
if u1 > uProj
    uLim1(2) = min(u1,uLim1(2));
    uLim2(1) = max(uProj,uLim2(1));
else
    uLim1(1) = max(u1,uLim1(1));
    uLim2(2) = min(uProj,uLim2(2));
end

end