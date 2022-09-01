function [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDir,toolSp,varargin)
% to calculate the residual height among the adjacent tool points.
%
% Usage:
%
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,sp,
%   toolPt2,toolNorm2,toolCutDir2,toolPt3,toolNorm3,toolCut3)
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
% [res,peakPt,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,sp,
%   ind1,ind2,ind3)
% Inputs:
%   toolPt          (3,:) double    the current point on the tool path
%   toolNorm        (3,:) double    the spindle direction of the current point
%   toolCutDirect   (3,:) double    the cutting direction of the current point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   ind1
%   ind2
%   ind3
% Outputs: same as the above

switch nargin
    case 7
        toolPt1 = toolPt(:,varargin{1});
        toolPt2 = toolPt(:,varargin{2});
        toolPt3 = toolPt(:,varargin{3});
        toolNorm1 = toolNorm(:,varargin{1});
        toolNorm2 = toolNorm(:,varargin{2});
        toolNorm3 = toolNorm(:,varargin{3});
        toolCutDir1 = toolCutDir(:,varargin{1});
        toolCutDir2 = toolCutDir(:,varargin{2});
        toolCutDir3 = toolCutDir(:,varargin{3});
    case 10
        toolPt1 = toolPt;
        toolPt2 = varargin{1};
        toolPt3 = varargin{4};
        toolNorm1 = toolNorm;
        toolNorm2 = varargin{2};
        toolNorm3 = varargin{5};
        toolCutDir1 = toolCutDir;
        toolCutDir2 = varargin{3};
        toolCutDir3 = varargin{6};
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

% to project the closest point and its tool edge to the current one
toolPtProj = projectionOnPlane(toolPt1,toolPt2,toolCutDir1);
% Van der Waals: 主流的做法是对位置插补 然后刀轴是跟随运动；位置已经运动的长度/总长度=已经转过的角度/总角度
% 是否可以直接用向量插值、或者四元数来代替上述公式里的角度？
% 这里的问题：对于螺旋线，这种投影方式是否失效？
t = norm(toolPtProj - toolPt2)/norm(toolPt3 - toolPt2);
toolNormProj = (1 - t)*toolNorm2 + t*toolNorm3;
toolCutDirProj = (1 - t)*toolCutDir2 + t*toolCutDir3;

% rigid transform of tool edge from the standard place to the corresponding
R1 = axesRot(toolNorm1,toolCutDir1,[0;0;1],[1;0;0]);
toolSp1 = toolSp;
toolSp1.coefs = R1*toolSp.coefs + toolPt1;
RProj = axesRot(toolNormProj,toolCutDirProj,[0;0;1],[1;0;0]);
toolSpProj = toolSp;
toolSpProj.coefs = RProj*toolSp.coefs + toolPtProj;

% 


iterSolveSp(toolSp1);

















end