function [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp( ...
    interpU,toolPt,toolNorm,toolCutDir,toolU,varargin)
% to linear interpolate the exact tool path between two known tool points
%
% Usage:
%
% [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp(
%   toolPt,toolNorm,toolCutDir,toolU,ind1,ind2,indInterp)
% Inputs:
%   toolPt          (3,:) double    the list of tool path points
%   toolNorm        (3,:) double    the list of the spindle direction
%   toolCutDir      (3,:) double    the list of the cutting direction
% Outputs:
%   toolPtInterp    (3,1) double    the current point on the tool path
%   toolNormInterp  (3,1) double    the spindle direction of the current point
%   toolCutDirInterp(3,1) double    the cutting direction of the current point
%
% [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp(
%   toolPt1,toolNorm1,toolCutDir1,toolU1,toolPt2,toolNorm2,toolCutDir2,toolU2,toolUInterp)
% Inputs:
%   toolPt1         (3,1) double    the 1st point on the tool path
%   toolNorm1       (3,1) double    the spindle direction of the 1st point
%   toolCutDir1     (3,1) double    the cutting direction of the 1st point
%   toolU1          (1,1) double    the parameter of the 1st tool path
%   toolPt2         (3,1) double    the 2nd point on the tool path
%   toolNorm2       (3,1) double    the spindle direction of the 2nd point
%   toolCutDir2     (3,1) double    the cutting direction of the 2nd point
%   toolU2          (1,1) double    the parameter of the 2nd tool path
%   toolUInterp     (1,1) double    the parameter of the current tool path
% Outputs: same as the above

switch nargin
    case 8
        toolPt1 = toolPt(:,varargin{1}); % the outer tool point
        toolPt2 = toolPt(:,varargin{2}); % the 1st inner tool point
        toolPt3 = toolPt(:,varargin{3}); % the 2nd inner tool point
        toolNorm1 = toolNorm(:,varargin{1});
        toolNorm2 = toolNorm(:,varargin{2});
        toolNorm3 = toolNorm(:,varargin{3});
        toolCutDir1 = toolCutDir(:,varargin{1});
        toolCutDir2 = toolCutDir(:,varargin{2});
        toolCutDir3 = toolCutDir(:,varargin{3});
        toolU1 = toolU(varargin{1});
        toolU2 = toolU(varargin{2});
        toolUInterp = interpU;
    case 12
        toolPt1 = toolPt;
        toolPt2 = varargin{1};
        toolPt3 = varargin{5};
        toolNorm1 = toolNorm;
        toolNorm2 = varargin{2};
        toolNorm3 = varargin{6};
        toolCutDir1 = toolCutDir;
        toolCutDir2 = varargin{3};
        toolCutDir3 = varargin{7};
        toolU1 = toolU;
        toolU2 = varargin{4};
        toolUInterp = interpU;
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

%% projection
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

% now toolPtInterp should be the interpolation between the toolPt1 and the
% toolPtProj

%% interpolation
u = (toolU1 - toolUInterp)/(toolU1 - toolU2); % p.s. toolU2 = toolU3 = toolUProj
toolPtInterp = toolPt1 - u*(toolPt1 - toolPtProj);
Rot12 = axesRot(toolNorm1,toolCutDir1,toolNormProj,toolCutDirProj,'zx');
quat12 = rotm2quat(Rot12);
quatInterp = slerp([1,0,0,0],quat12,u);
toolNormInterp = quat2rotm(quatInterp)*toolNorm1;
toolNormInterp = toolNormInterp./norm(toolNormInterp);
toolCutDirInterp = quat2rotm(quatInterp)*toolCutDir1;
toolCutDirInterp = toolCutDirInterp./norm(toolCutDirInterp);

end