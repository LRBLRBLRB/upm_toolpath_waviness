function [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp( ...
    toolPt,toolNorm,toolCutDir,toolU,varargin)
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
    case 7
        toolPt1 = toolPt(:,varargin{1});
        toolPt2 = toolPt(:,varargin{2});
        toolNorm1 = toolNorm(:,varargin{1});
        toolNorm2 = toolNorm(:,varargin{2});
        toolCutDir1 = toolCutDir(:,varargin{1});
        toolCutDir2 = toolCutDir(:,varargin{2});
        toolU1 = toolU(varargin{1});
        toolU2 = toolU(varargin{2});
        toolUInterp = toolU(varargin{3});
    case 10
        toolPt1 = toolPt;
        toolPt2 = varargin{1};
        toolNorm1 = toolNorm;
        toolNorm2 = varargin{2};
        toolCutDir1 = toolCutDir;
        toolCutDir2 = varargin{3};
        toolU1 = toolU;
        toolU2 = toolU(varargin{4});
        toolUInterp = toolU(varargin{5});
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

%%
t = (toolUInterp - toolU1)/(toolU2 - toolU1);
toolPtInterp = toolPt1 + t*(toolPt2 - toolPt1);
Rot12 = axesRot(toolNorm1,toolCutDir1,toolNorm2,toolCutDir2,'zx');
quat12 = rotm2quat(Rot12);
quatInterp = slerp([1,0,0,0],quat12,t);
toolNormInterp = quat2rotm(quatInterp)*toolNorm1;
toolNormInterp = toolNormInterp./norm(toolNormInterp);
toolCutDirInterp = quat2rotm(quatInterp)*toolCutDir1;
toolCutDirInterp = toolCutDirInterp./norm(toolCutDirInterp);

end