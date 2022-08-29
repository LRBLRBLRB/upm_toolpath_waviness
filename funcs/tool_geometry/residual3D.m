function [res,distInd,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,sp,toolPath,toolPathNorm)
% to calculate the residual height among the adjacent tool points.
%
% Usage:
%
% [res,distInd,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,sp,
% toolPath,toolPathNorm)
% Inputs:
%   toolPt          (3,1) double    the current point on the tool path
%   toolNorm        (3,1) double    the spindle direction of the current point
%   toolCutDirect   (3,1) double    the cutting direction of the current point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolPath        (3,:) double    the whole tool path
%   toolPathNorm    (3,:) double    the spindle direction of the whole tool path
% Outputs:
%   res             (1,1) double    the residual height of the two point
%   distInd         (1,1) double    the index of the one cloest to the
%                                   current point
%   uLim1           (1,2) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,2) double    the interval of the parameter for the
%                                   cloest tool point

ptNum = size(toolPath,2);
% to obtain those within the vicinity of the current tool point
% isAdj = (vecnorm(toolPath - ndgrid(toolPt,1:ptNum),2,1) - 2*toolRadius) <= 0;

% to find those in the tool path which is cloest to the current tool point
% 有问题：需要把两边的最近点都分别找到！！不然会出现找漏
toolPathProj = toolPath - ...
    dot(toolPath - ndgrid(toolPt,1:ptNum),ndgrid(toolCutDirect,1:ptNum),1) ...
    ./norm(toolCutDirect).*toolCutDirect;
dist = vecnorm(toolPathProj - ndgrid(toolPt,1:ptNum),2,1);
[~,distInd] = min(dist);

[res,~,uLim1,uLim2] = residual2D(toolPt,toolPath(distInd,:),toolNorm,toolPathNorm(distInd,:),sp);

end