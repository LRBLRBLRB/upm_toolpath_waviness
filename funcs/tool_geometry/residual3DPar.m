function [res,uLim1,uLim2] = residual3DPar(toolPt,toolNorm,toolCutDirect,sp,ind1,ind2,ind3)
% To calculate the residual height among the adjacent tool points.
% Different from the function "residual3D", this function is more siutable
% for parallel computing.
%
% Usage:
%
% [res,distInd,uLim1,uLim2] = residual3D(toolPt,toolNorm,toolCutDirect,sp,
% ind1,ind2,ind3)
% Inputs:
%   toolPt          (3,:) double    the current point on the tool path
%   toolNorm        (3,:) double    the spindle direction of the current point
%   toolCutDirect   (3,:) double    the cutting direction of the current point
%   sp              (1,1) struct    the B-form spline struct of the tool edge
%   toolPt2         (3,1) double    the closest point on the tool path
%   toolNorm2       (3,1) double    the spindle direction of "toolPt2"
% Outputs:
%   res             (1,1) double    the residual height of the two point
%   distInd         (1,1) double    the index of the one cloest to the
%                                   current point
%   uLim1           (1,2) double    the interval of the parameter for the
%                                   current tool point
%   uLim2           (1,2) double    the interval of the parameter for the
%                                   cloest tool point

% to project the closest point and its tool edge to the current one
toolPtProj = projectionOnPlane(toolPt(:,ind1),toolPt(:,ind2),toolCutDirect(:,ind1));
% Van der Waals: 主流的做法是对位置插补 然后刀轴是跟随运动；位置已经运动的长度/总长度=已经转过的角度/总角度
% 是否可以直接用向量插值、或者四元数来代替上述公式里的角度？
t = norm(toolPtProj - toolPt(:,ind2))/norm(toolPt(:,ind3) - toolPt(:,ind2));
toolNormProj = (1-t)*toolNorm(:,ind2) + t*toolNorm(:,ind3);

[res,~,uLim1,uLim2] = residual2D( ...
    toolPt(:,ind1),toolPtProj,toolNorm(:,ind1),toolNormProj,sp);

end