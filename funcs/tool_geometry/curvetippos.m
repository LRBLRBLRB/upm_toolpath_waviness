function [toolPathPt,toolQuat,toolContactU] = curvetippos(toolEdge,curvePt,curveNorm,toolPathNorm,toolPathFeed)
% usage: [toolPos,toolCutDirect,log] = toolPos(toolEdge,surfPt,surfNorm,designDirect,designNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface, with the tool axis unchanged.
% Notice: now the tool has not been tilted
% Notice: 刀具现在是固定住朝向的，即认为刀轴方形恒为[0;0;1]不变。
% Input:
%   toolEdge (structure) the edge model of the tool
%   curvePt (3,1) pose of the tool tip (the surface as well)
%   toolPathNorm (3,1) the normal vector of the toolpath
%   designNorm (3,1)
% Outputs:
%   toolPathPt (3,1) the position of the tool center
%   toolQuat    
%   toolContactU


if size(curvePt,1) ~= 3 || size(toolPathNorm,1) ~= 3
    error('Invalid input: the dimension goes wrong.');
end

% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
toolRot = axesRot(toolEdge.toolDirect,toolEdge.toolEdgeNorm, ...
    toolPathFeed,toolPathNorm,'yz');
toolEdge2 = toolRigid(toolEdge,toolRot,[0;0;0]);

% find the contact point on the tool edge
surfToolAng = vecAng(toolPathNorm,curveNorm,1);
if abs(surfToolAng) > toolEdge2.openAngle/2
%     varargout{1} = true;
%     toolVec = nan(3,1);
%     toolQuad = nan(4,1);
    error("tool path error: collision between tool and workpiece will happen at the point");
end
toolQuat = rotm2quat(toolRot);

% calculate the actual contact point on the tool edge
% [~,toolContactPt] = toolPtInv(toolEdge.toolBform,surfToolAng,1e-3, ...
%     "Type",'PolarAngle',"IncludedAng",toolEdge.openAngle);

curveNorm = toolRot'*curveNorm;
[toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,curveNorm,1e-3, ...
    "Type",'TangentPlane',"Radius",toolEdge.radius);

% surfNorm(1) = 0;
% [toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
%     "Type",'TangentLine',"Radius",toolEdge.radius);
toolPathPt = curvePt - toolRot*toolContactPt;


end