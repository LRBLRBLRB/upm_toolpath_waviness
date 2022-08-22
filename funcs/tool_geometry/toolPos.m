function [toolPos,toolCutDirect] = toolPos(toolEdge,surfPt,surfNorm,surfDirect)
% usage: [toolPos,toolNorm] = toolPos(toolEdge,surfPt,surfNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface, with the tool axis unchanged.
% Notice: now the tool has not been tilted
% Notice: 刀具现在是固定住朝向的，即认为刀轴方形恒为[0;0;1]不变。
% Input:
%   toolEdge (structure) the edge model of the tool
%   surfPt (3,1) pose of the tool tip (the surface as well)
%   surfNorm (3,1) the normal vector of the surface
% Outputs:
%   toolPos (3,1) the position of the tool center

arguments
    toolEdge
    surfPt (3,:)
    surfNorm (3,:)
    surfDirect (3,:)
end

%% method: rotation transform
% adjust the normal vector of the tool
if nnz(toolEdge.toolEdgeNorm - [0;0;1])
    toolRot = vecRot(toolEdge.toolEdgeNorm,[0;0;1]);
    toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);
end

% rotate the cutting direction and orientation of the tool
toolRot = vecRot(toolEdge.toolDirect,surfDirect);
toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);

% find the contact point on the tool edge
surfToolAng = vecAng(surfNorm,toolEdge.toolEdgeNorm,1);
if abs(surfToolAng) > toolEdge.includedAngle/2
    log = find(abs(surfToolAng) > toolEdge.includedAngle/2);
    error("tool path error: collision between tool and workpiece will happen at point %d",log);
end

% calculate the actual contact point on the tool edge
[~,toolContactPt] = bSplineParam(toolEdge.toolBform,surfToolAng,1e-3, ...
    "Type",'OpenAngle',"IncludedAng",toolEdge.includedAngle);
toolPos = toolEdge.center + surfPt - toolContactPt;
toolCutDirect = toolEdge.toolDirect;

%% method

% a faster method to calculate the angle of two rigid is needed.


end