function [toolPos,toolCutDirect,log] = toolPos(toolEdge,surfPt,surfNorm,surfDirect)
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
log = 0;
if nnz(toolEdge.toolEdgeNorm - [0;0;1])
    toolRot = vecRot(toolEdge.toolEdgeNorm,[0;0;1]);
    toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);
end

% rotate the cutting direction and orientation of the tool
toolRot = vecRot(toolEdge.cutDirect,surfDirect);
toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);

% find the contact point on the tool edge
surfToolAng = vecAng(surfNorm,toolEdge.toolDirect,1);
if abs(surfToolAng - pi/2) > toolEdge.includedAngle/2
    warning("tool path error: collision between tool and workpiece will happen at the point");
    log = 1;
    toolPos = nan(3,1);
    toolCutDirect = nan(3,1);
    return;
end

% calculate the actual contact point on the tool edge
[~,toolContactPt] = bSplinePtInv(toolEdge.toolBform,surfToolAng,1e-3, ...
    "Type",'PolarAngle',"IncludedAng",toolEdge.includedAngle);
toolContactPt = [toolContactPt(1);0;toolContactPt(2)];
toolContactPt = toolRot*toolContactPt;
toolPos = toolEdge.center + surfPt - toolContactPt;% 检查！！！！！！！！！！！！！！！！！！！！！
toolCutDirect = toolEdge.cutDirect;

%% method

% a faster method to calculate the angle of two rigid is needed.


end