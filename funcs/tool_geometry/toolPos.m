function [toolQuad,toolVec,isCollision] = toolPos(toolEdge,surfPt,surfNorm,designDirect,designNorm)
% usage: [toolPos,toolCutDirect,log] = toolPos(toolEdge,surfPt,surfNorm,designDirect,designNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface, with the tool axis unchanged.
% Notice: now the tool has not been tilted
% Notice: 刀具现在是固定住朝向的，即认为刀轴方形恒为[0;0;1]不变。
% Input:
%   toolEdge (structure) the edge model of the tool
%   surfPt (3,1) pose of the tool tip (the surface as well)
%   surfNorm (3,1) the normal vector of the surface
%   designDirect (3,1)
%   designNorm (3,1)
% Outputs:
%   toolPos (3,1) the position of the tool center

arguments
    toolEdge
    surfPt (3,:)
    surfNorm (3,:)
    designDirect (3,:)
    designNorm (3,:)
end

%% method: rotation transform
% adjust the normal vector of the tool
isCollision = false;
if nnz(toolEdge.toolEdgeNorm - [0;0;1])
    toolRot1 = vecRot(toolEdge.toolEdgeNorm,[0;0;1]);
    toolEdge = toolRigid(toolEdge,toolRot1,[0;0;0]);
end

% rotate the cutting direction and orientation of the tool
toolRot1 = vecRot(toolEdge.cutDirect,designDirect);
toolEdge = toolRigid(toolEdge,toolRot1,[0;0;0]);
toolRot2 = vecRot(toolEdge.toolEdgeNorm,designNorm);
toolEdge = toolRigid(toolEdge,toolRot2,[0;0;0]);

% find the contact point on the tool edge
surfToolAng = vecAng(surfNorm,toolEdge.toolDirect,1);
if abs(surfToolAng - pi/2) > toolEdge.includedAngle/2
    warning("tool path error: collision between tool and workpiece will happen at the point");
    isCollision = true;
    toolVec = nan(3,1);
    toolQuad = nan(4,1);
    return;
end
toolRot = toolRot2*toolRot1;
toolQuad = rotm2quat(toolRot);

% calculate the actual contact point on the tool edge
[~,toolContactPt] = toolPtInv(toolEdge.toolBform,surfToolAng,1e-3, ...
    "Type",'PolarAngle',"IncludedAng",toolEdge.includedAngle);
% [~,toolContactPt] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
%     "Type",'TangentPlane',"IncludedAng",toolEdge.includedAngle);
toolContactPt = [0;toolContactPt(1);toolContactPt(2)]; % the tool plane remains the yOz plane
toolVec = surfPt - toolRot*toolContactPt;

%% method
% a faster method to calculate the angle of two rigid is probably needed.

end