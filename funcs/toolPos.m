function [toolPos] = toolPos(toolEdge,surfPt,surfNorm,surfDirect)
% usage: [toolPos,toolNorm] = toolPos(toolEdge,surfPt,surfNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface, with the tool axis unchanged.
% Notice: now the tool has not been tilted
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
c = toolEdge.center;
% r = toolEdge.radius;
toolPt = toolEdge.toolPt; % (3,:) position of tool edge
toolNorm = toolEdge.toolEdgeNorm; % (3,:) normal vector of the tool edge

% adjust the normal vector of the tool
if nnz(toolNorm - [0;0;1])
    toolRot = vecRot(toolNorm,[0;0;1]);
    toolPt = toolRot*toolPt;
end

surfToolAng = vecAng(surfNorm,toolNorm,1);
if abs(surfToolAng) > toolEdge.includedAngle/2
    log = find(abs(surfToolAng) > toolEdge.includedAngle/2);
    error("tool path error: collision between tool and workpiece will happen at point %d",log);
end
surfRot = vecRot(surfNorm,toolNorm);
toolNorm = surfRot*toolNorm;
tmpC = toolRot*c;
[~,minInd] = min(toolPt(:,3));
tranVec = surfPt - toolPt(minInd,:);
toolPos = tmpC + tranVec;
% toolPos = tmpC/toolR; % tmpC*inv(R)

%% method

% a faster method to calculate the angle of two rigid is needed.




end