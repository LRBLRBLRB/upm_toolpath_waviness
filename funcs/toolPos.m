function [toolPos] = toolPos(toolEdge,surfPt,surfNorm,toolNorm)
% usage: [toolPos,toolNorm] = toolPos(toolTip,toolNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface
% Input:
%   toolEdge (structure) the edge model of the tool
%   surfPt (1,3) pose of the tool tip (the surface as well)
%   surfNorm (1,3) the normal vector of the surface
%   toolNorm (1,3) the normalized normal vector of the tool
% Outputs:
%   toolPos (1,3) the position of the tool center

arguments
    toolEdge 
    surfPt (:,3)
    surfNorm (:,3)
    toolNorm (:,3)
end


%% method: rotation transform
c = toolEdge.center;
% r = toolEdge.radius;
toolPt = toolEdge.toolPt; % (:,3) position of tool edge
tmpNorm = [0,0,1];
vecAng = atan2(norm(cross(surfNorm,toolNorm)),dot(surfNorm,toolNorm))*180/pi;
if abs(vecAng) > toolEdge.includedAngle/2
    log = find(abs(vecAng) > toolEdge.includedAngle/2);
    error("tool path error: collision between tool and workpiece will happen at point %d",log);
end
surfR = vecRot(surfNorm,tmpNorm);
toolNorm = toolNorm*surfR;
toolR = vecRot([0,1,0],toolNorm);
tmpPt = toolPt*toolR;
tmpC = c*toolR;
[~,minInd] = min(tmpPt(:,3));
tranVec = surfPt - tmpPt(minInd,:);
toolPos = tmpC + tranVec;
% toolPos = tmpC/toolR; % tmpC*inv(R)

%% method

% a faster method to calculate the angle of two rigid is needed.




end