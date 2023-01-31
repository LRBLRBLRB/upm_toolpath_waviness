function [surfPt,toolPathPt] = surfPos(toolEdge,toolPathPt,toolPathNorm, ...
    toolPathDir,surfFunc,surfNorm)
%   Solve the tool pose and the tool tip pose using the tangent relationship 
%   of tool tip and surface, with the tool axis unchanged.
%
% Input:
%   toolEdge (structure) the edge model of the tool
%   surfPt (3,1) pose of the tool tip (the surface as well)
%   surfNorm (3,1) the normal vector of the surface
%   designDirect (3,1)
%   designNorm (3,1)
% Outputs:
%   toolPos (3,1) the position of the tool center

R = toolEdge.toolRadius;
%% Initial value of iteration
    function F = ini2solve(x,toolPathXY,surfFunc,surfNorm,R)
        % Notice: x(3) actually equals to the z value of the tool
        % path point!!!
        normNorm = surfNorm/norm(surfNorm);
        F(1) = toolPathXY(1) - x(1) + R*normNorm(1);
        F(2) = toolPathXY(2) - x(2) + R*normNorm(1);
        F(3) = x(3) - surfFunc(x(1),x(2)) + R*normNorm(1);
    end

surfPt = fsolve(@(toolTip)ini2solve(toolTip,toolPathPt(1:2),surfFunc,surfNorm,R), ...
    [toolPathPt(1:2);surfFunc(toolPathPt(1),toolPathPt(2))]);
toolPathPt(3) = surfPt(3);
surfPt(3) = surfFunc(surfPt(1),surfPt(2));


% %% rotation transform
% % Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
% toolRot = axesRot(toolEdge.toolEdgeNorm,toolEdge.cutDirect, ...
%         toolPathNorm,toolPathDir,'');
% toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);
% 
% % find the contact point on the tool edge
% while true
%     toolEdge1 = toolRigid(toolEdge,eye(3),toolPathPt);
%     uQ = uLim(1):0.01:uLim(2);
%     QList = fnval(toolEdge1.toolBform,uQ);
%     QDist = dist2surf(QList,surfFunc, ...
%         -surfNorm(1)/surfNorm(3),-surfNorm(2)/surfNorm(3),'Lagrange-Multiplier');
%     
%     % when QDist only one equals to zero, then break!!!
% end





end