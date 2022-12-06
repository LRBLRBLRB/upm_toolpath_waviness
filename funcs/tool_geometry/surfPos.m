function [surfPt,toolPathPt] = surfPos(toolEdge,toolPathPt,toolPathNorm,toolPathDir,surfFunc,surfNorm)
%   Solve the tool pose and the tool tip pose using the tangent relationship 
%   of tool tip and surface, with the tool axis unchanged.





% rotation transform
% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
toolRot = axesRot(toolEdge.toolEdgeNorm,toolEdge.cutDirect, ...
        toolPathNorm,toolPathDir,'');
toolEdge = toolRigid(toolEdge,toolRot,[0;0;0]);

% find the contact point on the tool edge

[toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
    "Type",'TangentPlane',"Radius",toolEdge.radius);



end