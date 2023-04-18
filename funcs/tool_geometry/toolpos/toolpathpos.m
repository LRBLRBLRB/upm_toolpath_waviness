function [toolQuad,toolPathPt,toolContactU,surfPt,varargout] = toolpathpos(toolEdge, ...
    toolPathPt,toolPathNorm,toolPathDir,surfFunc,surfNormFunc,options)
% Description:
%   Solve the tool pose and the tool tip pose using the tangent relationship 
%   of tool tip and surface, with the tool axis unchanged.
%   The x & y coordinates of the tool path is given, while the z coordinate
%   as well as the surface contact point will be figured out.
% 
% Usage:
%   [toolQuad,toolPathPt,surfPt] = surfPos(toolEdge,toolPathPt, ...
%       toolPathNorm,toolPathDir,surfFunc,surfNorm)
% Inputs:
%   toolEdge (structure) the edge model of the tool
%   toolPathPt (3,1)    pose of the tool tip (the surface as well)+
%   toolPathNorm (3,1)  the normal vector at the tool path point
%   toolPathDir (3,1)   the cutting direction vector at the tool path point
%   surfFunc     function_handle    surface function
%   surfNormFunc function_handle    surface normal vector
% Outputs:
%   toolQuad    (1,4)   the quaternion of the tool
%   toolPathPt  (3,1)   the position of the tool center
%   surfPt      (3,1)   the position of the surface contact point

arguments
    toolEdge
    toolPathPt (3,1)
    toolPathNorm (3,1)
    toolPathDir (3,1)
    surfFunc function_handle
    surfNormFunc function_handle
    options.iniDisplay {mustBeMember(options.iniDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
    options.finalDisplay {mustBeMember(options.finalDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
    options.finalFuncTol {mustBePositive} = 1e-3
    options.finalStepTol {mustBePositive} = 1e-3
    options.uQTol {mustBePositive} = 1e-3
    options.useParallel logical = false
end

R = toolEdge.radius;

%% Initial value of iteration
    function F1 = ini2solve(x,toolPathXY)
        % Notice: x(3) actually equals to the z value of the tool
        % path point!!!
        surfNorm = surfNormFunc(toolPathXY(1),toolPathXY(2));
        F1(1) = toolPathXY(1) - x(1) + R*surfNorm(1);
        F1(2) = toolPathXY(2) - x(2) + R*surfNorm(2);
        F1(3) = x(3) - surfFunc(x(1),x(2)) + R*surfNorm(3);
    end

optimOpts1 = optimoptions('fsolve', ...
    'Display',options.iniDisplay, ...
    'UseParallel',options.useParallel);

surfPt = fsolve(@(toolTip)ini2solve(toolTip,toolPathPt(1:2)), ...
    [toolPathPt(1:2);surfFunc(toolPathPt(1),toolPathPt(2))], ...
    optimOpts1);
toolPathPt(3) = surfPt(3);

%% rotation transform
% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
toolRot = axesRot(toolEdge.toolEdgeNorm,toolEdge.cutDirect, ...
        toolPathNorm,toolPathDir,'');
toolEdgeRotate = toolRigid(toolEdge,toolRot,[0;0;0]);
toolQuad = rotm2quat(toolRot);

%% find the contact point on the tool edge
    function F2 = final2solve(x,toolPathXY)
        % surfNorm = surfNormFunc(toolPathXY(1),toolPathXY(2));
        toolEdge1 = toolRigid(toolEdgeRotate,eye(3),[toolPathXY;x]);
        toolEdgeList = fnval(toolEdge1.toolBform,uQ);
        toolEdgeDist = dist2surf(toolEdgeList,surfFunc, ...
            'CalculateType','Lagrange-Multiplier', ...
            'DisplayType','none');
        F2 = min(toolEdgeDist);
    end

uQ = linspace(0,1,ceil(1/options.uQTol));
optimOpts2 = optimoptions('fsolve', ...
    'Display',options.finalDisplay, ...
    'FunctionTolerance',options.finalFuncTol, ...
    'StepTolerance',options.finalStepTol, ...
    'UseParallel',options.useParallel);
toolPathZ = fsolve(@(x)final2solve(x,toolPathPt(1:2)),toolPathPt(3),optimOpts2);
toolPathPt(3) = toolPathZ;

surfNorm0 = surfNormFunc(toolPathPt(1),toolPathPt(2));
toolEdgeTrans = toolRigid(toolEdgeRotate,eye(3),toolPathPt);
toolEdgeList0 = fnval(toolEdgeTrans.toolBform,uQ);
[surfPtDistList,surfPtList] = dist2surf(toolEdgeList0,surfFunc, ...
    'CalculateType','Lagrange-Multiplier');
[varargout{1},Ind] = min(surfPtDistList);
surfPt = surfPtList(:,Ind);
[toolContactU,~,~] = toolPtInv(toolEdgeRotate.toolBform,surfNorm0,1e-3, ...
    "Type",'TangentPlane',"Radius",toolEdgeRotate.radius);
% toolContactU = 0.5;
end