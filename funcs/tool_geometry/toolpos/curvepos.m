function [toolPathPt,toolQuat,toolContactU,curvePt,varargout] = curvepos( ...
    curveFunc,curveFy,toolEdge,toolPathPt,toolPathNorm,toolPathFeed,options)
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
%   curveFunc
%   curveFy 
%   toolEdge struct the edge model of the tool
%   toolPathPt (3,1) pose of the tool tip (the surface as well)
%   toolPathNorm (3,1) the normal vector of the tool, whose default value 
%       is [0;0;-1], i.e., the tool points towards the workpiece
%   toolPathFeed (3,1) the feed direction of the tool. Ony [0;0;1] and
%       [0;0;-1] can be chosen, which represent Center-to-edge and
%       Edge-to-center feed, respectively
%   options
% Outputs:
%   toolPathPt  (3,1)   the position of the tool center
%   toolQuat    (1,4)   the quaternion of the tool
%   toolContactU (1,1)
%   surfPt      (3,1)   the position of the surface contact point

arguments
    curveFunc function_handle
    curveFy function_handle
    toolEdge struct
    toolPathPt (3,1)
    toolPathNorm (3,1) = [0;0;-1]
    toolPathFeed (3,1) = [0;-1;0]
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
    function F1 = ini2solve(x,toolpathxy)
        % the x and y coordinates of the toolpath has been given by the
        % input, while tooltip's x & y and the z value of the toolpath
        % need to be solved. so x = [tooltip_y; toolpath_z].
        curveNorm = [curveFy(x(1));-1];
        curveNorm = curveNorm./norm(curveNorm);
        F1(1) = toolpathxy - x(1) + R*curveNorm(1);
        % F1(2) = toolpathx(2) - x(2) + R*curveNorm(1);
        F1(2) = x(2) - curveFunc(x(1)) + R*curveNorm(2);
    end

optimOpts1 = optimoptions('fsolve', ...
    'Display',options.iniDisplay, ...
    'UseParallel',options.useParallel);

curvePt = fsolve(@(x)ini2solve(x,toolPathPt(2)), ...
    toolPathPt(2:3),optimOpts1);
toolPathPt(end) = curvePt(end);
% curvePt(end) = curveFunc(curvePt(1));

%% rotation transform
% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
toolRot = axesRot(toolEdge.toolDirect,toolEdge.toolEdgeNorm, ...
        toolPathFeed,toolPathNorm,'yz');
toolEdgeRotate = toolRigid(toolEdge,toolRot,[0;0;0]);
toolQuat = rotm2quat(toolRot);

%% find the contact point on the tool edge
    function F2 = final2solve(x,toolpathx)
        % surfNorm = surfNormFunc(toolPathXY(1),toolPathXY(2));
        toolEdge1 = toolRigid(toolEdgeRotate,eye(3),[toolpathx;x]);
        toolEdgeList = fnval(toolEdge1.toolBform,uQ);
        toolEdgeDist = dist2curve(toolEdgeList(2:3,:),curveFunc,curveFy, ...
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
toolPathZ = fsolve(@(x)final2solve(x,toolPathPt(1:2)),toolPathPt(end),optimOpts2);
toolPathPt(end) = toolPathZ;

curveNorm0 = [0;curveFy(curvePt(1));-1];
curveNorm0 = curveNorm0./norm(curveNorm0);
toolEdgeTrans = toolRigid(toolEdgeRotate,eye(3),toolPathPt);
toolEdgeList0 = fnval(toolEdgeTrans.toolBform,uQ);
[surfPtDistList,surfPtList] = dist2curve(toolEdgeList0(2:3,:), ...
    curveFunc,curveFy,'CalculateType','Lagrange-Multiplier');
[varargout{1},Ind] = min(surfPtDistList);
curvePt = [0;surfPtList(:,Ind)];

if abs(ini2solve([curvePt(2);toolPathPt(3)],toolPathPt(2))) > toolEdge.radius
    error('The curvature of the present on the surface is too large to be machined by the tool.');
end

[toolContactU,~,~] = toolPtInv(toolEdgeRotate.toolBform,curveNorm0,1e-3, ...
    "Type",'TangentPlane',"Radius",toolEdgeRotate.radius);
% toolContactU = 0.5;
end