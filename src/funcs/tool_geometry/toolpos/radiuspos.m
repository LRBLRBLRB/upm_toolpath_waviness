function [toolPathPt,toolQuat,toolContactU,curvePt] = radiuspos( ...
    curveFunc,curveFx,radius,toolPathPt,toolPathNorm,toolPathFeed,options)
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
    curveFx function_handle
    radius double
    toolPathPt (3,1)
    toolPathNorm (3,1) = [0;0;-1]
    toolPathFeed (3,1) = [0;-1;0]
    options.iniDisplay {mustBeMember(options.iniDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
    options.useParallel logical = false
end

%% Initial value of iteration
    function F1 = ini2solve(x,toolpathxy)
        % the x and y coordinates of the toolpath has been given by the
        % input, while tooltip's x & y and the z value of the toolpath
        % need to be solved. so x = [tooltip_y; toolpath_z].
        curveNorm = [curveFx(x(1));-1];
        curveNorm = curveNorm./norm(curveNorm);
        F1(1) = toolpathxy - x(1) + radius*curveNorm(1);
        % F1(2) = toolpathx(2) - x(2) + R*curveNorm(1);
        F1(2) = x(2) - curveFunc(x(1)) + radius*curveNorm(2);
    end

optimOpts1 = optimoptions('fsolve', ...
    'Display',options.iniDisplay, ...
    'UseParallel',options.useParallel);

curvePt = fsolve(@(x)ini2solve(x,toolPathPt(1)), ...
    toolPathPt(2:3),optimOpts1);
toolPathPt(end) = curvePt(end);
curvePt(2) = 0;
curvePt(3) = curveFunc(curvePt(1));

%% rotation transform
% Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
toolRot = axesRot([0;-1;0],[0;0;-1], ...
        toolPathFeed,toolPathNorm,'yz');
toolQuat = rotm2quat(toolRot);

curveNorm0 = [curveFx(curvePt(1));-1];
curveNorm0 = curveNorm0./norm(curveNorm0);
toolContactU = atan2(curveNorm0(end),curveNorm0(1));
end