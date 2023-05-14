function [toolPathPt,toolQuat,toolContactU] = radiustippos( ...
    toolData,curvePt,curveNorm,toolPathNorm,toolPathVec,options)
% usage: [toolPos,toolCutDirect,log] = toolPos(toolData,surfPt,surfNorm,designDirect,designNorm)
%   Solve the tool pose using the tangent relationship of tool tip and
%   surface, with the tool axis unchanged.
% Notice: now the tool has not been tilte                                                                                                                 d
% Notice: 刀具现在是固定住朝向的，即认为刀轴方形恒为[0;0;1]不变。
% Input:
%   toolData (structure) the edge model of the tool
%   curvePt (3,1) pose of the tool tip (the surface as well)
%   toolPathNorm (3,1) the normal vector of the toolpath
%   designNorm (3,1)
% Outputs:
%   toolPathPt (3,1) the position of the tool center
%   toolQuat    
%   toolContactU
arguments
    toolData struct
    curvePt (3,1) double
    curveNorm (3,1)
    toolPathNorm
    toolPathVec = []
    options.directionType {mustBeMember(options.directionType, ...
        {'quaternion','norm-cut','norm-feed'})} = 'norm-cut'
end

if size(curvePt,1) ~= 3 || size(toolPathNorm,1) ~= 3
    error('Invalid input: the dimension goes wrong.');
end

switch options.directionType
    case 'norm-cut'
        % Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
        toolRot = axesRot(toolData.toolEdgeNorm,toolData.cutDirect, ...
            toolPathNorm,toolPathVec,'zx');
    case 'norm-feed'
        % Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
        toolRot = axesRot(toolData.toolDirect,toolData.toolEdgeNorm, ...
            toolPathVec,toolPathNorm,'yz');
    case 'quaternion'
        % Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
        toolRot = quat2rotm(toolPathNorm);
end
toolQuat = rotm2quat(toolRot);

% calculate the actual contact point on the tool edge
toolPathPt = curvePt - toolData.toolRadius*curveNorm;
toolContactU = atan2(-curveNorm(end),-curveNorm(1));
end