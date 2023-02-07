function [toolQuad,toolVec,toolContactU,varargout] = toolPos( ...
    toolEdge,surfPt,surfNorm,toolPathNorm,varargin)
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

if nargin == 4 % 2D
    if size(surfPt,1) ~= 2 || size(surfNorm,1) ~= 2 || size(toolPathNorm,1) ~= 2
        error('Invalid input: the dimension goes wrong.');
    end

    varargout{1} = false;
    toolRot = rotz(vecAng([0;1],toolPathNorm,1));
    toolRot = toolRot(1:2,1:2);
    toolEdge.center = toolRot*toolEdge.center;
    toolEdge.toolEdgeNorm = toolRot*toolEdge.toolEdgeNorm;
    toolEdge.toolDirect = toolRot*toolEdge.toolDirect;

    surfToolAng = vecAng(surfNorm,toolEdge.toolDirect,1);
    if abs(surfToolAng - pi/2) > toolEdge.openAngle/2
    %     varargout{1} = true;
    %     toolVec = nan(3,1);
    %     toolQuad = nan(4,1);
        error("tool path error: collision between tool and workpiece will happen at the point");
    end

    surfNorm = toolRot'*surfNorm;
    [toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
        "Type",'TangentPlane',"Radius",toolEdge.radius);
    toolVec = surfPt - toolRot*toolContactPt;

else % 3D
    if size(surfPt,1) ~= 3 || size(surfNorm,1) ~= 3 || size(toolPathNorm,1) ~= 3
        error('Invalid input: the dimension goes wrong.');
    end

    toolPathDir = varargin{1};
    % adjust the normal vector of the tool
    varargout{1} = false;
    if nnz(toolEdge.toolEdgeNorm - [0;0;1])
        toolRot1 = vecRot(toolEdge.toolEdgeNorm,toolPathNorm);
        toolEdge1 = toolRigid(toolEdge,toolRot1,[0;0;0]);
    end

    % rotate the cutting direction and orientation of the tool
    % toolRot1 = vecRot(toolEdge.cutDirect,designDirect);
    % toolEdge = toolRigid(toolEdge,toolRot1,[0;0;0]);
    % toolRot2 = vecRot(toolEdge.toolEdgeNorm,designNorm);
    % toolEdge = toolRigid(toolEdge,toolRot2,[0;0;0]);
    % toolRot = toolRot2*toolRot1;

    % Def: toolNorm = [0;0;1]; cutDirect = [1;0;0];
    toolRot = axesRot(toolEdge.toolEdgeNorm,toolEdge.cutDirect, ...
        toolPathNorm,toolPathDir,'');
    toolEdge2 = toolRigid(toolEdge,toolRot,[0;0;0]);

    % find the contact point on the tool edge
    surfToolAng = vecAng(surfNorm,toolEdge2.toolDirect,1);
    if abs(surfToolAng - pi/2) > toolEdge2.openAngle/2
    %     varargout{1} = true;
    %     toolVec = nan(3,1);
    %     toolQuad = nan(4,1);
        error("tool path error: collision between tool and workpiece will happen at the point");
    end
    toolQuad = rotm2quat(toolRot);

    % calculate the actual contact point on the tool edge
    % [~,toolContactPt] = toolPtInv(toolEdge.toolBform,surfToolAng,1e-3, ...
    %     "Type",'PolarAngle',"IncludedAng",toolEdge.openAngle);
    
    surfNorm = toolRot'*surfNorm;
    [toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
        "Type",'TangentPlane',"Radius",toolEdge.radius);

    % surfNorm(1) = 0;
    % [toolContactU,toolContactPt,~] = toolPtInv(toolEdge.toolBform,surfNorm,1e-3, ...
    %     "Type",'TangentLine',"Radius",toolEdge.radius);
    toolVec = surfPt - toolRot*toolContactPt;
end


end