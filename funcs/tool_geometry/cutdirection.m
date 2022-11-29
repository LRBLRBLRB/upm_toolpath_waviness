function toolDirect = cutdirection(varargin)
% to solve the turning direction of the tool edge, provide the tool
% direction of the tool edge at the exact points that are given.
% usage: 
%   toolDirect = cutDirection(surfPt,surfCenter)
% Inputs:
%   surfPt (3,:) the surface point array
%   surfCenter (3,1) the rotation center or the point on the rotation axis
% Outputs:
%   toolDirect (3,:) the cutting direction of each point in the array
%
%   toolDirect = cutDirection(surfPt)

switch nargin
    case 2
        surfPt = varargin{1};
        surfCenter = varargin{2};
        % 思路一：对于回转体的同心圆算法（concentric），可以直接求回转轴指向曲面点的向量，然后投影到回转轴垂直的表面，再取x,y后右转90度。
        Rot = rotz(pi/2);
        toolDirect = Rot*(surfPt - surfCenter);
        toolDirect(3,:) = 0;
        toolDirect = toolDirect ./ vecnorm(toolDirect,2,1);
    case 1
        surfPt = varargin{1};
        % 思路二：直接拿surfMesh中的下一个点和上一个点的坐标，然后算方向：
        surfPt = [surfPt(:,1),surfPt,surfPt(:,1)];
        toolDirect = 0.5*(surfPt(:,3:end) - surfPt(:,1:end-2));
        toolDirect = toolDirect ./ vecnorm(toolDirect,2,1);
    otherwise
        error("Invalid input");
end

end