function projVec = vecOnPlane(vec,varargin)
% to figure out the projection of the given point on the given plane
%
% usage:
%
% projVec = vecOnPlane(vec,param)
%   the plane is described with parameters of the general equation
%   vec         (3,:)   the vector that will be projected
%   param       (4,1)   the parameters of the equation, i.e., A, B, C & D,
%                       respectively
%   projVec     (3,:)   the projection of the point Pt
%
% projVec = vecOnPlane(vec,surfPt,surfNorm)
%   the plane is described with a point on it and its normal vector
%   vec         (3,:)   the vector that will be projected
%   surfPt      (3,1)   the point on the plane
%   surfNorm    (3,1)   the normal vector of the plane
%   projVec     (3,:)   the projection of the point Pt
%
% Reference: 
%   https://www.cnblogs.com/nobodyzhou/p/6145030.html
%   https://juejin.cn/post/6913106721355661319
switch nargin
    case 2
        param = varargin{1};
        surfNorm = [param(1);param(2);param(3)];
    case 3
        % surfPt = varargin{1};
        surfNorm = varargin{2};
    otherwise
        error('Invalid input format. Only 2 or 3 inputs are allowed.');
end

vecNum = size(vec,2);
projVec = vec - (dot(vec,ndgrid(surfNorm,1:vecNum),1) ...
    ./norm(surfNorm)).*surfNorm;

end