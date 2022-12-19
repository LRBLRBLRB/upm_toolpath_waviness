function dist = dist2Plane(Pt,varargin)
% to figure out the projection of the given point on the given plane
%
% usage:
%
% dist = dist2Plane(Pt,param)
%   the plane is described with parameters of the general equation
%   Pt          (3,:)   the point that will be projected
%   param       (4,1)   the parameters of the equation, i.e., A, B, C & D,
%                       respectively
%   dist        (1,1)   the distance from the point Pt to the plane
%
% dist = dist2Plane(Pt,surfPt,surfNorm)
%   the plane is described with a point on it and its normal vector
%   Pt          (3,:)   the point that will be projected
%   surfPt      (3,1)   the point on the plane
%   surfNorm    (3,1)   the normal vector of the plane
%   dist        (1,1)   the distance from the point Pt to the plane

switch nargin
    case 2
        param = varargin{1};
        abc = param(1:3);
        d = param(4);
        dist = (abc'*Pt + d)./sqrt(sum(abc.^2));
    case 3
        surfPt = varargin{1};
        surfNorm = varargin{2};
        ptNum = size(Pt,2);
        if ptNum == 1
            dist = dot(Pt - surfPt,surfNorm)./norm(surfNorm);
        else
            dist = dot(Pt - ndgrid(surfPt,1:ptNum),ndgrid(surfNorm,1:ptNum),1)./norm(surfNorm);
        end
    otherwise
        error('invalid input format. Only 2 or 3 inputs are allowed.');
end

end