function projPt = projectionOnPlane(Pt,varargin)
% to figure out the projection of the given point on the given plane
%
% usage:
%
% projPt = projectionOnPlane(Pt,param)
%   the plane is described with parameters of the general equation
%   Pt          (3,:)   the point that will be projected
%   param       (4,1)   the parameters of the equation, i.e., A, B, C & D,
%                       respectively
%   projPt      (3,:)   the projection of the point Pt
%
% projPt = projectionOnPlane(Pt,surfPt,surfNorm)
%   the plane is described with a point on it and its normal vector
%   Pt          (3,:)   the point that will be projected
%   surfPt      (3,1)   the point on the plane
%   surfNorm    (3,1)   the normal vector of the plane
%   projPt      (3,:)   the projection of the point Pt
%
% Reference: 
%   https://www.cnblogs.com/nobodyzhou/p/6145030.html
%   https://juejin.cn/post/6913106721355661319
switch nargin
    case 1
        param = varargin{1};
        abc = param(1:3);
        d = param(4);
        t = (abc'*Pt + d)./sum(abc.^2);
        projPt = Pt - t.*abc;
    case 2
        surfPt = varargin{1};
        surfNorm = varargin{2};
        ptNum = size(Pt,2);
        if ptNum == 1
            projPt = Pt ...
                - dot(Pt - surfPt,surfNorm)./norm(surfNorm).*surfNorm;
        else
            projPt = Pt - ...
                dot(Pt - ndgrid(surfPt,1:ptNum),ndgrid(surfNorm,1:ptNum),1) ...
                ./norm(surfNorm).*surfNorm;
        end
    otherwise
        error('invalid input format. Only 2 or 3 inputs are allowed.');
end

end