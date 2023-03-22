function [res,u1,u2] = residual2D_multi(s1,s2,eps,surfFunc,varargin)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,peakPt,ind1,ind2] = residual2D_numeric(s1,s2,eps,p1,p2,method)
%   s1 & s2 struct  the B-form struct of the two tool edge
%   eps (1,1)  the discretization of the parameter u
%   p1 & p2 () 
%   res (1,1)  the residual within the two position
%   
%
% [res,peakPt,ind1,ind2] = residual2D_numeric(sp1,sp2,...)
%   all the same except that the sp1 and sp2 remain the B-form spline
%   struct of the tool edge

%% input pre-processing
switch nargin
    case 4
        method = 'DSearchn';
    case 5
        method = varargin{1};
    otherwise
        error('Too many parameters input!');
end

%% to ensure that the form of the input point remains the 3*1 vector
if isstruct(s1)
    u = 0:eps:1;
    sp1 = fnval(s1,u);
    sp2 = fnval(s2,u);
elseif size(s1,1) == 1
    sp1 = s1';
    sp2 = s2';
end
dim = size(sp1,1);

%% solve the intersection point of the two spline curve
switch method
    case 'pt2surf'
        if dim == 2
            dist1 = dist2curve(sp1,surfFunc);
            dist2 = dist2curve(sp2,surfFunc);
        else
            dist1 = dist2surf(sp1,surfFunc);
            dist2 = dist2surf(sp2,surfFunc);
        end
%         u1 = dist1 <= dist2; % the closer one will be selected
%         u2 = dist1 > dist2;
          cmp = dist1 <= dist2;
    case 'BoundingBox'
        [peakPt,ind1,ind2,epsCross] = boundingBox(sp1,sp2);
    otherwise
        error('Invalid method input!');
end

end