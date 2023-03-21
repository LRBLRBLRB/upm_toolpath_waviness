function [res,peakPt,varargout] = residual2D_multi(s1,s2,eps,varargin)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,peakPt,ind1,ind2] = residual2D_numeric(s1,s2,eps,c1,c2,vec1,vec2,r,method)
%   s1 (2,:) 1st scatters of the tool edge
%   c1 (2,1) 1st center of the tool edge
%   vec1 (2,1) 1st normal vector of the tool edge
%   c2,vec2,s2 are parameters of the other tool edge
%   res (1,1) the residual within the two position
%   interPt (2,1)
%
% [res,peakPt,ind1,ind2] = residual2D_numeric(s1,s2,eps,p1,p2,method)
%   s1 & s2 struct  the B-form struct of the two tool edge
%   eps (1,1)  the discretization of the parameter u
%   p1 & p2 () 
%   res (1,1)  the residual within the two position
%   
%
% [res,peakPt,ind1,ind2] = residual2D_numeric(sp1,sp2,...)
%   all the same except that the sp1 and sp2 remain the B-form spline
%   struct pf the tool edge

%% input pre-processing
switch nargin
    case {5,8}
        method = 'DSearchn';
    case 6
        method = varargin{3};
    case 9
        method = varargin{6};
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

%% solve the intersection point of the two spline curve
switch method
    case 'DSearchn'
        [peakPt,ind1,ind2,epsCross] = pcCrossDN(sp1,sp2);
    case 'BoundingBox'
        [peakPt,ind1,ind2,epsCross] = boundingBox(sp1,sp2);
    otherwise
        warning('Invalid method input! default method ''dsearchn'' will be used');
        [peakPt,ind1,ind2,epsCross] = pcCrossDN(sp1,sp2);
end

%% calculate the residual height on the 2-D plane
switch nargin
    case {5,6}
        cutPt1 = varargin{1};
        cutPt2 = varargin{2};
        if size(cutPt1,1) == 1
            cutPt1 = cutPt1';
            cutPt2 = cutPt2';
        end
        % use the formula of calculating the area of triangle to calculate
        % the height of the triangle, the vertices of which are cutPt1,
        % cutPt2, peakPt.
        res = norm(cross(cutPt1 - peakPt,cutPt2 - peakPt)) ...
            /norm(cutPt2 - cutPt1); 
    case {8,9}
        peakPt(1) = [];
        cen1 = varargin{1};
        cen2 = varargin{2};
        vec1 = varargin{3};
        vec2 = varargin{4};
        r = varargin{5};
        if size(cen1,1) == 1
            cen1 = cen1';
            cen2 = cen2';
            vec1 = vec1';
            vec2 = vec2';
        end
        vec = (vec1 + vec2)/norm(vec1 + vec2);
        cen1 = cen1 + r*vec;
        cen2 = cen2 + r*vec;
        line(1) = cen2(2) - cen1(2);
        line(2) = cen1(1) - cen2(1);
        line(3) = cen2(1)*cen1(2) - cen1(1)*cen2(2);
        res = pt2Line(peakPt,line);
end

varargout{1} = u(ind1);
varargout{2} = u(ind2);
if epsCross >= norm(sp1(:,1) - sp1(:,2))
    % no intersection
    res = inf;
    peakPt = inf(size(sp1,1),1);
end

end