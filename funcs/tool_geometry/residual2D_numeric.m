function [res,peakPt,varargout] = residual2D_numeric(s1,s2,cutPt1,cutPt2,options)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,peakPt,ind1,ind2] = residual2D_numeric(c1,c2,vec1,vec2,s1,s2,method)
%   vec1 (2,1) 1st normal vector of the tool edge
%   s1 (2,:) 1st scatters of the tool edge
%   c2,vec2,s2 are parameters of the other tool edge
%   res (1,1) the residual within the two position
%   interPt (2,1)
%
% [res,peakPt,ind1,ind2] = residual2D_numeric(c1,c2,vec1,vec2,sp1,sp2,method)
%   all the same except that the sp1 and sp2 remain the B-form spline
%   struct pf the tool edge

arguments
    s1
    s2
    cutPt1
    cutPt2
    options.method {mustBeMember(options.method,['DSearchn','BoundingBox'])}
    options.eps {mustBePositive} = 1e-3
end

%% to ensure that the form of the input point remains the 3*1 vector
if isstruct(s1)
    u = 0:options.eps:1;
    sp1 = fnval(s1,u);
    sp2 = fnval(s2,u);
elseif size(cutPt1,1) == 1
    sp1 = s1';
    sp2 = s2';
end

if size(cutPt1,1) == 1
    cutPt1 = cutPt1';
    cutPt2 = cutPt2';
end

% solve the intersection point of the two spline curve
switch options.method
    case 'DSearchn'
        [peakPt,ind1,ind2,epsCross] = pcCrossDN(sp1,sp2);
    case 'BoundingBox'
        [peakPt,ind1,ind2,epsCross] = boundingBox(sp1,sp2);
end

res = norm(cross(cutPt1 - peakPt,cutPt2 - peakPt)) ...
    /norm(cutPt2 - cutPt1); % ?????
varargout{1} = u(ind1);
varargout{2} = u(ind2);
if epsCross >= norm(sp1(:,1) - sp1(:,2))
    % no intersection
    res = nan;
end

end