function [res,peakPt,varargout] = residual2D_numeric(c1,c2,vec1,vec2,s1,s2,options)
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
    c1 {mustBeFinite}
    c2 {mustBeFinite}
    vec1 {mustBeFinite}
    vec2 {mustBeFinite}
    s1
    s2
    options.method {mustBeMember(options.method,['DSearchn','BoundingBox'])}
    options.eps {mustBePositive} = 1e-3
end

if size(c1,1) ~= 1
    c1 = c1';
    c2 = c2';
    vec1 = vec1';
    vec2 = vec2';
end

if isstruct(s1)
    u = 0:options.eps:1;
    sp1 = s1;
    s1 = (fnval(sp1,u))';
    sp2 = s2;
    s2 = (fnval(sp2,u))';
elseif size(s1,1) ~= 1
    s1 = s1';
    s2 = s2';
end

% solve the intersection point of the two spline curve
switch options.method
    case 'DSearchn'
        [peakPt,ind1,ind2,epsCross] = pcCrossDN(s1,s2);
    case 'BoundingBox'
        [peakPt,ind1,ind2,epsCross] = curveCross(s1,s2);
end

res = norm(cross(cutPt1 - peakPt,cutPt2 - peakPt)) ...
    /norm(cutPt2 - cutPt1); % ?????
varargout{1} = ind1;
varargout{2} = ind2;
if epsCross >= norm(s1(:,1) - s1(:,2))
    % no intersection
    res = nan;
end

end


%% 计算两端B样条曲线s1、s2的交点 interPt
% 最简单的算法：包围盒法
function interPt = curveCross(s1,s2)
if size(s1,1) <= 2 || size(s2,1) <= 2
    interPt = s1(1,:); % 这是比较粗略的处理，就把s1的最初一个输出，这里的递归终止条件还要深入讨论以提高精度
    return;
end
% min1 = min(s1,[],1);
% max1 = max(s1,[],1);
% min2 = min(s2,[],1);
% max2 = max(s2,[],1);
% isinter = recInter(min1,max1,min2,max2);

% 如果包围盒有共同区域就一分为二递归下去
n1 = length(s1);
n2 = length(s2);

interPt = curveCross(s1(),s2());
end

function isinter = recInter(min1,max1,min2,max2)
% 判断包围盒是否相交：如果两个矩形相交，即两矩形中心点间x、y方向距离分别小于或等于x、y边长和的一半。
mid = abs((min1 + max1) - (min2 + max2));
edge = max1 - min1 + max2 - min2;
if mid-edge>0 %包围盒不相交，交点为空
    isinter = 0;
else
    isinter = 1;
end
end