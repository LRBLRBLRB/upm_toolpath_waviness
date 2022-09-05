function [res,interPt,varargout] = residual2D(c1,c2,vec1,vec2,varargin)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,interPt] = residual2D(c1,c2,vec1,vec2,s1,s2)
%   vec1 (2,1) 1st normal vector of the tool edge
%   s1 (2,:) 1st scatters of the tool edge
%   c2,vec2,s2 are parameters of the other tool edge
%   res (1,1) the residual within the two position
%   interPt (2,1)
%
% [res,interPt] = residual2D(c1,c2,vec1,vec2,sp)
%   all the same except that the sp1 and sp2 remain the B-form spline
%   struct pf the tool edge

if nargin == 5 % or isstruct(sp1)
    % 继续写这个！！！！
else
    sp1 = varargin{1};
    sp2 = varargin{2};
    % 求两个刀位的轮廓交点
    % [interPt,~] = bsplineCross(sp1,sp2);
    [interPt,~] = pcCrossDN(sp1,sp2);
    
    % 求刀尖点：刀尖方向上最远点
    [~,cutPtIndex1] = max(abs(dot((sp1-c1),ndgrid(vec1,sp1(1,:)))/norm(vec1)));
    [~,cutPtIndex2] = max(abs(dot((sp2-c2),ndgrid(vec2,sp2(1,:)))/norm(vec2)));
    cutPt1 = sp1(:,cutPtIndex1);
    cutPt2 = sp2(:,cutPtIndex2);
    % 求刀尖点：两个刀尖的公切线
    
    if size(interPt,1) == 2
        interPt = [0;interPt];
        cutPt1 = [0;cutPt1];
        cutPt2 = [0;cutPt2];
    end
    res = 2*abs(cross(interPt-cutPt1,interPt-cutPt2))/norm(cutPt1-cutPt2);
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