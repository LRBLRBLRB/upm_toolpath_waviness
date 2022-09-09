function interPt = boundingBox(s1,s2)
%% 计算两端B样条曲线s1、s2的交点 interPt
% 最简单的算法：包围盒法
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

interPt = boundingBox(s1(),s2());
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