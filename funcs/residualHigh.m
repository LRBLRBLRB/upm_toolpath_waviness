% 计算两个刀位间的残高
% 输入两个刀位的位置x、方向vec、刀尖轮廓s
% 假设残高就是两个刀位的刀尖轮廓交点到两个刀位的刀尖点连线的距离
function [res,interPt] = residualHigh(x1,vec1,s1,x2,vec2,s2)
% 求两个刀位的轮廓交点
% [interPt,~] = bsplineCross(s1,s2);
[interPt,~] = bsplineCrossDN(s1,s2);

% 求刀尖点：刀尖方向上最远点
[~,cutPtIndex1] = min(dot((s1-x1),meshgrid(vec1,s1(:,1)))/norm(vec1));
[~,cutPtIndex2] = min(dot((s2-x2),meshgrid(vec2,s2(:,1)))/norm(vec2));
cutPt1 = s1(cutPtIndex1,:);
cutPt2 = s2(cutPtIndex2,:);
% 求刀尖点：两个刀尖的公切线


if size(interPt,2) == 2
    interPt = [interPt,0];
    cutPt1 = [cutPt1,0];
    cutPt2 = [cutPt2,0];
end
res = 2*abs(cross(interPt-cutPt1,interPt-cutPt2))/norm(cutPt1-cutPt2);
end

%% 计算两端B样条曲线s1、s2的交点 interPt
% 最简单的算法：包围盒法
function interPt = bsplineCross(s1,s2)
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

interPt = bsplineCross(s1(),s2());
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

%% MATLAB自带的求距离函数
function [interPt,eps] = bsplineCrossDN(s1,s2)
[index,dist]=dsearchn(s1,s2);
[eps,I2]=min(dist);
I1 = index(I2);
interPt = 0.5*(s1(I1,:)+s2(I2,:));
end