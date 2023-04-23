function mid = findSpan(n,p,u,U)
% 确定参数u在节点矢量U中所在的区间的下标
% 输入：
%   控制点数：n
%   基函数次数：p
%   参数：u
%   节点矢量：U
% 输出：
%   u所在的节点区间的左侧节点下标：mid
% 注意——MATLAB的数组下标从1开始，因此涉及数组下标的都要加1

if u == U(n+1)
    mid = n;
    return;
end
low = p+1;
high = n+1;
mid = floor((low + high)/2);
while u < U(mid) || u >= U(mid + 1)
    if u < U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low + high)/2);
end
end

% 方法二
% function lo = findSpan(n,p,u,U)
% if u == U(n+1)
%     lo = n;
%     return;
% end
% lo = p + 1;
% hi = n + 1;
% while lo < hi
%     mid = floor((lo + hi)/2);
%     if u < U(mid)
%         hi = mi;
%     else
%         lo = mi + 1;
%     end
% end
% lo = lo - 1;
% end
