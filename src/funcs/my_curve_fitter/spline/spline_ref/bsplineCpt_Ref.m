% 由给定插值点Q反求B样条曲线控制点P
% Q：n*3
% k：B样条次数，且为k次非均匀B样条
% P：n*3
function P = bsplineCpt_Ref(Q,k)
% Q = [0 0 0;3 4 0;-1 4 0;-4 0 0;-4 -3 0;0 -3 0];
% k = 3; % order of the B spline
n = size(Q,1); % 插值点个数（即控制点个数）
m = n+k+1; % 节点个数
if n-k-1<0
    printf("插值点个数太少或样条次数太高。");
    return;
end

%% 求节点矢量
l = sum((Q(1:end-1,:)-Q(2:end,:)).^2,2).^(1/2); % 相邻插值点的距离，(n-1,1)
L = sum(l);
% Riesenfeld
UTemp = zeros(n-k-1,1);
if mod(k,2) == 0
    tmp = sum(l(1:k/2));
    for i = 1:n-k-1
        tmp = tmp + l((k+2*i)/2);
        UTemp(i) = tmp/L;
        tmp = tmp + l((k+2*i+2)/2);
    end
else
    tmp = sum(l(1:(k-1)/2));
    for i = 1:n-k-1
        tmp = tmp + l((k+2*i-1)/2);
        UTemp(i) = tmp/L;
    end
end

% Hartley-Judd
% deno = 
% for i = 1:n-k-1
%     
% end

% 弦长参数化的非均匀B样条节点矢量
U = [zeros(k+1,1); ...
    UTemp; ...
    ones(k+1,1)];

%% 列反求控制点的方程
A = zeros(n,n);
dist = zeros(n,1);
for ii = 2:n
    dist(ii) = dist(ii-1) + l(ii-1);
end
dist = dist/L;
for i=1:n
   span = FindSpan(n,k,dist(i),U); 
   N = BasicFuns(span,dist(i),k,U);
   for j=0:k
      A(i,span-k+j)=N(j+1); 
   end
end
P = A\Q;

end

%% 相关函数
% 获取下标位置
function mid = FindSpan(n,k,u,U)
if u == U(n+1)
    mid = n;
    return;
end
low = k;
high = n+1;
mid = floor((low+high)/2);
while u<U(mid) || u>=U(mid+1)
    if u<U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
end

% 计算B样条基函数：类似于树
function N = BasicFuns(span,u,k,U)
N = zeros(k+1,1);
left = zeros(k,1);
right = zeros(k,1);
N(1) = 1;
for j = 1:k
    left(j) = u-U(span+1-j);
    right(j) = U(span+j)-u;
    saved = 0;
    for r = 0:j-1
        temp = N(r+1)/(right(r+1)+left(j-r));
        N(r+1) = saved + right(r+1)*temp;
        saved = left(j-r)*temp;
    end
    N(j+1) = saved;
end
end
% end