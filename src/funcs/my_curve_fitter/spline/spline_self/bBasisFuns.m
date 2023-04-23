function N = bBasisFuns(i,k,u,U)
% usage: N = bBasisFuns(i,u,p,U)
% 求出在u处所有的非零p次B样条基函数
% Inputs：
%   参数所处节点区间下标：i
%   order of B-spline basis function：k
%   node paramters：u
%   node vector: U
% Outputs：
%   所有非零基函数的序列：N（在基函数序列中的下标为i-p+1,...,i+1，下标从1开始）

N = zeros(k+1,1);
left = zeros(k,1);
right = zeros(k,1);
N(1) = 1; % p=0: Ni=1
for j = 1:k
    left(j) = u - U(i+1-j);
    right(j) = U(i+j) - u;
    saved = 0; % to store the replicant value
    for r = 1:j
        temp = N(r)/(left(j-r+1)+right(r));
        N(r) = saved + right(r)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1) = saved;
end
end
