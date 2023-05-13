function sp=bspline_ls_clampedder(q,n,u,U,p,der)
% 首末端点固定且带一阶导矢约束的B样条最小二乘拟合
% cpt：控制点
% q：拟合的离散点
% n: 控制点个数
% u：离散点对应的参数
% U：节点矢量
% p：B样条曲线次数
% der：两端一阶导矢约束

M=size(q,1);
A=zeros(M-2,n);%非约束
B=spcol(U,p+1,[u(1) u(1) u(end) u(end)]);%端点位置和一阶导矢约束
C=q(2:end-1,:);
D=[q(1,:);der(1,:);q(end,:);der(end,:)];

for i=2:M-1
   span=FindSpan(n,p,u(i),U);
   N=BBasicFuns(span,u(i),p,U);
   for j=0:p
      A(i-1,span-p+j)=N(j+1); 
   end
end

left=[A.'*A,B.';B,zeros(4,4)];
right=[A.'*C;D];
solution=left\right;
cpt=solution(1:n,:);
sp=spmak(U,cpt.');

end