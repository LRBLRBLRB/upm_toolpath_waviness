function [cpts,U,V] = bSplineSurfCpts(Q,k,l,options)
% [cpts,U] = bSplineCpts(Q,k,option)
% To solve the k-order B-spline control pts 'cpts' with data points Q
% Input:
%   matrix of data points: Q (m,n,3)
%   order of B spline: k, l
%   interpolation method: option.InterpMethod
% Output:
%   matrix of control points：cpts(n,3)
%   node vector: U (m+k+1,1), V(n+l+1,1)
% E.g. Q = [0 0 0;3 4 0;-1 4 0;-4 0 0;-4 -3 0;0 -3 0]; k = 3;

arguments
    Q {mustBeFinite}
    k (1,1) {mustBeInteger} = 3
    l (1,1) {mustBeInteger} = 3
    options.InterpMethod {mustBeMember(options.InterpMethod, ...
        ['uniform','concentric','chord',''])} = 'chord'
end

[m,n,~] = size(Q);
% test if the cpts are enough to generate the quasi-uniform B-spline
if m-k-1<0 || n-l-1<0
    error("control points are not enough, or the required order is too large");
end

%% global interpolation
% solve node params corresponding to Q
[uQ,vQ] = interpParamSurf(Q,option.InterpMethod);

% node vector U & V generation
U = nodeVector(k,m,'uQ',uQ,'nodeMethod','Interpolation');
V = nodeVector(l,n,'uQ',vQ,'nodeMethod','Interpolation');

%% equation to solve the cpts
cpts = zeros(m,n,3);
R = zeros(m,n,3);
for q = 1:n
    A = zeros(m,m);
    for i = 1:m
       span = findSpan(m,k,uQ(i),U); 
       N = bBasisFuns(span,k,uQ(i),U);
       for j = 0:k
          A(i,span-k+j) = N(j+1); 
       end
    end
    tmpR = A\reshape(Q(:,q,:),[m 3]);
    R(:,q,:) = reshape(tmpR,[m 1 3]);
end

for p = 1:m
    A = zeros(n,n);
    for i = 1:n
       span = findSpan(n,l,vQ(i),V); 
       N = bBasisFuns(span,l,vQ(i),V);
       for j = 0:l
          A(i,span-l+j) = N(j+1); 
       end
    end
    tmpCpt = A\reshape(R(p,:,:),[n 3]);
    cpts(p,:,:) = reshape(tmpCpt,[1,n,3]);
end

end
