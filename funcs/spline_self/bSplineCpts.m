function [cpts,U] = bSplineCpts(Q,k,interpMethod)
% [cpts,U] = bSplineCpts(Q,k,option)
% To solve the k-order B-spline control pts 'cpts' with data points Q
% Input:
%   matrix of data points: Q (n,3)
%   order of B spline: k
%   interpolation method: option.InterpMethod
% Output:
%   matrix of control points：cpts(n,3)
%   node vector: U (n+k+1,1)
% E.g. Q = [0 0 0;3 4 0;-1 4 0;-4 0 0;-4 -3 0;0 -3 0]; k = 3;

arguments
    Q {mustBeFinite}
    k (1,1) {mustBeInteger} = 3
    interpMethod {mustBeMember(interpMethod, ...
        ['uniform','concentric','chord',''])} = 'chord'
end

[n,~] = size(Q); % n, number of control pts; dim, dimension of the curve
% m = n+k+1; % the number of the nodes
% test if the cpts are enough to generate the quasi-uniform B-spline
if n-k-1<0
    error("control points are not enough, or the required order is too large");
end

%% global interpolation: solve node params corresponding to Q
if strcmp(interpMethod,'uniform') % uniform parameterization

else
    if strcmp(interpMethod,'concentric') % concentric parameterization
        l = sum((Q(1:end-1,:)-Q(2:end,:)).^2,2).^(1/4); % sqrt of dist of adjacent Q, size(n-1,1)
    else % chord parameterization
        l = sum((Q(1:end-1,:)-Q(2:end,:)).^2,2).^(1/2); % dist of adjacent Q, size(n-1,1)
    end
    L = sum(l);
    uQ = zeros(n,1);
    for ii = 2:n-1
        uQ(ii) = uQ(ii-1) + l(ii-1);
    end
    uQ = uQ/L; % control points are on the spline curve
    uQ(end) = 1;
end

%% node vector U generation
U = nodeVector(k,n,'nodeMethod','Interpolation','uQ',uQ);

%% equations to solve the cpts
A = zeros(n,n);
for i = 1:n
   span = findSpan(n,k,uQ(i),U); 
   N = bBasisFuns(span,k,uQ(i),U);
   for j = 0:k
      A(i,span-k+j) = N(j+1); 
   end
end
cpts = A\Q;

end
