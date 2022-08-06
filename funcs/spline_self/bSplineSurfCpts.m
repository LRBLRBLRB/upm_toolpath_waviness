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

%% global interpolation: solve node params corresponding to Q
if strcmp(options.InterpMethod,'uniform') % uniform parameterization

else
    if strcmp(options.InterpMethod,'concentric') % concentric parameterization
        % sqrt of dist of adjacent Q, size(n-1,1)
        uDist = sum((Q(2:end,:,:)-Q(1:end-1,:,:)).^2,3).^(1/4); % (m-1,n)
        vDist = sum((Q(:,2:end,:)-Q(:,1:end-1,:)).^2,3).^(1/4); % (m,n-1)
    else % chord parameterization
        % dist of adjacent Q, size(n-1,1)
        uDist = sum((Q(2:end,:,:)-Q(1:end-1,:,:)).^2,3).^(1/2); % (m-1,n)
        vDist = sum((Q(:,2:end,:)-Q(:,1:end-1,:)).^2,3).^(1/2); % (m,n-1)
    end
    uSum = sum(uDist,1); % (1,n)
    vSum = sum(vDist,2); % (m,1)
    % 1st param of Q
    uQp = zeros(m,n);
    for ii = 2:m
        uQp(ii,:) = uQp(ii-1,:) + uDist(ii-1,:);
    end
    uQp = uQp./uSum;
    uQp(isnan(uQp)) = 0; 
    uQ = mean(uQp,2);
    uQ(end) = 1; % elliminate calcuation error
    % 2nd param of Q
    vQq = zeros(m,n);
    for jj = 2:n
        vQq(:,jj) = vQq(:,jj-1) + vDist(:,jj-1);
    end
    vQq = vQq./vSum;
    vQq(isnan(vQq)) = 0;
    vQ = transpose(mean(vQq,1));
    vQ(end) = 1;
end

%% node vector U & V generation
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
