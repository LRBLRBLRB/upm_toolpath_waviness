function pt = bSplinePt(cpts,k,u,U)
% usage: pt = bSplinePt(cpts,k,u,U)
% solve the point of B-spline curve with parameter u
% Input：
%   order of B-spline: k
%   node parameter: u
%   node vector: U (n+k+1,1)
%   control points: cpt (n,dimention)
% Output:
%   B-spline curve point: pt

[n,dim] = size(cpts); % n: the number of control pts
% dim: dimension of the spline curve
ii = findSpan(n,k,u,U);
N = bBasisFuns(ii,k,u,U);
pt = zeros(1,dim);

% calculate the spline coordinate based on the definition
for jj = 1:k+1
    % indices of non-zero basis functions are "(i-p):(i)"
    pt = pt + N(jj).*cpts(jj+ii-1-k,:);
end

end
