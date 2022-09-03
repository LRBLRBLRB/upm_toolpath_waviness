function x0 = spSolveIter(sp,x0,eps,maxIter,f)
% traverse the spline based on to find out the particular parameter 'u' to 
% saatisfy the given equation 'eqs = 0'
%
% Usage:
%   u = splineTrav(

iterNum = zeros(1,length(x0));

% case Newton
if nargin == 4
    f = @(x) x - fnval(sp,x)./fnval(fnder(sp,1),x);
end

while max(iterNum) < maxIter
    logi = (abs(fnval(sp,x0)) <= eps);
    if all(logi), return; end
    tmp = ~logi.*f(x0);
    tmp(isnan(tmp)) = 0;
    x0 = tmp + logi.*x0;
    iterNum = iterNum + 1*(~logi);
end

end
