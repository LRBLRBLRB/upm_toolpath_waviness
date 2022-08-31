function x0 = splineTrav(eq,sp,x0,eps,maxIter)
% traverse the spline based on to find out the particular parameter 'u' to 
% saatisfy the given equation 'eqs = 0'
%
% Usage:
%   u = splineTrav(

iterNum = zeros(1,length(x0));
% case Newton
f = @(x) x - eps./diff(eq);
f = @(x) x - fnval(sp,x)./fnval(fnder(sp,1),x);

while max(iterNum) < maxIter
    logi = (abs(eq(x0)) <= eps);
    if all(logi), return; end
    tmp = ~logi.*f(x0);
    tmp(isnan(tmp)) = 0;
    x0 = tmp+logi.*x0;
    iterNum = iterNum+1*~logi;
end

end
