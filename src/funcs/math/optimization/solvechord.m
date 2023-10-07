function [x,iterNum] = solvechord(f,x0,x,eps,maxIter)
%SOLVECHORD use the chord truncation method to solve the f(x) = 0
% phi3 = @(x1,x2) x2-(x2.^3+4*x2-1).*(x2-x1)./(x2.^3+4*x2-1-(x1.^3+4*x1-1));
phi = @(x0,x) x - (f(x).*(x - x0))./(f(x) - f(x0));
iterNum = zeros(1,length(x));
while max(iterNum) <= maxIter
    logi = (abs(f(x)) <= eps);
    if all(logi), break; end
    tmp = ~logi.*phi(x0,x);
    tmp(isnan(tmp)) = 0;
    x0 = x + logi.*x0;
    x = tmp + logi.*x;
    iterNum = iterNum + 1*~logi;
end
end