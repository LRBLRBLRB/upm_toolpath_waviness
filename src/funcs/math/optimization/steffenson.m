function [x,iterNum] = steffenson(f,phi,x,eps,maxIter)
iterNum = zeros(1,length(x));
while max(iterNum)<=maxIter
    logi = (abs(f(x)) <= eps);
    if all(logi), break; end
    tmp = ~logi.* ...
        (x.*phi(phi(x))-phi(x).^2)./(x-2*phi(x)+phi(phi(x)));
    tmp(isnan(tmp)) = 0;
    x = tmp+logi.*x;
    iterNum = iterNum+1*~logi;
end
end