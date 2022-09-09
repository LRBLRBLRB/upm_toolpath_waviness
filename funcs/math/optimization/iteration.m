function varargout = iteration(f,phi,x,eps,maxIter)
% Iteration process to solve the non-linear equation
%
%
%   f the equation that satisfies "f(x) == 0"
%   phi the iterative function, equals to x-f/f'
%   x the interative initial value

iterNum = zeros(1,length(x));
while max(iterNum) < maxIter
    % store whether each value satisfies the equation f==0
    logi = (abs(f(x)) <= eps);
    % if all remain 0, then the iteration ends.
    if all(logi), break; end
    % else, 
    tmp = ~logi.*phi(x);
    tmp(isnan(tmp)) = 0;
    x = tmp + logi.*x;
    iterNum = iterNum + 1*~logi;
end
varargout{1} = x;
varargout{2} = iterNum;
end