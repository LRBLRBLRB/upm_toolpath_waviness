function [iter,x] = isBisection(f,a,b,eps)
% fprintf("Bisection process begins: \n");
place = -log10(eps);
iter = 0; % 迭代次数
intDist = (b-a)/2; % 区间长度
fa = f(a);
fb = f(b);
if(fa*fb > 0)
    fprintf("Wrong input: there is no root in the interval [a,b]!\n");
    iter = -1; x = nan;
    return;
elseif(fa==0)
    x = a;
%     fprintf("The value of x=a remains 0.\nBisection ends!\n");
    return;
elseif(fb==0)
    x = b;
%     fprintf("The value of x=b remains 0.\nBisection ends!\n");
    return;
else
    c = (a+b)/2;
%     fprintf("k\t%*s\t%*s\n", place+4, 'xk', place+4, "eps"); %xk\t\t\t\t
%     fprintf("0\t%*.*f\t%.3e\n", place+4, place, (a+b)/2, intDist);
    while intDist > eps
        iter = iter+1;
        fc = f(c);
        if(fa*fc<0)
            b = c;
        else
            a = c;
        end
        c = (a+b)/2;
        intDist = (b-a)/2;
%         fprintf("%d\t%*.*f\t%.3e\n", iter, place+4, place, c, intDist);
    end
%     fprintf("Bisection ends!\n");
    x = (a+b)/2;
end
end
