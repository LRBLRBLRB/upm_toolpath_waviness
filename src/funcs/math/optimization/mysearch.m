function [x,iter] = mysearch(f,x0,h,x1,eps)
%MYSEARCH calculate the root of a monotonically increasing function
%   此处显示详细说明

%% initial interval
% fprintf('initial linear search begins: \n');
tmp = f(x0);
if tmp == 0
    x = x0;
    iter = 0;
%     fprintf('f(%f) = 0, and the iteratoin ends.\n',x0);
    return;
elseif tmp > 0
    b = x0;
    a = x0;
    ii = 0;
%     fprintf('No.\tx\tf(x)\n');
    tmp = f(a);
%     fprintf('%d\t%f\t%f\n',ii,a,tmp);
    while tmp > 0
        ii = ii + 1;
        b = a;
        a = a - h;
        tmp = f(a);
%         fprintf('%d\t%f\t%f\n',ii,a,tmp);
        if a > x1(2) 
            b = a + h;
            a = x1(2);
%             fprintf('Linear search failed. all the f(x) > 0.');
            break;
        elseif a < x1(1)
            b = a + h;
            a = x1(1);
%             fprintf('Linear search failed. all the f(x) > 0.');
            break;
        end
    end
else
    a = x0;
    b = x0;
    ii = 0;
%     fprintf('No.\tx\tf(x)\n');
    tmp = f(b);
%     fprintf('%d\t%f\t%f\n',ii,b,tmp);
    while tmp < 0
        ii = ii + 1;
        a = b;
        b = b + h;
%         tmp0 = tmp;
        tmp = f(b);
%         fprintf('%d\t%f\t%f\n',ii,b,tmp);
        if b > x1(2)
            a = b - h;
            b = x1(2);
%             fprintf('Linear search failed. all the f(x) < 0.');
            break;
        elseif b < x1(1)
            a = b - h;
            b = x1(1);
%             fprintf('Linear search failed. all the f(x) < 0.');
            break;
        end
    end
end
if a > b
    tmp = a;
    a = b;
    b = tmp;
end
% fprintf('bisection interval found by 1-dim search: [%f,%f]\n',a,b);

%% bisection
[iter,x] = isBisection(f,a,b,eps);

end