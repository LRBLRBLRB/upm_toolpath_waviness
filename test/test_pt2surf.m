addpath(genpath('test'));
EPS = 1e-6;

% surface function
func = @(x,y) y - 3.*x;

% surface dimesion
f_dim = 2;

xps = [0.5,2];

x0 = 1.5;

[x_opt,dist] = cal_distance2surface(func,f_dim,xps,x0,EPS);

% plot the result




%% 计算高维空间中点xps离函数func(x) = 0构成的曲面的距离和最近点坐标
function [x_opt,dist] = cal_distance2surface(func,f_dim,xps,x0,EPS)

% calculate the partial differentiates
funcG = @(x,y,z) func(x,y) - z;
syms x y;
funcG_x = matlabFunction(diff(funcG,x));
funcG_y = matlabFunction(diff(funcG,y));

% non-linear function to be solved
    function eqs = eqs2solve(x)
        pd_values
        eqs = zeros(1,f_dim);
        for i = 1:f_dim
            if pd_values(end) == 0
                pd_values(end) = EPS;
            end
            if pd_values(i) == 0
                pd_values(i) = EPS;
            end
            eqs(i) = (xps(i) - x(i)) / (pd_values(i) / pd_values(end)) - xps(end) + x(end);
        end
    end

root = fsolve(eqs2solve,x0)



end