function [dist,ptProj] = dist2surf(pt,surfFunc,surfFunc_x,surfFunc_y,options)
%DIST2SURF calculate the distance between the given point to the surface
% pt         (3,:) the points set, the distance of which would be calculated 
% surfFunc   sym   the surface function, i.e., surfFunc = f(x,y,z) = 0
% surfFunc_x 
% surfFunc_y
% options    

% 问题：这里针对并行运算（多个初值同时求解）的优化是否正确？

arguments
    pt (3,:)
    surfFunc
    surfFunc_x = []
    surfFunc_y = []
    options.CalculateType {mustBeMember(options.CalculateType, ...
        {'Taylor-Expand','Lagrange-Multiplier'})} = 'Lagrange-Multiplier'
end

optimOpt = optimoptions('fsolve','Display','final-detailed', ...
    'FunctionTolerance',1e-6,'MaxIterations',500);

switch options.CalculateType
    case 'Taylor-Expand'
        
    case'Lagrange-Multiplier'
        if isempty(surfFunc_x)
            syms x y;
            surfFunc_x = matlabFunction(diff(surfFunc,x),'Vars',{x,y});
            surfFunc_y = matlabFunction(diff(surfFunc,y),'Vars',{x,y});
        end
        ptProj = fsolve(@(Q) eqs2solve(Q,pt,surfFunc,surfFunc_x,surfFunc_y),pt(1:2,:),optimOpt);
        ptProj(3,:) = surfFunc(ptProj(1,:),ptProj(2,:));
        dist = vecnorm(pt - ptProj,2,1);
end

end

function F = eqs2solve(Q,P,func,func_x,func_y)
    F(1,:) = (P(1,:) - Q(1,:)) + func_x(Q(1,:),Q(2,:)) .* (P(3,:) - func(Q(1,:),Q(2,:)));
    F(2,:) = (P(2,:) - Q(2,:)) + func_y(Q(1,:),Q(2,:)) .* (P(3,:) - func(Q(1,:),Q(2,:)));
end