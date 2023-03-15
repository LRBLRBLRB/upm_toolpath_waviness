function [dist,ptProj] = dist2curve(pt,curveFunc,curveFunc_y,options)
%DIST2SURF calculate the distance between the given point to the surface.
% The distance is named as the oriented distance from the point pt to the
% surface surfFunc. If the pt is on the same side of the surface based on
% the normal vector, then the distance will be positive.
%
% pt         (3,:) the points set, the distance of which would be calculated 
% curveFunc  sym   the surface function, i.e., surfFunc = f(x,y,z) = 0
% curveFunc_y 
% options    

% 问题：这里针对并行运算（多个初值同时求解）的优化是否正确？

arguments
    pt (2,:)
    curveFunc
    curveFunc_y = []
    options.CalculateType {mustBeMember(options.CalculateType, ...
        {'Taylor-Expand','Lagrange-Multiplier'})} = 'Lagrange-Multiplier'
    options.DisplayType {mustBeMember(options.DisplayType,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
end

optimOpt = optimoptions('fsolve','Display',options.DisplayType, ...
    'FunctionTolerance',1e-6,'MaxIterations',500,'UseParallel',false);

switch options.CalculateType
    case 'Taylor-Expand'
        if isempty(curveFunc_y)
            syms y;
            curveFunc_y = matlabFunction(diff(curveFunc,y),'Vars',y);
        end
        % ptProj
        % dist
    case'Lagrange-Multiplier'
        if isempty(curveFunc_y)
            syms y;
            curveFunc_y = matlabFunction(diff(curveFunc,y),'Vars',y);
        end
        ptProj = fsolve(@(Q) eqs2solve(Q,pt,curveFunc,curveFunc_y),pt(1,:),optimOpt);
        ptProj(2,:) = curveFunc(ptProj(1,:));
        dist = vecnorm(pt - ptProj,2,1);

        % calculate the direction
        ifSame = pt(end,:) - ptProj(end,:) < 0; % QP & norm are opposite
        dist(ifSame) = -dist(ifSame);
end

end

function F = eqs2solve(Q,P,func,func_x)
    F = func_x(Q(1,:)).*(P(1,:) - Q(1,:)) + (P(2,:) - func(Q(1,:)));
end