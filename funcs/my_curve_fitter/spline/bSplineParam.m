function [u,Q] = bSplineParam(sp,known,eps,options)
% usage 
% to solve the parameter of the B-spline curve with various types of
% inputs.

arguments
    sp
    known
    eps {mustBeFinite} = 1e-3
    options.Type {mustBeMember(options.Type, ...
        ['OpenAngle','PolarAngle',''])} = 'PolarAngle'
    options.IncludedAng {mustBeFinite}
end

switch options.Type
    case 'PolarAngle'
        u = 0.5 - (known - pi/2)/options.IncludedAng;
        while true
            Q = fnval(sp,u);
            delta = atan2(Q(2),Q(1)) - known;
            if  abs(delta) < eps 
                break;
            end
            u = u - delta/options.IncludedAng;
        end
    case 'OpenAngle'
        u = 
end


end