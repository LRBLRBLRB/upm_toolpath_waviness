function L = arclengthoncurve(curvFunc,t1,t2,options)
%ARCLENGTHONCURVE Calculate the arc length of the curve curvFunc between t1
% and t2 
%
% L = arclengthoncurve(curvFunc,t1,t2,options)
% Inputs:
%   curvFunc    function_handle the equation of the curve
%   t1          the beginning point
%   t2          the end point
%   option
%       funcType the equation type of the curve
% Output
%   L           the arc length

arguments
    curvFunc function_handle
    t1
    t2
    options.funcType {mustBeMember(options.funcType,{'parametric','polar','cartesian'})} = 'polar'
end

syms t;

switch options.funcType
    case 'parametric' % [x,y,z,...] = f(t)
        dr = diff(curvFunc,t);
        intFunc = sqrt(sum(dr.^2));
        L = vpa(int(intFunc,t,t1,t2));
        % intFunc = matlabFunction(sqrt(sum(dr.^2)),'Vars',t);
        % L = integral(intFunc,t1,t2);
    case 'polar' % r = r(t)
        dr = diff(curvFunc,t);
        intFunc = sqrt(sym(curvFunc).^2 + dr.^2);
        L = vpa(int(intFunc,t,t1,t2));
        % intFunc = matlabFunction(sqrt(sym(curvFunc).^2 + dr.^2),'Vars',t);
        % L = integral(intFunc,t1,t2);
    case 'cartesian'
        L = 0;
end

end