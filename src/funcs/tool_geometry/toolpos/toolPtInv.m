function [u,Q,varargout] = toolPtInv(sp,known,eps,maxIter,options)
% usage 
% to solve the parameter of the B-spline curve with various types of
% inputs
%   sp struct the B-form of the b-spline curve or surface
%   options.Type
%       OpenAngle the open angle of the point is given
%       PolarAngle the polar angle of the point is given
%   the above two are not applicable for spline surfaces
% outputs
%   u (1,1) the parameter of the given point
%   Q(2,1) the coordinate same as Cpts of the given point

arguments
    sp
    known
    eps {mustBeFinite} = 1e-3
    maxIter {mustBeFinite} = 100
    options.Type {mustBeMember(options.Type, ...
        ['OpenAngle','PolarAngle','Cartesian', ...
        'TangentPlane','TangentLine',''])} = 'PolarAngle'
    options.method {mustBeMember(options.method, ...
        ['Newton-Iteration','Traverse'])}
    options.IncludedAng {mustBeFinite}
    options.Radius {mustBeFinite}
end

switch options.Type
    case 'PolarAngle' % DEBUG: the iteration may end to a endless loop!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if strcmp(options.method,'Traverse')
            
        else
            u = 0.5 - (known - pi/2)/options.IncludedAng;
    %         u = fsolve(@(u) funPolar(u,sp),u0);
    %         Q = fnval(sp,u);
            iter = 0;
            while iter <= maxIter
                Q = fnval(sp,u);
                delta = atan2(Q(end),Q(end - 1)) - known;
                if  abs(delta) < eps 
                    break;
                end
                u = u - delta/options.IncludedAng;
                iter = iter + 1;
            end
        end
    case 'OpenAngle'
        u = 0.5 + known/options.IncludedAng;
        iter = 0;
        while iter <= maxIter
            Q = fnval(sp,u);
            delta = atan2(Q(end),Q(end - 1)) - known;
            if  abs(delta) < eps 
                break;
            end
            u = u - delta/options.IncludedAng;
            iter = iter + 1;
        end
    case 'TangentPlane'
        % The input known is the normal vector of a plane that lies outside
        % the tool edge, and the point closest to the plane will be worked
        % out.
        % simply traverse the parameter throughout all the interval
        surfPt = 2*options.Radius*known/norm(known);
        u = 0:eps:1;
        Pt = fnval(sp,u);
        dist = abs(dist2Plane(Pt,surfPt,known));
        [varargout{1},minI] = min(dist);
        u = u(minI);
        Q = fnval(sp,u);
    case 'TangentLine'
        % The input known is the normal vector of a plane that lies outside
        % the tool edge, and the point closest to the plane will be worked
        % out.
        % simply traverse the parameter throughout all the interval
        surfPt = 2*options.Radius*known/norm(known);
        u = 0:eps:1;
        Pt = fnval(sp,u);
        dist = abs(dist2Plane(Pt,surfPt,known));
        [varargout{1},minI] = min(dist);
        u = u(minI);
        Q = fnval(sp,u);
    case 'Cartesian' % solve the parameter of the given point coordinate
        % use the algorithm that stated in the ob
        u = bSplinePtInv(sp,known,eps,maxIter);
end

end


%% 
function theta = funPolar(u,sp)
% sp the B-form of the spline curve
% u the independent variable of the function
P = fnval(sp,u);
theta = atan2(P(2),P(1));
end