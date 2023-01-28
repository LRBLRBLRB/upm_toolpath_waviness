function [dist,eps] = dist2surf(pt,surfFunc,options)
%DIST2SURF calculate the distance between the given point to the surface
% pt        (3,:) the points set, the distance of which would be calculated 
% surfFunc  sym   the surface function, i.e., surfFunc = f(x,y,z) = 0

arguments
    pt (3,:)
    surfFunc
    options.CalculateType {mustBeMember(options.CalculateType, ...
        {'Taylor-Expand','Lagrange-Multiplier'})}
end

switch options.CalculateType
    case 'Taylor-Expand'
        
    case'Lagrange-Multiplier'
end


end

