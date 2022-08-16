function ang = vecAng(vec1,vec2,dim)
% to solve the angle between two vectors vec1 and vec2
% dim refers to the dimension of vec1 that is regarded as the vector, e.g.,
% if vec1 = [0,0,1], vec2 = [0,1,0], then dim = 2.
arguments
    vec1 {mustBeFinite}
    vec2 {mustBeFinite}
    dim {mustBeInteger}
end


ang = acos( ...
    dot(vec1,vec2,dim)./vecnorm(vec1,2,dim)./vecnorm(vec2,2,dim) ...
    );

% ang2 = atan2( ...
%     vecnorm(cross(vec1,vec2,dim),2,dim), ...
%     dot(vec1,vec2,dim) ...
%     );

end
