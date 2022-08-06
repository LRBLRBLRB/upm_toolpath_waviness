function ang = vecAng(vec1,vec2,dim)
% to solve the angle between two vectors vec1 and vec2

arguments
    vec1 {mustBeFinite}
    vec2 {mustBeFinite}
    dim {mustBeInteger}
end

ang = acos( ...
    dot(vec1,vec2,dim)./norm(vec1)./norm(vec2) ...
    );
end
