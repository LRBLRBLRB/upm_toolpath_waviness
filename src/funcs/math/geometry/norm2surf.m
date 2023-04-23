function normVec = norm2surf(pts,surfType)
% to calculate the normal vector of the surface


arguments
    pts
    surfType {mustBeMember(surfType, ...
        ['plane','xyz_surface','tri_patch'])} = 'tri_patch'
end

switch surfType
    case 'plane'
    case 'xyz_surface'
    case 'tri_patch'
        vec2 = pts(2,:) - pts(1,:);
        vec3 = pts(3,:) - pts(1,:);
        normVec = cross(vec2,vec3);
end
normVec = normVec./norm(normVec);
end