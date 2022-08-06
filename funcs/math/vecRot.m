%% calculate the rotation matrix from vec1 to vec2
function R = vecRot(vec1,vec2)
u = cross(vec1,vec2);
if u == 0
    R = eye(3);
else
    u = u./norm(u);
    if size(u,1) ~= 1
        u = u';
    end
    theta = vecAng(vec1,vec2,1);
    R = transpose(rotationVectorToMatrix(theta*u));
end
end
