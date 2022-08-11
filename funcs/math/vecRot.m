%% calculate the rotation matrix from vec1 to vec2
function R = vecRot(vec1,vec2)
u = cross(vec1,vec2);
if u == 0
    R = eye(3);
else
    u = u./norm(u);
    if size(u,2) ~= 3
        u = u';
    end
    theta = vecAng(vec1,vec2,1);
    q = [cos(theta/2),sin(theta/2).*u];
    R = [2*q(1).^2 + 2*q(2).^2 - 1,  2*(q(2).*q(3) + q(1).*q(4)), 2*(q(2).*q(4) - q(1).*q(3));
        2*(q(2).*q(3) - q(1).*q(4)), 2*q(1).^2 - 1 + 2*q(3).^2,   2*(q(3).*q(4) + q(1).*q(2));
        2*(q(2).*q(4) + q(1).*q(3)), 2*(q(3).*q(4) - q(1).*q(2)), 2*q(1).^2 - 1 + 2*q(4).^2];
    % R = transpose(rotationVectorToMatrix(theta*u));
end
end