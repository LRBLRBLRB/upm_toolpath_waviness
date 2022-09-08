function R = vecRot(vec1,vec2)
% calculate the rotation matrix from vec1 to vec2
arguments
    vec1 {mustBeFinite}
    vec2 {mustBeFinite}
end
vec1 = reshape(vec1,1,3); % to ensure that vec1 is a row vector
vec2 = reshape(vec2,1,3);
u = cross(vec1,vec2);
if u == 0
    R = eye(3);
else
    u = u./norm(u);
    theta = vecAng(vec1,vec2,2);
    q = [cos(theta/2),sin(theta/2).*u];
    R = [2*q(1).^2 + 2*q(2).^2 - 1,  2*(q(2).*q(3) - q(1).*q(4)), 2*(q(2).*q(4) + q(1).*q(3));
        2*(q(2).*q(3) + q(1).*q(4)), 2*q(1).^2 - 1 + 2*q(3).^2,   2*(q(3).*q(4) - q(1).*q(2));
        2*(q(2).*q(4) - q(1).*q(3)), 2*(q(3).*q(4) + q(1).*q(2)), 2*q(1).^2 - 1 + 2*q(4).^2];
end
end