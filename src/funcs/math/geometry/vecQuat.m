function q = vecQuat(vec1,vec2)
% calculate the unit quarternion from vec1 to vec2
arguments
    vec1 {mustBeFinite}
    vec2 {mustBeFinite}
end
vec1 = reshape(vec1,1,3); % to ensure that vec1 is a row vector
vec2 = reshape(vec2,1,3);
u = cross(vec1,vec2);
if u == 0
    q = [1,0,0,0];
else
    u = u./norm(u);
    theta = vecAng(vec1,vec2,2);
    q = [cos(theta/2),sin(theta/2).*u];
end
end