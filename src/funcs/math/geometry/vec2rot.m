function R = vec2rot(v)
% to solve the rotation matrix R corresponding to the rotation vector v
% inputs:
%   v (1,3) or (3,1) rotation vector, the length of which represents the
%   rotation angle, and the orientation of which stands the rotation axis
% outputs:
%   R (3,3) rotation matrix

% use the Rodrigues formula
theta = norm(v);
v = v ./ theta;
V = zeros(3,3);
V(2) = v(3);
V(3) = -v(2);
V(4) = -v(3);
V(6) = v(1);
V(7) = v(2);
V(8) = -v(1);
R = eye(3) + sin(theta)*V + (1 - cos(theta))*V*V;
% RR = transpose(rotationVectorToMatrix(v*theta));
% from Computer Vision Toolbox
end
