function R = quat2rot(q)
% to solve the rotation matrix corresponding to the quarternion
if length(q) ~= 4
    error('Error: invalid quaternion.');
end
R = [2*q(1).^2-1+2*q(2)^2, 2*(q(2)*q(3)+q(1)*q(4)), 2*(q(2)*q(4)-q(1)*q(3)); 
    2*(q(2)*q(3)-q(1)*q(4)), 2*q(1)^2-1+2*q(3)^2, 2*(q(3)*q(4)+q(1)*q(2)); 
    2*(q(2)*q(4)+q(1)*q(3)), 2*(q(3)*q(4)-q(1)*q(2)), 2*q(1)^2-1+2*q(4)^2];
end
