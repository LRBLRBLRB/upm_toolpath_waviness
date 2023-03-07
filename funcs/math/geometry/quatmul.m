function quat12 = quatmul(quat1,quat2)
%QUATMUL multiplication of two quaternions
%   此处显示详细说明

num = size(quat1,1);
quat12 = zeros(num,4);

quat12(:,1) = quat1(:,1).*quat2(:,1) - quat1(:,2:4)*(quat2(:,2:4)).';
quat12(:,2:4) = quat1(:,1).*quat2(:,2:4) + quat2(:,1).*quat1(:,2:4) ...
    + cross(quat1(:,2:4),quat2(:,2:4),2);

end

