function q = slerp(q1,q2,t)
% To solve the interpolation process of vectors / quaternion based on 
% spherical linear interpolation, short for slerp.

q1 = q1./vecnorm(q1,2,1);
q2 = q2./vecnorm(q2,2,1);
if size(q1,1) == 1
    cosa = dot(q1,q2);
    if cosa < 0
        q2 = -q2;
        cosa = -cosa;
    end
    if cosa > cosd(5)
        q = (1 - t)*q1 + t*q2;
        return;
    end
    sina = sqrt(1 - cosa.^2);
    a = atan2(sina,cosa);
    k1 = sin((1 - t).*a)./sina;
    k2 = sin(t.*a)./sina;
    
    q(cosa <= cosd(5),:) = k1*q1 + k2*q2;
else
    cosa = dot(q1,q2,2);
    q2(cosa<0,:) = -q2(cosa<0,:);
    q = zeros(size(q1,1),4);
    q(cosa > cosd(5),:) = (1 - t)*q1 + t*q2;

    cosa = -cosa;
    sina = sqrt(1 - cosa.^2);
    a = atan2(sina,cosa);
    k1 = sin((1 - t).*a)./sina;
    k2 = sin(t.*a)./sina;
    
    q(cosa <= cosd(5),:) = k1*q1 + k2*q2;
end

end