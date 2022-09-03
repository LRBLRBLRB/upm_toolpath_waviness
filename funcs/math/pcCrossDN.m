function [interPt,I1,I2,eps] = pcCrossDN(s1,s2)
s1 = s1';
s2 = s2';
[index,dist] = dsearchn(s1,s2);
[eps,I2] = min(dist);
I1 = index(I2);
interPt = (0.5*(s1(I1,:) + s2(I2,:)))';
end