function [interPt,I1,I2,eps] = pcCrossDN(s1,s2)
% To search the nearest point within the two point sets
% 
% Usage:
%   


s1 = s1'; % 3*n -> n*3
s2 = s2'; % 3*n -> n*3

% index: the ith element represents the index of the point in s1 that is
% nearest to the ith point in s2
% dist: the distance between the nearest point in s1 and every point in s2
[index,dist] = dsearchn(s1,s2);
[eps,I2] = min(dist);
I1 = index(I2);
interPt = (0.5*(s1(I1,:) + s2(I2,:)))';
end