function [ind1,ind2] = getInnerLoopToolPathIndex(angle,angleDel)
% to get the index of the closest tool path point of the current one within
% the inner loop
%
% Inputs:
%   angle    (1,n) the concentric angle of the points on the closest loop
%   angleDel (1,n) the difference between angle and the aimed point

angleDel1 = angleDel;
if isempty(angle(angleDel1 >= 0))
    % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
    angleDel1 = angleDel1 + 2*pi;
end
ind1 = find(angle == min(angle(angleDel1 >= 0)));

angleDel2 = angleDel;
if isempty(angle(angleDel2 < 0))
    angleDel2 = angleDel2 - 2*pi;
end
ind2 = find(angle == max(angle(angleDel2 < 0)));


end