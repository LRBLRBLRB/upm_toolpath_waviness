function [ind1,ind2] = getInnerLoopToolPathIndex(angle,angleDel)
% to get the index of the closest tool path point of the current one within
% the inner loop

if isempty(angle(angleDel >= 0))
    % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
    angleDel = angleDel + 2*pi;
end
ind1 = find(angle == min(angle(angleDel >= 0)));
if isempty(angle(angleDel < 0))
    angleDel = angleDel - 2*pi;
end
ind2 = find(angle == max(angle(angleDel < 0)));
end