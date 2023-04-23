function [ind1,ind2] = getOuterLoopToolPathIndex(angle,angleDel,loopPtNumLast)
% to get the index of the closest tool path point of the current one within
% the outer loop

if isempty(angle(angleDel >= 0))
    % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
    angleDel = angleDel + 2*pi;
end
ind1 = loopPtNumLast + find(angle == min(angle(angleDel >= 0)));
if isempty(angle(angleDel < 0))
    angleDel = angleDel - 2*pi;
end
ind2 = loopPtNumLast + find(angle == max(angle(angleDel < 0)));
end