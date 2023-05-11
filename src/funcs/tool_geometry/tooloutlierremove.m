function tool = tooloutlierremove(toolOri)
%TOOLINI 此处显示有关此函数的摘要
%   此处显示详细说明

% toolOri (2,:)

toolAng = atan2(toolOri(2,:),toolOri(1,:));

diffAng = diff(toolAng);
diffAng = abs(diffAng);

standard = abs(diffAng(1));

% the index of the pt which is the reverse one in the whole tooldata
duplicateInd = find(diffAng > 5*standard,1) + 1;

if isempty(duplicateInd)
    tool = toolOri;
    return;
end

% the index of the pt which is the last one to be resumed in the left part
resumeInd = find(toolAng(1:duplicateInd - 1) - toolAng(duplicateInd) < 0,1) - 1;

% toolOri(:,resumeInd:duplicateInd - 1) = [];

% the right part, with a offset to connect with th left
% toolLeftInd = resumeInd - 1;
% toolRightInd = duplicateInd;


tool = [toolOri(:,1:resumeInd - 1) ...
    + ndgrid(toolOri(:,duplicateInd) - toolOri(:,resumeInd - 1),1:resumeInd - 1), ...
    toolOri(:,duplicateInd:end)];

end