function [tool1,tool2] = toolAdj(tool1,tool2,toolCenPt)

% 输入刀路和刀具模型，输出的两个刀具模型中，节点参数做了限制和张角都分别重新做了限制，从而实现删除片段

