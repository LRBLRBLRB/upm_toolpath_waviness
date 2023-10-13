function cncData = cncpreprocess(cncMode,cncData,options)
%CNCPREPROCESS 此处显示有关此函数的摘要
%   此处显示详细说明

arguments
    cncMode {mustBeMember(cncMode,{'CXZ'})}
    cncData double
    options.cStart logical = false
end

% zero point of the C axis
if options.cStart
    cInd = strfind(cncMode,'C');
    cncData(cInd,:) = cncData(cInd,:) - cncData(cInd,1);
    cncData(cInd,:) = wrapTo360(cncData(cInd,:));
    cncData(cInd,abs(cncData(cInd,:) - 360) < 1e-3) = 0;
end

% zero point of z axis
zInd = strfind(cncMode,'Z');
if cncData(zInd,end) ~= 0
    cncData(zInd,:) = cncData(zInd,:) - cncData(zInd,end);
end

% direction correction
%         if strcmp(startDirection,'X Plus')
% axisX = -1.*axisX;
% axisC = wrapTo360(-1.*axisC);

end

