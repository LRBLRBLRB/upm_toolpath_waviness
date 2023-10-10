function [cncData,spiralAngle,spiralPath] = moore650ikine( ...
    cncMode,spiralAngle,spiralPath,options)
%IKINE inverse kinematics for moore 650 FG V2
%
% Notice: the unit of rotation axes must be converted to degree!!!

arguments
    cncMode {mustBeMember(cncMode,{'CXZ'})}
    spiralAngle (3,:)
    spiralPath (3,:)
    options.cStart double = []
    options.shift logical = false
end

axNum = strlength(cncMode);
ptNum = size(spiralAngle,2);
cncData = zeros(axNum,ptNum);

%         axisC = atan2(2*(spiralQuat(:,2).*spiralQuat(:,3) - spiralQuat(:,1).*spiralQuat(:,4)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,3).^2) - 1);
%         axisB = atan2(-2*(spiralQuat(:,2).*spiralQuat(:,4) + spiralQuat(:,1).*spiralQuat(:,3)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,4).^2) - 1);
%         axisZ = (spiralPath1(3,:).')./cos(axisB);
%         axisX = axisZ.*sin(axisB).*cos(axisC) - spiralPath1(1,:).';
%         axisY = axisZ.*sin(axisB).*sin(axisC) - spiralPath1(2,:).';

switch cncMode
    case 'CXZ'
        cncData(1,:) = (wrapTo360(spiralAngle));
        cncData(3,:) = spiralPath(3,:);
        cncData(2,:) = sign(spiralPath(1,1))*vecnorm(spiralPath(1:2,:),2,1);
end

% zero point of the C axis
if ~isempty(options.cStart)
    cInd = strfind(cncMode,'C');
    cncData(cInd,:) = cncData(cInd,:) - cncData(cInd,1);
    cncData(cInd,:) = wrapTo360(cncData(cInd,:));
end

% zero point of z axis
zInd = strfind(cncMode,'C');
if cncData(zInd,end) ~= 0
    cncData(zInd,:) = cncData(zInd,:) - cncData(zInd,end);
end

% direction correction
%         if strcmp(startDirection,'X Plus')
% axisX = -1.*axisX;
% axisC = wrapTo360(-1.*axisC);
cncData(cInd,abs(cncData(cInd,:) - 360) < 1e-3) = 0;

end