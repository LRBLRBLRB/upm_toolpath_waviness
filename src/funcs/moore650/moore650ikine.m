function cncData = moore650ikine( ...
    cncMode,spiralAngle,spiralPath,options)
%IKINE inverse kinematics for moore 650 FG V2
%
% Notice: the unit of rotation axes must be converted to degree!!!

arguments
    cncMode {mustBeMember(cncMode,{'CXZ'})}
    spiralAngle (1,:)
    spiralPath (3,:)
    options.cStart logical = false
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

cncData = cncpreprocess(cncMode,cncData,'cStart',options.cStart);

end