function [axisC,axisX,axisZ,varargout] = moore650ikine(varargin)
%IKINE inverse kinematics for moore 650 FG V2

unitList = {'m','mm','\mum','nm'};
inputUnit = varargin{end - 1};
outputUnit = varargin{end};
presUnit = find(strcmp(unitList,inputUnit),1);
aimUnit = find(strcmp(unitList,outputUnit),1);

switch nargin
    case 4 % spiralAngle spiralPath inputUnit outputUnit
        spiralAngle = varargin{1};
        spiralPath = varargin{2};
        spiralAngle1 = 180/pi*spiralAngle; % rad -> deg
        spiralPath1 = 1000^(aimUnit - presUnit)*spiralPath; % um -> mm

%         axisC = atan2(2*(spiralQuat(:,2).*spiralQuat(:,3) - spiralQuat(:,1).*spiralQuat(:,4)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,3).^2) - 1);
%         axisB = atan2(-2*(spiralQuat(:,2).*spiralQuat(:,4) + spiralQuat(:,1).*spiralQuat(:,3)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,4).^2) - 1);
%         axisZ = (spiralPath1(3,:).')./cos(axisB);
%         axisX = axisZ.*sin(axisB).*cos(axisC) - spiralPath1(1,:).';
%         axisY = axisZ.*sin(axisB).*sin(axisC) - spiralPath1(2,:).';

        axisC = (wrapTo360(spiralAngle1));
        axisZ = spiralPath1(3,:);
        axisX = sign(spiralPath1(1,1))*vecnorm(spiralPath1(1:2,:),2,1);
    case 5 % axisC axisX axisZ inputUnit outputUnit
        axisC = pi/180*varargin{1};
        axisX = 1000^(aimUnit - presUnit)*varargin{2};
        axisZ = 1000^(aimUnit - presUnit)*varargin{3};
end

% zero point of the C axis
if axisC(1) ~= 0
    axisC = axisC - axisC(1);
    axisC = wrapTo360(axisC);
end

% zero point of z axis
if axisZ(end) ~= 0
    axisZ = axisZ - axisZ(end);
end

% direction correction
%         if strcmp(startDirection,'X Plus')
% axisX = -1.*axisX;
% axisC = wrapTo360(-1.*axisC);
axisC(find(abs(axisC - 360) < 1e-3)) = 0;

if nargin == 4
    varargout{1} = spiralAngle1;
    varargout{2} = spiralPath1;
end

end