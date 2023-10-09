function [spiralAngle,spiralPath,spiralQuat,spiralNorm,spiralCut] = moore650kine(cncMode,cncData, ...
    inputLengthUnit,outputLengthUnit,inputAngleUnit,outputAngleUnit,options)
%MOORE650KINE kinematics of Moore Nanotech 650FG V2

arguments
    cncMode {mustBeMember(cncMode,{'CXZ'})}
    cncData double
    inputLengthUnit {mustBeMember(inputLengthUnit,{'m','mm','\mum','nm'})}
    outputLengthUnit {mustBeMember(outputLengthUnit,{'m','mm','\mum','nm'})}
    inputAngleUnit {mustBeMember(inputAngleUnit,{'deg','rad',''})} = []
    outputAngleUnit {mustBeMember(outputAngleUnit,{'deg','rad',''})} = []
    options.toolData struct = []
    options.shift logical = false
end

spiralPtNum = size(cncData,2);
spiralPath = zeros(3,spiralPtNum);
spiralAngle = zeros(1,spiralPtNum);
spiralQuat = zeros(spiralPtNum,4);
switch cncMode
    case 'CXZ'
        spiralPath(:,1) = [cncData(2,1);0;cncData(3,1)];
        spiralQuat(1,:) = [1,0,0,0];
        for ii = 1:spiralPtNum
            spiralAngle(ii) = cncData(1,ii);
            spiralQuat(ii,:) = rotm2quat(rotz(180/pi*spiralAngle(ii)));
            spiralPath(1,ii) = cncData(2,ii)*cos(spiralAngle(ii));
            spiralPath(2,ii) = cncData(2,ii)*sin(spiralAngle(ii));
            spiralPath(3,ii) = cncData(3,ii);
        end
        for ii = 2:spiralPtNum % only "X Minus" & "Edge-to-Center" is testified
            if (spiralAngle(ii) - spiralAngle(ii - 1))*(conThetaBound(2) - conThetaBound(1)) < 0
                spiralAngle(ii:end) = spiralAngle(ii:end) + conThetaBound(2) - conThetaBound(1);
            end
        end
end

%% spiral normal vector & cut direction vector
if nargout == 5 && ~isempty(options.toolData)
    spiralNorm = zeros(3,spiralPtNum);
    spiralCut = zeros(3,spiralPtNum);
    for ii = 1:spiralPtNum
            spiralNorm(:,ii) = quat2rotm(spiralQuat(ii,:))*toolData.toolEdgeNorm;
            spiralCut(:,ii) = quat2rotm(spiralQuat(ii,:))*toolData.cutDirect;
    end
end

if shift
    if isempty(options.toolData)
        warning("No enough input values for tool path shifting, options.toolData is needed.\n");
        return;
    else
        toolSp = toolData.toolBform;
        toolSp.coefs = quat2rotm(spiralQuat(end,:))*toolSp.coefs + spiralPath(:,end);
        toolSpPt = fnval(toolSp,0:0.0001:1);
        spiralPath(3,:) = spiralPath(3,:) - min(toolSpPt(3,:));
    end
end