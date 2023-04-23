function varargout = spiralcontactu(spiralPath,spiralQuat,var1,var2,var3,options)
%SPIRALCONTACTU calculate the contact u parameter from the tool to the surface 
%   此处显示详细说明

arguments
    spiralPath (3,:)
    spiralQuat (:,4)
    var1
    var2
    var3 = []
    options.method {mustBeMember(options.method, ...
        {'distance-to-surf','linear-interp'})} = 'distance-to-surf'
    options.interpType {mustBeMember(options.interpType, ...
        {'linear','nearest','next','previous','pchip','cubic,' ...
        'quadratic'})} = 'linear'
    options.uQTol {mustBePositive} = 1e-2
end

switch options.method
    case 'distance-to-surf'
        toolSp = var1;
        surfFunc = var2;
        % tool contact u calculation
        spiralNum = size(spiralPath,2);
        uDist = zeros(1,spiralNum);
        toolContactU = zeros(1,spiralNum);
        for kk = 1:spiralNum
            uQ = linspace(0,1,ceil(1/options.uQTol));
            toolSp1 = toolSp;
            toolSp1.coefs = quat2rotm(spiralQuat(kk,:))*toolSp.coefs + spiralPath(:,kk);
            uVal = fnval(toolSp1,uQ);
            dist = dist2surf(uVal,surfFunc, ...
                'CalculateType','Lagrange-Multiplier', ...
                'DisplayType','none');
            [uDist(kk),toolContactU(kk)] = min(dist);
        end
    case 'linear-interp'
        % inplemented in the function arclengthparam.m
%         spiralAngle0 = var1;
%         spiralAngle = var2;
%         maxAng = var3;
%         angDist = abs((spiralAngle0(2:end) - spiralAngle0(1:end - 1)) - maxAng);
%         options.uQTol = 1e-6; % the epsilon to judge whether the angle difference is equal or not
%         angularMaxInd = find(angDist > options.uQTol,1);
%         % tool contact u calculation with linear interpolation
%         spiralNum = size(spiralPath,2);
%         toolContactU = zeros(1,spiralNum);
%         for ii = angularMaxInd + 1:numEq
%             ind1 = find(spiralAngle(ii) >= spiralAngle0);
%             ind1 = ind1(end); % the closest parameter smaller than tEq
%             ind2 = find(spiralAngle(ii) < spiralAngle0,1); % the closest parameter larger than tEq
% 
%         end

end

if nargout == 1
    varargout{1} = toolContactU;
else
    varargout{1} = toolContactU;
    varargout{2} = uDist;
end

end

