function [circ2D,RMSE] = arcFit2D(scatterOri,param0,options)
% usage: [c,r,ang,RMSE,startV,endV] = arcFit2D(scatterOri,options)
%
% fit the circle from a cluster of 2D points
%
% Inputs: 
%   scatterMat: original 2D points (2,n)
%   fitMethod: the method of curve fitting
% Outputs: 
%   r: radius of the arc (1,1)
%   ang: the open angle of the tool (1,1)
%   c: center of the arc (2,1)
%   RMSE: the square-mean-root error of curve fitting
%   startV: start edge vector of the arc
%   endV: end edge vector of the arc
%
% method:
%   least square method by normal equation solving

arguments
    scatterOri (2,:) {mustBeFinite}
    param0
    options.arcFitMethod {mustBeMember(options.arcFitMethod, ...
        {'gradient-decent','normal-equation','levenberg-marquardt'})} ...
        = 'levenberg-marquardt'
    options.displayType {mustBeMember(options.displayType, ...
        ['off','none','iter','iter-detailed','final','final-detailed'])} = 'final'
end

x = (scatterOri(1,:))';
y = (scatterOri(2,:))';
n = length(x);

%% circle fitting
switch options.arcFitMethod
    case 'gradient-decent' %% method one: gradient decent
        M = [sum(x.^2), sum(x.*y), sum(x);
             sum(x.*y), sum(y.^2), sum(y);
             sum(x),    sum(y),    n];
        b = [sum(x.*x.*x) + sum(x.*y.*y);
             sum(x.*x.*y) + sum(y.*y.*y);
             sum(x.*x)    + sum(y.*y)];
        param = -M\b;
    case 'normal-equation' %% method two: nomal equation 
        M = [x,y,ones(n,1)];
        b = -x.^2 - y.^2;
        param = (M'*M)\(M'*b);
        % The above two ways gives the same result.

    case 'levenberg-marquardt' %% method three: Levenberg-Marquardt
        if startsWith(options.displayType,'i') || startsWith(options.displayType,'f')
            fprintf('The arc is fitted with the method [%s].\n',options.arcFitMethod);
        end
        F = @(p,x) x(:,1).^2 + x(:,2).^2 + p(1)*x(:,1) + p(2)*x(:,2) + p(3);
        optimOpt = optimoptions(...
            'lsqcurvefit',...
            'Algorithm',options.arcFitMethod, ...
            'MaxIterations',2000, ...
            'Display',options.displayType);
        lb = [];
        ub = [];
        param = lsqcurvefit(F,param0,scatterOri',zeros(n,1),lb,ub,optimOpt);
        % lsqnonlin lsqcurvefit

        %% accuracy
        F = x.^2 + y.^2 + param(1)*x + param(2)*y + param(3);
        MSE = sum(F.^2)/n;
        RMSE = sqrt(MSE);
end

%% params of the arc
circ2D.center = zeros(2,1);
circ2D.center(1) = -param(1)/2;
circ2D.center(2) = -param(2)/2;
circ2D.radius = 0.5*sqrt(param(1)^2 + param(2)^2 - 4*param(3));
if imag(circ2D.radius) ~= 0
    error('Error: cannot fit a circle!');
end
circ2D.startV = scatterOri(:,1) - circ2D.center;
circ2D.endV = scatterOri(:,end) - circ2D.center;
circ2D.openAng = vecAng(circ2D.startV,circ2D.endV,1);
circ2D.arcVec = 0.5*(circ2D.startV/norm(circ2D.startV) + circ2D.endV/norm(circ2D.endV));
circ2D.arcVec = circ2D.arcVec/norm(circ2D.arcVec);

end
