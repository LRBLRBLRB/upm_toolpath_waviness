function [c,r,ang,scatterDst,RMSE] = toolFit2D(scatterOri,fitMethod)
% usage: [c,r,ang,scatterDst,RMSE] = toolFit2D(scatterOri)
%
% solve the edge sharpness of a arc turning tool: 
% calculate the radius of the tool tip arc, 
% as well as transforming the tool coordinate
%
% Inputs: 
%   scatterMat: original tool edge point from measuring (n,2)
%   options: 
%       fitMethod: the method of curve fitting
% Outputs: 
%   r: radius of the arc (1,1)
%   ang: the open angle of the tool (1,1)
%   scatterDst: tool edge points with pose adjustment (n,2)
%   c: center of the arc (1,2)
%   RMSE: the square-mean-root error of curve fitting
% method:
%   least square method by normal equation solving

arguments
    scatterOri (:,2) {mustBeFinite}
    fitMethod {mustBeMember(fitMethod, ...
        ['Gradient-decent','Normal-equation','Levenberg-Marquardt'])} ...
        = 'Levenberg-Marquardt'
end

n = size(scatterOri,1);

%% circle fitting
[c,r,ang,RMSE,startV,endV] = circleFit2D(scatterOri,fitMethod);

%% rigid transform
mid = 0.5*(startV + endV);
rotAng = pi/2 - atan2(mid(2),mid(1));
rotMat = [cos(rotAng), -sin(rotAng);
          sin(rotAng),  cos(rotAng)];
scatterDst = (scatterOri - meshgrid(c,1:n))*(rotMat');

end