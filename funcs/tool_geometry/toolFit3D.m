function [r,ang,scatterDst,varargout] = toolFit3D(scatterOri,options)
% usage: [c,r,ang,scatterDst,RMSE] = toolFit3D(scatterOri,options)
%
% solve the edge sharpness of a arc turning tool: 
% calculate the radius of the tool tip arc, 
% as well as transforming the tool coordinate
%
% Inputs: 
%   scatterMat: original tool edge point from measuring (3,n)
%   options: 
%       fitMethod: the method of curve fitting
% Outputs: 
%   r: radius of the arc (1,1)
%   ang: the open angle of the tool (1,1)
%   scatterDst: tool edge points with pose adjustment (3,n)
%   c: center of the arc (3,1)
%   RMSE: the square-mean-root error of curve fitting
%
% method:
%   least square method by normal equation solving

arguments
    scatterOri (3,:) {mustBeFinite}
    options.fitMethod {mustBeMember(options.fitMethod, ...
        ['Gradient-decent','Normal-equation','Levenberg-Marquardt'])} ...
        = 'Levenberg-Marquardt'
    options.displayType {mustBeMember(options.displayType, ...
        ['off','none','iter','iter-detailed','final','final-detailed'])} = 'final'
end

[circ3D,RMSE,scatterPlane,c,startV,endV] = arcFit3D(scatterOri, ...
    'fitMethod',options.fitMethod,'displayType',options.displayType);
r = circ3D{2};
ang = circ3D{3};

%% rigid transform: to get the standardized tool profile
n = size(scatterOri,2);
mid = 0.5*(startV + endV);
rotAng = pi/2 - atan2(mid(2),mid(1));
rotMat = rotz(rotAng);
rotMat = rotMat(1:2,1:2);
scatterDst = rotMat*(scatterPlane(1:2,:) - ndgrid(c,1:n));

%% optional output
varargout{1} = RMSE;

end