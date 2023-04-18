function [circ3D,scatterDst,varargout] = toolFit3D(scatterOri,options)
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
    options.arcFitMethod {mustBeMember(options.arcFitMethod, ...
        ['gradient-decent','normal-equation','levenberg-marquardt'])} ...
        = 'levenberg-marquardt'
    options.arcFitdisplayType {mustBeMember(options.arcFitdisplayType, ...
        ['off','none','iter','iter-detailed','final','final-detailed'])} = 'final'
end

[circ3D,RMSE,scatterPlane,circ2D] = arcFit3D(scatterOri, ...
    'arcFitMethod',options.arcFitMethod,'arcdisplayType',options.arcFitdisplayType);

% debug
% hold on;
% scatter3(circ3D{1}(1),circ3D{1}(2),circ3D{1}(3));
% plot3([circ3D{1}(1),circ3D{1}(1) + circ3D{2}*circ3D{4}(1)], ...
%     [circ3D{1}(2),circ3D{1}(2) + circ3D{2}*circ3D{4}(2)], ...
%     [circ3D{1}(3),circ3D{1}(3) + circ3D{2}*circ3D{4}(3)]);
% RLeft = vec2rot(circ3D{3}/2*circ3D{5});
% tmpLeft = RLeft*circ3D{4};
% plot3([circ3D{1}(1),circ3D{1}(1) + circ3D{2}*tmpLeft(1)], ...
%     [circ3D{1}(2),circ3D{1}(2) + circ3D{2}*tmpLeft(2)], ...
%     [circ3D{1}(3),circ3D{1}(3) + circ3D{2}*tmpLeft(3)]);
% RRight = vec2rot(-circ3D{3}/2*circ3D{5});
% tmpRight = RRight*circ3D{4};
% plot3([circ3D{1}(1),circ3D{1}(1) + circ3D{2}*tmpRight(1)], ...
%     [circ3D{1}(2),circ3D{1}(2) + circ3D{2}*tmpRight(2)], ...
%     [circ3D{1}(3),circ3D{1}(3) + circ3D{2}*tmpRight(3)]);

%% rigid transform: to get the standardized tool profile
n = size(scatterOri,2);
rotAng = -pi/2 - atan2(circ2D.arcVec(2),circ2D.arcVec(1));
rotMat = rotz(rotAng);
rotMat = rotMat(1:2,1:2);
scatterDst = rotMat*(scatterPlane(1:2,:) - ndgrid(circ2D.center,1:n));

%% optional output
varargout{1} = RMSE;

end