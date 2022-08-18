function [pts,sp] = bsplinePts_spap2(Q,k,u,options)
% usage: [pts,cpts] = bsplinePts_spapi(Q,k,u,options)
%        [pts,cpts,sp] = bsplinePts_spapi(Q,k,u,options)
% to approximate the given points using segmented B-spline curves in MATLAB
% curve fitting toolbox.
% Inputs:
%   Q (n,3) or (n,2) the points array that is to be fitted
%   k (1,1) the degree of the B-spline curves
%   u (:,1) the parameter vector of the spline curve to be solved after
%   approximation
%   options:
%       paramMethod {'uniform','concentric','chord',''} the method of
%       parameterization during curve approximation
% Outputs:
%   pts (2,:) or (3,:) the points on the spline curve to be solved after
%   approximation
%   varargout{1} struct the struct of the spline, the form of which
%   remains the B-form.

arguments
    Q 
    k {mustBePositive}
    u (:,1)
    options.paramMethod {mustBeMember(options.paramMethod, ...
    ['uniform','concentric','chord',''])} = 'chord'
end

if size(Q,2) == 2 || size(Q,2) == 3
    Q = Q';
end
[~,n] = size(Q);
if n-k-1<0
    error("control points are not enough, or the required degree is too large");
end

% chord parameterization
uQ = approxParam(Q',"paramMethod",options.paramMethod);

% node vector U generation
U = nodeVector(k,n,'nodeMethod','Interpolation','uQ',uQ);

% 拟合的U可以微调？
sp = spap2(U,k+1,uQ,Q); % spline fitting
pts = fnval(sp,u);

end