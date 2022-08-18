function [pts,varargout] = bsplinePts_spapi(Q,k,u,options)
% usage: [pts,cpts] = bsplinePts_spapi(Q,k,u,options)
%        [pts,cpts,sp] = bsplinePts_spapi(Q,k,u,options)
% to interpolate the given points using segmented B-spline curves in MATLAB
% curve fitting toolbox.
% Inputs:
%   Q (n,3) or (n,2) the points array that is to be interpolated
%   k (1,1) the degree of the B-spline curves
%   u (:,1) the parameter vector of the spline curve to be solved after interpolation
%   options:
%       paramMethod {'uniform','centripetal','chord',''} the method of
%       parametertization during interpolation
% Outputs:
%   pts (2,:) or (3,:) the points on the spline curve to be solved after
%   interpolation
%   varargout{1} struct the struct of the spline, the form of which
%   remains the B-form.

arguments
    Q 
    k {mustBePositive}
    u (:,1)
    options.paramMethod {mustBeMember(options.paramMethod, ...
        ['uniform','centripetal','chord',''])} = 'chord'
end

if size(Q,2) == 2 || size(Q,2) == 3
    Q = Q';
end
[~,n] = size(Q);
if n-k-1<0
    error("control points are not enough, or the required degree is too large");
end

% chord parameterization
uQ = interpParam(Q','paramMethod',options.paramMethod);

% node vector U generation
U = nodeVector(k,n,'nodeMethod','Interpolation','uQ',uQ);

sp = spapi(U,uQ,Q); % spline interpolation
pts = fnval(sp,u); % some questions here
varargout{1} = sp;
end