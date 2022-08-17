function [pts,cpts,varargout] = bsplinePts_spapi(Q,k,u,options)

arguments
    Q 
    k {mustBePositive}
    u (:,1)
    options.interpMethod {mustBeMember(options.interpMethod, ...
        ['uniform','concentric','chord',''])} = 'chord'
end

if size(Q,2) == 2 || size(Q,2) == 3
    Q = Q';
end
[~,n] = size(Q);
if n-k-1<0
    error("control points are not enough, or the required order is too large");
end

% chord parameterization
uQ = interpParam(Q','interpMethod',options.interpMethod);

% node vector U generation
U = nodeVector(k,n,'nodeMethod','Interpolation','uQ',uQ);

sp = spapi(U,uQ,Q); % spline interpolation
cpts = sp.coefs;
pts = fnval(sp,u); % some questions here
varargout{1} = sp;
end