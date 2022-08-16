function [pts,sp] = bsplinePts_spapi(Q,k,u)

arguments
    Q 
    k {mustBePositive}
    u (:,1)
end

if size(Q,2) == 2 || size(Q,2) == 3
    Q = Q';
end
[~,n] = size(Q);
if n-k-1<0
    error("control points are not enough, or the required order is too large");
end

% chord parameterization
uQ = interpParam(Q','chord');
U = nodeVector(k,n,'nodeMethod','Interpolation','uQ',uQ);

sp = spapi(U,uQ,Q); % spline interpolation
pts = fnval(sp,u); % some questions here

end