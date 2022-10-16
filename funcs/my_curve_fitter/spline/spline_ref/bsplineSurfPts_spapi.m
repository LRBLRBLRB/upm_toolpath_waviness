function [ptsMesh,varargout] = bsplineSurfPts_spapi(Q,k,l,u,v,options)
% to generate the B-spline surface with interpolation of the given mesh
% points
% Usage:
%   


arguments
    Q (:,:,3) {mustBeFinite} 
    k {mustBePositive}
    l {mustBePositive}
    u (:,1)
    v (:,1)
    options.paramMethod {mustBeMember(options.paramMethod, ...
        ['uniform','centripetal','chord',''])} = 'chord'
    options.cptsType {mustBeMember(options.cptsType, ...
        ['Cartesian','Polar','Spherical',''])} = 'Cartesian'
end

[m,n,~] = size(Q);
% test if the cpts are enough to generate the quasi-uniform B-spline
if m-k-1<0 || n-l-1<0
    error("control points are not enough, or the required order is too large");
end

switch options.cptsType
    case 'Cartesian'
        Qcart = Q;
    case 'Polar'
        Qcart = Q;
        Qcart(:,:,1) = sqrt(Q(:,:,1).^2 + Q(:,:,2).^2);
        Qcart(:,:,2) = atan2(Q(:,:,2),Q(:,:,1));
    case 'Spherical'
        Qcart(:,:,1) = sqrt(Q(:,:,1).^2 + Q(:,:,2).^2 + Q(:,:,3).^2); % polar radius
        Qcart(:,:,2) = atan2(Q(:,:,2),Q(:,:,1)); % longitude
        Qcart(:,:,3) = atan2(sqrt(Q(:,:,1).^2 + Q(:,:,2).^2),Q(:,:,3)); % latitude
end

% solve node params corresponding to Q
[uQ,vQ] = interpParamSurf(Qcart,'paramMethod',options.paramMethod);

% node vector U & V generation
U = nodeVector(k,m,'uQ',uQ,'nodeMethod','Interpolation');
V = nodeVector(l,n,'uQ',vQ,'nodeMethod','Interpolation');

% B-spilne interpolation
Qtmp = permute(Q,[3,1,2]); % to use the function spapi, Q must be reshaped as 3*m*n
sp = spapi({U,V},{uQ,vQ},Qtmp);
ptsMesh = fnval(sp,{u,v});
varargout{1} = sp;
end