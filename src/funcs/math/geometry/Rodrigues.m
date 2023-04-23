function v_rot = Rodrigues(v,k,theta,eps)
% usage: 
% to solve the v_rot, i.e., rotate vector v by theta around k axis
% inputs:
%   v (:,3) original vectors
%   k (:,3) rotation axis, the length of which remains 1
%   theta (:,1) rotation angle in radian, which satisfies right-hand rule
%   eps (1,1)
% outputs:
%   v_rot (:,3) vectors after rotating

arguments
    v (:,3) {mustBeFinite}
    k (:,3) {mustBeFinite}
    theta (:,1) {mustBeFinite}
    eps (1,1) {mustBeNonnegative} = 1e-2
end

% to ensure that every axis vector is unit vector
if any(abs(vecnorm(k,2,2) - 1) > eps,"all")
    ind = find(abs(vecnorm(k,2,2) - 1) > eps);
    error(['Error: axis No.',ind,' remains a non-unit vector']);
end

numV = size(v,1);
numK = size(k,1);
if numV ~= 1 && numK == 1
    k = meshgrid(k,1:numV);
elseif numV ~= numK
    error('Error: the size of vector should be the same of that of axis.')
end

% rotated vector (Rodrigues' formula)
v_rot = v.*cos(theta) ...
    + cross(k,v,2).*sin(theta) ...
    + k.*(dot(k,v,2)).*(1 - cos(theta)); 
end
