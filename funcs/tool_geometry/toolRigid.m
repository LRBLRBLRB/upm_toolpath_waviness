function tool2 = toolRigid(tool1,varargin)
% solve the rigid transform of the tool struct
%
% usage: 
% tool2 = toolRigid(tool1,Rot,Vec)
%   solve the rigid transform of the tool edge, satisfying the fomula
%   tool2 = Rot*tool1 + Vec;
% tool2 = toolRigid(tool1,……

Rot = varargin{1};
Vec = varargin{2};
if any(size(Rot) - [3,3]) || any(size(Vec) - [3,1])
    error('Incorrect dimension of the rotation matrix or the traslation vector.');
end

% H = [Rot,Vec;0,0,0,1];
tool2 = tool1;
tool2.center = Rot*tool1.center + Vec;
tool2.toolEdgeNorm = Rot*tool1.toolEdgeNorm; % (3,:) normal vector of the tool edge
tool2.cutDirect = Rot*tool1.cutDirect; % (3,1) cutting direction of the tool edge
tool2.toolDirect = Rot*tool1.toolDirect; % (3,1) tool edge direction of the tool edge
tool2.toolBform.coefs = Rot*tool1.toolBform.coefs + Vec;
% tool2.toolPt = Rot*tool1.toolPt + Vec; % (3,:) position of tool edge

end