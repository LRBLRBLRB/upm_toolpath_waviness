function tool2 = toolRigid(tool1,Rot,Vec)
% solve the rigid transform of the tool struct

arguments
    tool1
    Rot (3,3) {mustBeFinite}
    Vec (3,1) {mustBeFinite}
end

% H = [Rot,Vec;0,0,0,1];
tool2.center = Rot*tool1.center + Vec;
tool2.toolEdgeNorm = Rot*tool1.toolEdgeNorm; % (3,:) normal vector of the tool edge
tool2.toolDirect = Rot*tool1.toolDirect; % (3,1) cutting direction of the tool edge
% tool2.Bform.coefs = Rot(1:2,1:2)*tool1.Bform.coefs + Vec(1:2);
% tool2.toolPt = Rot*tool1.toolPt + Vec; % (3,:) position of tool edge

end