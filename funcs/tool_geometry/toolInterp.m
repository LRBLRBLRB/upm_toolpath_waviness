function [toolPtInterp,toolUInterp,varargout] = toolInterp( ...
    interpR,toolR,ind1,ind2,ind3,toolPt,toolContactU,var1,var2)
% to linear interpolate the exact tool path between two known tool points
%
% Usage:
%
% [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp(
%   toolPt,toolNorm,toolCutDir,toolU,ind1,ind2,indInterp)
% Inputs:
%   toolPt          (3,:) double    the list of tool path points
%   toolNorm        (3,:) double    the list of the spindle direction
%   toolCutDir      (3,:) double    the list of the cutting direction
% Outputs:
%   toolPtInterp    (3,1) double    the current point on the tool path
%   toolNormInterp  (3,1) double    the spindle direction of the current point
%   toolCutDirInterp(3,1) double    the cutting direction of the current point
%
% [toolPtInterp,toolNormInterp,toolCutDirInterp] = toolInterp(
%   toolPt1,toolNorm1,toolCutDir1,toolU1,toolPt2,toolNorm2,toolCutDir2,toolU2,toolUInterp)
% Outputs: same as the above

if size(var1,2) == 4
    toolQuat = var1;
    toolCutDir1 = var2;
    method = 1;
elseif  size(var2,2) == 4
    toolCutDir1 = var1;
    toolQuat = var2;
    method = 1;
else
    toolCutDir1 = var2(:,ind1);
    % toolCutDir2 = var2(:,ind2);
    % toolCutDir3 = var2(:,ind3);
    method = 2;
end

%% tool point projection & interpolation
toolPt1 = toolPt(:,ind1); % the outer tool point
toolPt2 = toolPt(:,ind2); % the 1st inner tool point
toolPt3 = toolPt(:,ind3); % the 2nd inner tool point
toolR1 = toolR(ind1);
toolR2 = toolR(ind2);
toolU1 = toolContactU(ind1);
toolU2 = toolContactU(ind2);
toolU3 = toolContactU(ind3);

toolPtProj = lineIntersectPlane(toolPt2,toolPt3,toolPt1,toolCutDir1);
t = norm(toolPtProj - toolPt2)/norm(toolPt3 - toolPt2);
toolUProj = (1 - t)*toolU2 + t*toolU3;

toolRInterp = (toolR1 - interpR)/(toolR1 - toolR2); % p.s. toolU2 = toolU3 = toolUProj
toolPtInterp = toolPt1 - toolRInterp*(toolPt1 - toolPtProj);
toolUInterp = (1 - toolRInterp)*toolU1 + toolRInterp*toolUProj;

%% tool direction projection & interpolation
switch method
    case 1
        toolQuat1 = toolQuat(ind1,:);
        toolQuat2 = toolQuat(ind2,:);
        toolQuat3 = toolQuat(ind3,:);

        toolQuatProj = slerp(toolQuat2,toolQuat3,t);

        toolQuatInterp = slerp(toolQuat1,toolQuatProj,toolRInterp);

        varargout{1} = toolQuatInterp;

    case 2
        toolNorm1 = var1(:,ind1);
        toolNorm2 = var1(:,ind2);
        toolNorm3 = var1(:,ind3);

        % theta = t*vecAng(toolNorm2,toolNorm3,1);
        % axis = cross(toolNorm2,toolNorm3);
        % RInterp = expm(theta*axis);
        
        toolQuatProj = vecQuat(toolNorm2,toolNorm3);
        toolNormInterp = quat2rotm(toolQuatProj)*toolNorm2;
        % toolCutDirProj = quat2rotm(qInterp)*toolCutDir2;
        toolNormProj = vecOnPlane(toolNormInterp,toolPtProj,toolCutDir1);
        toolNormProj = toolNormProj./norm(toolNormProj);
        % 讲道理两个应该是相等的但是并不等，是否意味着我的刃口并不在对应平面内？
        % toolCutDirProj1 = vecRot(toolNorm2,toolNormProj)*quat2rotm(qInterp)*toolCutDir2; 
        toolCutDirProj = toolCutDir1;

        Rot12 = axesRot(toolNorm1,toolCutDir1,toolNormProj,toolCutDirProj,'zx');
        quat12 = rotm2quat(Rot12);
        toolQuatInterp = slerp([1,0,0,0],quat12,toolRInterp);
        toolNormInterp = quat2rotm(toolQuatInterp)*toolNorm1;
        toolNormInterp = toolNormInterp./norm(toolNormInterp);
        toolCutDirInterp = quat2rotm(toolQuatInterp)*toolCutDir1;
        toolCutDirInterp = toolCutDirInterp./norm(toolCutDirInterp);

        varargout{1} = toolNormInterp;
        varargout{2} = toolCutDirInterp;
    otherwise
        error('Invalid input. Not enough or tool many input parameters');
end

end