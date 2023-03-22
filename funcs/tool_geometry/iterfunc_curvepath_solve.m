function [toolPathPt,toolQuat,toolContactU,surfPt,res,peakPt,uLim] = ...
    iterfunc_curvepath_solve(surfFuncr,surfFxr,toolData,toolPathPt, ...
    toolQuat,toolContactU,surfPt,delta0,aimRes,rRange,options)
%ITERFUNC_CURVEPATH the objective function of the iterative process to
%calculate the concentric toolpath of the aspheric surface, i.e., the
%toolpath of the curve, which is the generatrix of the surface.

% no need to use the global variables

arguments
    surfFuncr function_handle
    surfFxr function_handle
    toolData struct
    toolPathPt (3,1) double
    toolQuat (1,4) double
    toolContactU (1,1) double
    surfPt (3,1) double
    delta0 (1,1) double
    aimRes (1,1) double
    rRange (1,2) double
    options.algorithm {mustBeMember(options.algorithm, ...
        {'trust-region-dogleg','trust-region','levenberg-marquardt', ...
        'fzero','search-bisection','genetic','particle-swarm'})} = 'search-bisection'
    options.directionType {mustBeMember(options.directionType, ...
        {'quaternion','norm-cut','norm-feed'})} = 'norm-cut'
    options.radiusDisplay {mustBeMember(options.radiusDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'iter'
    options.optimopt
end

toolSp = toolData.toolBform; % B-form tool tip arc
if rRange(2) > rRange(1) % to ensure the rake face on top
    % if range(2) > rRange(1), then feed direction is center-to-edge, and
    % the toolDirect is [0;1;0]. So cut direction is [-1;0;0]
    cutDirect = [0;1;0];
else
    cutDirect = [0;-1;0];
end

uLim = [0;1]; % the interval of each toolpath
peakPt = zeros(3,1);
res = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes


    function diffRes = iterfunc(r)
        surfPt(:,ind) = [r;0;surfFuncr(r)];
        surfNorm = [surfFxr(r);0;-1];
        surfNorm = surfNorm./norm(surfNorm);
        % calculate the surfPt and toolpathPt from center to edge
        [toolPathPt(:,ind),toolQuat(ind,:),toolContactU(ind)] = ...
            curvetippos(toolData,surfPt(:,ind),surfNorm,[0;0;-1],cutDirect, ...
            "directionType",options.directionType);
        % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

        % calculate the residual height of the loop and the inner nearest loop
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(ind,:))*toolSp1.coefs + toolPathPt(:,ind);
        toolContactPt1 = fnval(toolSp1,toolContactU(ind));
        toolSp2 = toolSp;
        toolSp2.coefs = quat2rotm(toolQuat(ind - 1,:))*toolSp2.coefs + toolPathPt(:,ind - 1);
        toolContactPt2 = fnval(toolSp2,toolContactU(ind - 1));
        [res(ind),peakPt(:,ind),uLim1,uLim2] = residual2D_numeric(toolSp1,toolSp2,1e-4, ...
            toolContactPt1,toolContactPt2,'DSearchn');
        if isinf(res(ind))
            fprintf('No intersection of the current tool path.\n');
            return;
        end
        if cutDirect(2) > 0
            uLim(2,ind) = uLim1;
            uLim(1,ind - 1) = uLim2;
        else
            uLim(1,ind) = uLim1;
            uLim(2,ind - 1) = uLim2;
        end
        diffRes = res(ind) - aimRes;

        scatter(toolPathPt(1,ind),toolPathPt(3,ind),36,[0.4940,0.1840,0.5560]);
        scatter(toolPathPt(1,ind - 1),toolPathPt(3,ind - 1),36,[0.4940,0.1840,0.5560]);
        toolPt1 = fnval(toolSp1,0:0.001:1);
        plot(toolPt1(1,:),toolPt1(3,:),'Color',[0.7,.7,.70]);
        toolPt2 = fnval(toolSp2,0:0.001:1);
        plot(toolPt2(1,:),toolPt2(3,:),'Color',[0.7,.7,.70]);
        scatter(toolContactPt1(1),toolContactPt1(3),18,[0.929,0.694,0.1250],"filled");
        scatter(toolContactPt2(1),toolContactPt2(3),18,[0.929,0.694,0.1250],"filled");
        scatter(peakPt(1,ind),peakPt(3,ind),18,[0.850,0.325,0.0980],"filled");
    end

r = rRange(1);
ind = 1;
while (r - rRange(2))*delta0 < 0
    ind = ind + 1;
    tic
    % initial value
    delta = 2*sqrt((2*toolData.radius*aimRes - aimRes^2))*cos(atan(surfFxr(r)));
    % direction to iterate
    delta = sign(delta0)*delta;
    r0 = r + delta;

    % iteration process
    switch options.algorithm
        case {'trust-region-dogleg','trust-region','levenberg-marquardt'}
            options.optimopt = optimoptions(options.optimopt,'Algorithm',options.algorithm);
            r = fsolve(@iterfunc,r0,options.optimopt);
        case 'fzero'
            [r,diffRes1] = fzero(@iterfunc,r0,options.optimopt);
        case 'search-bisection'
            h = delta/50;
            r = mysearch(@iterfunc,r0,h,[r0 + delta,r],options.xTol);
        case 'genetic'
            iniMat
            options.optimopt = optimoptions('ga','UseParallel',true, ...
                'InitialPopulationMatrix',iniMat);
            r = ga(@iterfunc,1,[],[],[],[],r0 + delta,r,options.optimopt);
        case 'particle-swarm'
            r = particleswarm(@iterfunc,1,r0 + delta,r,options.optimopt);
    end
    % [r,diffRes2] = fminsearch(@iterfunc,r0,opt1);

%     if isempty(r) || isnan(r)
%         r = r0;
%         iterfunc(r0);
%     end

%     if uLim(1,ind-1) > uLim(2,ind-1)
%        tmp = uLim(1,ind-1);
%        uLim(1,ind-1) = uLim(2,ind-1);
%        uLim(2,ind-1) = tmp;
%     end

%     scatter(toolPathPt(1,ind),toolPathPt(3,ind),36,[0.4940,0.1840,0.5560]);
%     toolSp0 = toolSp;
%     toolSp0.coefs = quat2rotm(toolQuat(ind,:))*toolSp0.coefs + toolPathPt(:,ind);
%     toolPt0 = fnval(toolSp0,0:0.001:1);
%     toolContactPt0 = fnval(toolSp0,toolContactU(ind));
%     plot(toolPt0(1,:),toolPt0(3,:),'Color',[0.7,.7,.70]);
%     scatter(toolContactPt0(1),toolContactPt0(3),18,[0.929,0.694,0.1250],"filled");
%     scatter(peakPt(1,ind),peakPt(3,ind),18,[0.850,0.325,0.0980],"filled");

    fprintf('-----\nNo.%d\t toolpath point at [r = %f] is calculated within %fs.\n-----\n',ind,r,toc);
    % pause(0.5);
end

% if uLim(1,ind-1) > uLim(2,ind-1)
%    tmp = uLim(1,ind-1);
%    uLim(1,ind-1) = uLim(2,ind-1);
%    uLim(2,ind-1) = tmp;
% end

end