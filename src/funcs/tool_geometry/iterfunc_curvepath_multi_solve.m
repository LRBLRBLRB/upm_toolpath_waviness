function [curvePathPt,curveQuat,curveContactU,curvePt,curveRes, ...
    curvePeakPt,curveInterPt,curveULim] = iterfunc_curvepath_multi_solve( ...
    surfFunc,surfFx,toolData,curvePathPt,curveQuat,curveContactU,curvePt, ...
    delta0,aimRes,rRange,options)
%ITERFUNC_CURVEPATH_MULTI_SOLVE the objective function of the iterative 
%process to calculate the concentric toolpath of the aspheric surface, 
%i.e., the toolpath of the curve, which is the generatrix of the surface.

% no need to use the global variables

arguments
    surfFunc function_handle
    surfFx function_handle
    toolData struct
    curvePathPt (3,1) double
    curveQuat (1,4) double
    curveContactU (1,1) double
    curvePt (3,1) double
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

% the rest
curveULim = {[0;1]}; % the interval of each toolpath
curvePeakPt = zeros(5,1);
curveInterPt = {zeros(3,1)};
curveRes = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes


    function diffRes = iterfunc(r)
        curveULim{ind} = [0;1];
        curvePt(:,ind) = [r;0;surfFunc(r)];
        if abs(curvePt(1,ind) -  curvePt(1,ind - 1)) < 1e-3
            diffRes = - aimRes;
            return;
        end
        surfNorm = [surfFx(r);0;-1];
        surfNorm = surfNorm./norm(surfNorm);
        % calculate the surfPt and toolpathPt from center to edge
        [curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind)] = ...
            curvetippos(toolData,curvePt(:,ind),surfNorm,[0;0;-1],cutDirect, ...
            "directionType",'norm-cut');
        % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;
    
        % calculate the residual height of the loop and the inner nearest loop
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(curveQuat(ind,:))*toolSp1.coefs + curvePathPt(:,ind);
        % toolContactPt1 = fnval(toolSp1,curveContactU(ind));
        toolSp2 = toolSp;
        toolSp2.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp2.coefs + curvePathPt(:,ind - 1);
        % toolContactPt2 = fnval(toolSp2,curveContactU(ind - 1));
    
        [curveRes(ind),curvePeakPt(:,ind),curveInterPt{ind},curveULim1, ...
            curveULim2] = residual2D_multi(toolSp1,toolSp2,1e-5, ...
            curvePt(:,ind),curvePt(:,ind - 1),curveULim{ind - 1});
        % curvePeakPt(5,ind) = curvePeakPt(5,ind) + ind;
    
        if isinf(curveRes(ind))
            fprintf('No intersection of the current tool path.\n');
            return;
        end

        diffRes = curveRes(ind) - aimRes;

%         scatter(curvePathPt(1,ind),curvePathPt(3,ind),36,[0.4940,0.1840,0.5560]);
%         scatter(curvePathPt(1,ind - 1),curvePathPt(3,ind - 1),36,[0.4940,0.1840,0.5560]);
%         toolPt1 = fnval(toolSp1,0:0.001:1);
%         plot(toolPt1(1,:),toolPt1(3,:),'Color',[0.7,.7,.70]);
%         toolPt2 = fnval(toolSp2,0:0.001:1);
%         plot(toolPt2(1,:),toolPt2(3,:),'Color',[0.7,.7,.70]);
%         scatter(toolContactPt1(1),toolContactPt1(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(toolContactPt2(1),toolContactPt2(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(curvePeakPt(1,ind),curvePeakPt(3,ind),18,[0.850,0.325,0.0980],"filled");
    end


r = rRange(1);
ind = 1;
while (r - rRange(2))*delta0 < 0
    curveULim1 = [];
    curveULim2 = [];
    ind = ind + 1;
    tic
    % initial value
    delta = 2*sqrt((2*toolData.radius*aimRes - aimRes^2))*cos(atan(surfFx(r)));
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
            r = mysearch(@iterfunc,r0,h,[r0 + delta,r],options.optimopt.XTol);
        case 'traverse-bisection'
        case 'genetic'
            % 'InitialPopulationMatrix',iniMat
%             options.optimopt = optimoptions(options.optimopt, ...
%                 'PlotFcn',{'gaplotscorediversity','gaplotselection','gaplotgenealogy'});
            r = ga(@(x)abs(iterfunc(x)),1,[],[],[],[],r0,r - aimRes,[],options.optimopt);
        case 'particle-swarm'
            options.optimopt = optimoptions(options.optimopt,'PlotFcn','pswplotbestf');
            r = particleswarm(@(x)abs(iterfunc(x)),1,r0 + delta,r - aimRes,options.optimopt);
    end
    % [r,diffRes2] = fminsearch(@iterfunc,r0,opt1);

    curveULim{ind} = curveULim1;
    curveULim{ind - 1} = curveULim2;

%         scatter(curvePathPt(1,ind),curvePathPt(3,ind),36,[0.4940,0.1840,0.5560]);
%         scatter(curvePathPt(1,ind - 1),curvePathPt(3,ind - 1),36,[0.4940,0.1840,0.5560]);
%         toolSp10 = toolSp;
%         toolSp10.coefs = quat2rotm(curveQuat(ind,:))*toolSp10.coefs + curvePathPt(:,ind);
%         toolContactPt10 = fnval(toolSp10,curveContactU(ind));
%         toolSp20 = toolSp;
%         toolSp20.coefs = quat2rotm(curveQuat(ind - 1,:))*toolSp20.coefs + curvePathPt(:,ind - 1);
%         toolContactPt20 = fnval(toolSp20,curveContactU(ind - 1));
%         toolPt1 = fnval(toolSp10,0:0.001:1);
%         plot(toolPt1(1,:),toolPt1(3,:),'Color',[0.7,.7,.70]);
%         toolPt2 = fnval(toolSp20,0:0.001:1);
%         plot(toolPt2(1,:),toolPt2(3,:),'Color',[0.7,.7,.70]);
%         scatter(toolContactPt10(1),toolContactPt10(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(toolContactPt20(1),toolContactPt20(3),18,[0.929,0.694,0.1250],"filled");
%         scatter(curveInterPt{ind}(1,:),curveInterPt{ind}(3,:),18,[0.850,0.325,0.0980],"filled");

    fprintf('No.%d\t toolpath %f\t[r = %f] is calculated within %fs.\n-----\n',ind,curvePathPt(1,ind),r,toc);
end

end