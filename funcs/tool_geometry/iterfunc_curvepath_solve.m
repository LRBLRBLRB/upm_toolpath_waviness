function [toolPathPt,toolQuat,toolContactU,surfPt,res,peakPt,uLim] = ...
    iterfunc_curvepath_solve(surfFuncr,surfFyr,toolData,toolPathPt, ...
    toolQuat,toolContactU,surfPt,delta0,aimRes,rRange,options)
%ITERFUNC_CURVEPATH the objective function of the iterative process to
%calculate the concentric toolpath of the aspheric surface, i.e., the
%toolpath of the curve, which is the generatrix of the surface.

% no need to use the global variables

arguments
    surfFuncr function_handle
    surfFyr function_handle
    toolData struct
    toolPathPt (3,1) double
    toolQuat (1,4) double
    toolContactU (1,1) double
    surfPt (3,1) double
    delta0 (1,1) double
    aimRes (1,1) double
    rRange (1,2) double
    options.radiusDisplay {mustBeMember(options.radiusDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'iter'
    options.maxIter (1,1) double = 50
    options.funcTol (1,1) double = 1e-6
    options.xTol (1,1) double = 1e-6
    options.useParallel logical = false
end

toolSp = toolData.toolBform; % B-form tool tip arc
if rRange(2) > rRange(1)
    toolDirect = [0;1;0];
else
    toolDirect = [0;-1;0];
end

uLim = [0;1]; % the interval of each toolpath
peakPt = zeros(3,1);
res = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes


    function diffRes = iterfunc(r)
        surfPt(:,ind) = [0;r;surfFuncr(r)];
        surfNorm = [0;surfFyr(r);-1];
        surfNorm = surfNorm./norm(surfNorm);
        % calculate the surfPt and toolpathPt from center to edge
        [toolPathPt(:,ind),toolQuat(ind,:),toolContactU(ind)] = ...
            curvetippos(toolData,surfPt(:,ind),surfNorm,[0;0;-1],toolDirect);
        % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

        % calculate the residual height of the loop and the inner nearest loop
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(ind,:))*toolSp1.coefs + toolPathPt(:,ind);
        toolContactPt1 = fnval(toolSp1,toolContactU(ind));
        toolSp2 = toolSp;
        toolSp2.coefs = quat2rotm(toolQuat(ind - 1,:))*toolSp2.coefs + toolPathPt(:,ind - 1);
        toolContactPt2 = fnval(toolSp2,toolContactU(ind - 1));
        [res(ind),peakPt(:,ind),uLim1,uLim2] = residual2D_numeric(toolSp1,toolSp2,1e-3, ...
            toolContactPt1,toolContactPt2,'DSearchn');
        if toolDirect(2) > 0
            uLim(2,ind) = uLim1;
            uLim(1,ind - 1) = uLim2;
        else
            uLim(1,ind) = uLim1;
            uLim(2,ind - 1) = uLim2;
        end
        diffRes = res(ind) - aimRes;
    end

opt1 = optimset('Display',options.radiusDisplay,'PlotFcns',{@optimplotx,@optimplotfval}, ...
    'MaxIter',options.maxIter,'TolFun',options.funcTol,'TolX',options.xTol);

opt2 = optimoptions('fsolve','Algorithm','levenberg-marquardt', ...
    'Display',options.radiusDisplay, ...
    'MaxIterations',options.maxIter,'FunctionTolerance',options.funcTol,'StepTolerance',options.xTol);
% 'PlotFcn',{@optimplotx,@optimplotfval},

r0 = rRange(1) + delta0;
r = rRange(1) + delta0;
ind = 2;
while (r - rRange(2))*delta0 < 0
    tic
    % iteration process
%     [r,diffRes1] = fzero(@iterfunc,r0,opt1);

    % [r,diffRes2] = fminsearch(@iterfunc,r0,opt1);

    r = fsolve(@iterfunc,r0,opt2);
    if uLim(1,ind-1) > uLim(2,ind-1)
       tmp = uLim(1,ind-1);
       uLim(1,ind-1) = uLim(2,ind-1);
       uLim(2,ind-1) = tmp;
    end

    scatter(toolPathPt(2,ind),toolPathPt(3,ind));
    toolSp0 = toolSp;
    toolSp0.coefs = quat2rotm(toolQuat(ind,:))*toolSp0.coefs + toolPathPt(:,ind);
    toolPt0 = fnval(toolSp0,uLim(1,ind):0.001:uLim(2,ind));
    toolContactPt0 = fnval(toolSp0,toolContactU(ind));
    plot(toolPt0(2,:),toolPt0(3,:));
    fprintf('No.%d\t toolpath point at r = %f is calculated within %fs.\n-----\n',ind,r,toc);

    ind = ind + 1;
    r0 = r + delta0;
end

if uLim(1,ind-1) > uLim(2,ind-1)
   tmp = uLim(1,ind-1);
   uLim(1,ind-1) = uLim(2,ind-1);
   uLim(2,ind-1) = tmp;
end

end