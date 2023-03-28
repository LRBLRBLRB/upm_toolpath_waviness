function [toolPathPt,toolQuat,toolContactU,surfPt,res,peakPt,uLim] = ...
    iterfunc_curvepath_iter(surfFuncr,surfFx,toolData, ...
    toolPathPt,toolQuat,toolContactU,surfPt,delta0,aimRes,rMax,options)
%ITERFUNC_CURVEPATH the objective function of the iterative process to
%calculate the concentric toolpath of the aspheric surface, i.e., the
%toolpath of the curve, which is the generatrix of the surface.

% no need to use the global variables

arguments
    surfFuncr function_handle
    surfFx function_handle
    toolData struct
    toolPathPt (3,1) double
    toolQuat (1,4) double
    toolContactU (1,1) double
    surfPt (3,1) double
    delta0 (1,1) double
    aimRes (1,1) double
    rMax (1,1) double
    options.radiusDisplay {mustBeMember(options.radiusDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'iter'
    options.maxIter (1,1) double = 50
    options.funcTol (1,1) double = 1e-6
    options.xTol (1,1) double = 1e-6
    options.useParallel logical = false
    options.iniDisplay {mustBeMember(options.iniDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
    options.finalDisplay {mustBeMember(options.finalDisplay,{'none','off', ...
        'iter','iter-detailed','final','final-detailed'})} = 'none'
    options.finalFuncTol {mustBePositive} = 1e-3
    options.finalStepTol {mustBePositive} = 1e-3
    options.uQTol {mustBePositive} = 1e-3
end

toolSp = toolData.toolBform; % B-form tool tip arc

uLim = [0;1]; % the interval of each toolpath
peakPt = zeros(3,1);
res = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes


    function diffRes = iterfunc(r)
        toolPathPt(:,ind) = [r;0;surfFuncr(r)];
        % calculate the surfPt and toolpathPt from center to edge
        [toolPathPt(:,ind),toolQuat(ind,:),toolContactU(ind),surfPt(:,ind)] = curvepos( ...
            surfFuncr,surfFx,toolData,toolPathPt(:,ind),[0;0;-1], ...
            'iniDisplay',options.iniDisplay,'finalDisplay',options.finalDisplay, ...
            'finalFuncTol',options.finalFuncTol,'finalStepTol',options.finalStepTol, ...
            'uQTol',options.uQTol,'useParallel',options.useParallel);
        % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

        % calculate the residual height of the loop and the inner nearest loop
        toolSp1 = toolSp;
        toolSp1.coefs = quat2rotm(toolQuat(ind,:))*toolSp1.coefs + toolPathPt(:,ind);
        toolContactPt1 = fnval(toolSp1,toolContactU(ind));
        toolSp2 = toolSp;
        toolSp2.coefs = quat2rotm(toolQuat(ind - 1,:))*toolSp2.coefs + toolPathPt(:,ind - 1);
        toolContactPt2 = fnval(toolSp2,toolContactU(ind - 1));
        [res(ind),peakPt(:,ind),uLim(1,ind),uLim(2,ind - 1)] = residual2D_numeric(toolSp1,toolSp2,1e-3, ...
            toolContactPt1,toolContactPt2,'DSearchn');

        diffRes = res(ind) - aimRes;
    end

opt1 = optimset('Display',options.radiusDisplay,'PlotFcns',{@optimplotx,@optimplotfval}, ...
    'MaxIter',options.maxIter,'TolFun',options.funcTol,'TolX',options.xTol);

opt2 = optimoptions('fsolve','Algorithm','levenberg-marquardt', ...
    'Display',options.radiusDisplay,'PlotFcn',{@optimplotx,@optimplotfval}, ...
    'MaxIterations',options.maxIter,'FunctionTolerance',options.funcTol,'StepTolerance',options.xTol);

r0 = delta0;
r = delta0;
ind = 2;
while r <= rMax
    tic
    % iteration process
%     [r,diffRes1] = fzero(@iterfunc,r0,opt1);

    % [r,diffRes2] = fminsearch(@iterfunc,r0,opt1);

    [r,~] = fsolve(@iterfunc,r0,opt2);

    fprintf('No.%d\t toolpath point at r = %f is calculated within %fs.\n-----\n',ind,r,toc);

    ind = ind + 1;
    r0 = r + delta0;
end

end