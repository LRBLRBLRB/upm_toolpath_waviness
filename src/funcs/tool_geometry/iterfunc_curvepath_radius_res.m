function [curvePathPt,curveQuat,curveContactU,curvePt,curveRes, ...
    curvePeakPt,curveULim] = iterfunc_curvepath_radius_res( ...
    surfFunc,surfFx,toolData,curvePathPt,curveQuat,curveContactU,curvePt, ...
    delta0,aimRes,rRange,options)
%ITERFUNC_CURVEPATH the objective function of the iterative process to
%calculate the concentric toolpath of the aspheric surface, i.e., the
%toolpath of the curve, which is the generatrix of the surface.

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

if abs(rRange(2)) > abs(rRange(1)) % to ensure the rake face on top
    % if range(2) > rRange(1), then feed direction is center-to-edge, and
    % the toolDirect is [0;1;0]. So cut direction is [-1;0;0]
    cutDirect = [0;1;0];
else
    cutDirect = [0;-1;0];
end

curveRes = 5*aimRes; % the residual height, initialized with 5 times the standard aimRes
curvePeakPt = zeros(3,1); % peak point coordinates
curveULim = [0;2*pi];

    function diffRes = iterfunc(r)
        curveULim(:,ind) = [0;2*pi];
        curvePathPt(:,ind) = [r;0;surfFunc(r)];
        if abs(curvePathPt(1,ind) -  curvePathPt(1,ind - 1)) < 1e-3
            diffRes = - aimRes;
            return;
        end
        % calculate the surfPt and toolpathPt from center to edge
        [curvePathPt(:,ind),curveQuat(ind,:),curveContactU(ind),curvePt(:,ind)] = ...
            radiuspos(surfFunc,surfFx,toolData,curvePathPt(:,ind),[0;0;-1],cutDirect, ...
            'directionType',options.directionType);
        % toolNormDirect(:,ind) = quat2rotm(toolQuat(ind,:))*toolData.toolEdgeNorm;

        % calculate the residual height of the loop and the inner nearest loop
        a = norm(curvePathPt(:,ind) - curvePathPt(:,ind - 1));
        curveRes(ind) = toolData.toolRadius - sqrt(toolData.toolRadius^2 - a^2/4);
        diffRes = curveRes(ind) - aimRes;

%         figure;
%         scatter(toolPathPt(1,ind),toolPathPt(3,ind));
%         rectangle('Position',[toolPathPt(1,ind)-toolData.toolRadius,toolPathPt(3,ind)-toolData.toolRadius, ...
%             2*toolData.toolRadius,2*toolData.toolRadius],'Curvature',[1,1]);
%         hold on;
%         scatter(toolPathPt(1,ind - 1),toolPathPt(3,ind - 1));
%         rectangle('Position',[toolPathPt(1,ind - 1)-toolData.toolRadius,toolPathPt(3,ind - 1)-toolData.toolRadius, ...
%             2*toolData.toolRadius,2*toolData.toolRadius],'Curvature',[1,1]);

%         scatter(toolPathPt(1,ind),toolPathPt(3,ind),36,[0.4940,0.1840,0.5560]);
%         scatter(toolPathPt(1,ind - 1),toolPathPt(3,ind - 1),36,[0.4940,0.1840,0.5560]);
%         toolThe = 0:0.01:2*pi;
%         toolPt1(1,:) = toolPathPt(1,ind) + toolData.toolRadius*cos(toolThe);
%         toolPt1(3,:) = toolPathPt(3,ind) + toolData.toolRadius*sin(toolThe);
%         plot(toolPt1(1,:),toolPt1(3,:),'Color',[0.7,.7,.7]);
%         toolPt2(1,:) = toolPathPt(1,ind - 1) + toolData.toolRadius*cos(toolThe);
%         toolPt2(3,:) = toolPathPt(3,ind - 1) + toolData.toolRadius*sin(toolThe);
%         plot(toolPt2(1,:),toolPt2(3,:),'Color',[0.7,.7,.7]);
    end

r = curvePathPt(1,1);
ind = 1;
while (r - rRange(2))*delta0 < 0
    ind = ind + 1;
    tic
    % initial value
    delta = 2*sqrt((2*toolData.toolRadius*aimRes - aimRes^2))*cos(atan(surfFx(r)));
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
            if delta > 0
                searchInt = [r,r0 + delta];
            else
                searchInt = [r0 + delta,r];
            end
            r = mysearch(@iterfunc,r0,h,searchInt,options.optimopt.XTol);
        case 'genetic'
            % 'InitialPopulationMatrix',iniMat
%             options.optimopt = optimoptions('ga','UseParallel',true, ...
%                 'PlotFcn',{'gaplotscorediversity','gaplotselection','gaplotgenealogy'}, ...
%                 'CreationFcn','gacreationuniform', ...
%                 'SelectionFcn','selectionstochunif', ...
%                 'CrossoverFcn','crossoversinglepoint', ...
%                 'MutationFcn',{@mutationgaussian}, ...
%                 'HybridFcn','fminunc', ...
%                 'Display','diagnose');

            options.optimopt = optimoptions(options.optimopt, ...
                'PlotFcn',{'gaplotscorediversity','gaplotselection','gaplotscores'});

%             defaultopt = struct('PopulationType', 'doubleVector', ...
%                 'PopInitRange', [], ... 
%                 'PopulationSize', '50 when numberOfVariables <= 5, else 200', ... 
%                 'EliteCount', '0.05*PopulationSize', ...  
%                 'CrossoverFraction', 0.8, ...
%                 'MigrationDirection','forward', ...
%                 'MigrationInterval',20, ...
%                 'MigrationFraction',0.2, ...
%                 'Generations', '100*numberOfVariables', ...
%                 'FitnessLimit', -inf, ...
%                 'StallTest', 'averageChange', ... 
%                 'StallGenLimit', 50, ...
%                 'TolFun', 1e-6, ...
%                 'TolCon', 1e-3, ...
%                 'FitnessScalingFcn', @fitscalingrank);

            r = ga(@(x)abs(iterfunc(x)),1,[],[],[],[],r0 + delta,r,[],options.optimopt);
        case 'particle-swarm'
            options.optimopt = optimoptions(options.optimopt,'PlotFcn','pswplotbestf');
            r = particleswarm(@(x)abs(iterfunc(x)),1,r0 + delta,r,options.optimopt);
    end
    % [r,diffRes2] = fminsearch(@iterfunc,r0,opt1);

    % 
    Norm = norm(curvePathPt(:,ind) - curvePathPt(:,ind - 1));
    curvePeakPt(:,ind) = 1/2*(curvePathPt(:,ind) + curvePathPt(:,ind - 1)) + sqrt(toolData.toolRadius^2 - Norm^2/4) ...
        *roty(-90)*(curvePathPt(:,ind) - curvePathPt(:,ind - 1))/Norm;    
    vec1 = curvePathPt(:,ind) - curvePeakPt;
    curveULim(1,ind) = atan2(vec1(end),vec1(1));
    vec2 = curvePathPt(:,ind - 1) - curvePeakPt;
    curveULim(2,ind - 1) = atan2(vec2(end),vec2(1));

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
%     toolThe0 = 0:0.01:2*pi;
%     toolPt0(1,:) = toolPathPt(1,ind) + toolData.toolRadius*cos(toolThe0);
%     toolPt0(3,:) = toolPathPt(3,ind) + toolData.toolRadius*sin(toolThe0);
%     plot(toolPt0(1,:),toolPt0(3,:),'Color',[0.7,.7,.7]);
%     plot(peakPt(1),peakPt(3),'.','MarkerSize',12);

    fprintf('-----\nNo.%d\t toolpath point at [r = %f] is calculated within %fs.\n-----\n',ind,r,toc);
    % pause(0.5);
end

% center point correction
curvePathPt(1,end) = 0;
[curvePathPt(:,end),curveQuat(end,:),curveContactU(end),curvePt(:,end)] = ...
    radiuspos(surfFunc,surfFx,toolData,curvePathPt(:,end),[0;0;-1],cutDirect, ...
    'directionType',options.directionType);
a = norm(curvePathPt(:,end) - curvePathPt(:,end - 1));
curveRes(end) = toolData.toolRadius - sqrt(toolData.toolRadius^2 - a^2/4);
end