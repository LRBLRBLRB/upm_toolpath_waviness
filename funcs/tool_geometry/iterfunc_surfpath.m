function iterfunc_surfpath(r0,toolData,delta0,surfFunc,surfNormFunc, ...
    angularDiscrete,conThetaBound,maxAngPtDist,arcLength,angularLength, ...
    aimRes,rMax,tRes0)
%     loopPtNum,accumPtNum,toolREach,surfPt,surfNorm,surfDirect, ...
%     toolQuat,toolVec,toolContactU,isCollision, ...
%     toolPathPt,toolCutDirect,toolNormDirect,uLim,peakPt,res)
%CIRCLETOOLPATH 

global loopPtNum accumPtNum toolREach toolRAccum toolNAccum surfPt surfNorm surfDirect ...
    toolQuat toolVec toolContactU isCollision toolPathAngle ...
    toolPathPt toolCutDirect toolNormDirect uLim peakPt res;

toolSp = toolData.toolBform; % B-form tool tip arc
toolRadius = toolData.radius; % fitted tool radius
% rStep = delta;
% iter = 1;

    function diffRes = iterfunc(r)
%     while true
        % calculate the discretization scatters
        switch angularDiscrete
            case 'Constant Arc'
                % if r*maxAngPtDist < arcLength, then discrete the circle with constant angle
                conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
                    ceil(2*pi/min(maxAngPtDist,arcLength/r)) + 1);
                conTheta(end) = [];
            case 'Constant Angle'
                conTheta = linspace(conThetaBound(1),conThetaBound(2), ...
                    ceil(2*pi/angularLength) + 1);
                conTheta(end) = [];
            otherwise
                error('Invalid angular discretization type.');
        end

        clear surfPtTmp surfNormTmp surfDirectTmp toolQuatTmp toolVecTmp ...
            toolContactUTmp isCollisionTmp toolPathPtTmp toolNormDirectTmp ...
            toolCutDirectTmp;
        clear uLimTmp peakPtOutTmp resTmp loopPtNumLast;
    
        % calculate the cutting pts
        loopPtNumTmp = length(conTheta); % the number of points in the loop
        loopPtNumLast = loopPtNum(end);  % the number of points in the former loop
        surfPtTmp(1,:) = r*cos(conTheta); % x coordinates of surface cutting points in the loop
        surfPtTmp(2,:) = r*sin(conTheta); % y coordinates of surface cutting points
        surfPtTmp(3,:) = surfFunc(surfPtTmp(1,:),surfPtTmp(2,:)); % z coordinates of surface cutting points
        surfNormTmp = surfNormFunc(surfPtTmp(1,:),surfPtTmp(2,:)); % normolization of the normal vector
        surfDirectTmp = cutdirection(surfPtTmp,[0;0;0],conTheta,'method','vertical'); % cutting direction of each points
        
        % calculate the tool path and residual height
        toolQuatTmp = zeros(loopPtNumTmp,4);
        toolVecTmp = zeros(3,loopPtNumTmp);
        toolContactUTmp = zeros(1,loopPtNumTmp);
        isCollisionTmp = zeros(1,loopPtNumTmp);
        toolPathPtTmp = zeros(3,loopPtNumTmp);
        toolCutDirectTmp = zeros(3,loopPtNumTmp);
        toolNormDirectTmp = zeros(3,loopPtNumTmp);
        for indTmp = 1:loopPtNumTmp
            [toolQuatTmp(indTmp,:),toolVecTmp(:,indTmp),toolContactUTmp(indTmp)] = tooltippos( ...
                toolData,surfPtTmp(:,indTmp),surfNormTmp(:,indTmp),[0;0;-1],surfDirectTmp(:,indTmp));
            toolPathPtTmp(:,indTmp) = quat2rotm(toolQuatTmp(indTmp,:))*toolData.center + toolVecTmp(:,indTmp);
            toolCutDirectTmp(:,indTmp) = quat2rotm(toolQuatTmp(indTmp,:))*toolData.cutDirect;
            toolNormDirectTmp(:,indTmp) = quat2rotm(toolQuatTmp(indTmp,:))*toolData.toolEdgeNorm;
    
            % debug
            % plot3(toolPathPtTmp(1,indTmp),toolPathPtTmp(2,indTmp),toolPathPtTmp(3,indTmp), ...
            %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]); hold on;
            % plot3(surfPtTmp(1,indTmp),surfPtTmp(2,indTmp),surfPtTmp(3,indTmp), ...
            %     '.','MarkerSize',36,'Color',[0,0.4470,0.7410]); hold on;
            % qui1 = quiver3(toolPathPtTmp(1,indTmp),toolPathPtTmp(2,indTmp),toolPathPtTmp(3,indTmp), ...
            %     toolNormDirectTmp(1,indTmp),toolNormDirectTmp(2,indTmp),toolNormDirectTmp(3,indTmp), ...
            %     300,'r');
            % qui2 = quiver3(toolPathPtTmp(1,indTmp),toolPathPtTmp(2,indTmp),toolPathPtTmp(3,indTmp), ...
            %     toolCutDirectTmp(1,indTmp),toolCutDirectTmp(2,indTmp),toolCutDirectTmp(3,indTmp), ...
            %     300,'b');
            % toolSp1 = toolData.toolBform;
            % toolSp1.coefs = quat2rotm(toolQuatTmp(indTmp,:))*toolData.toolBform.coefs + toolVecTmp(:,indTmp);
            % Q = fnval(toolSp1,0:0.01:1);
            % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5);
            % delete(qui1);
            % delete(qui2);
        end
    
        % restore the data that will be used in the residual height calculation
        toolPathPtRes = [toolPathPt(:,end - loopPtNumLast + 1:end),toolPathPtTmp];
        toolNormDirectRes = [toolNormDirect(:,end + 1 - loopPtNumLast:end),toolNormDirectTmp];
        toolCutDirectRes = [toolCutDirect(:,end + 1 - loopPtNumLast:end),toolCutDirectTmp];
        toolContactURes = [toolContactU(end + 1 - loopPtNumLast:end),toolContactUTmp];
        uLimTmp = [uLim(:,end - loopPtNumLast + 1:end), ...
            [zeros(1,loopPtNumTmp);ones(1,loopPtNumTmp)]]; % the interval of each toolpath
        resTmp = [res(:,end - loopPtNumLast + 1:end), ...
            5*aimRes*ones(2,loopPtNumTmp)];
        peakPtInTmp = [peakPt(1:3,end - loopPtNumLast + 1:end), ...
            zeros(3,loopPtNumTmp)];
        peakPtOutTmp = zeros(3,loopPtNumLast + loopPtNumTmp);
    
        % calculate the residual height of the loop and the inner nearest loop
        angle = atan2(toolPathPtRes(2,:),toolPathPtRes(1,:));
        % inner side of each point on the tool path
        angleIn = angle(1:loopPtNumLast);
        for indIn1 = loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp
            angleDel = angleIn - angle(indIn1);
            [indIn2,indIn3] = getInnerLoopToolPathIndex(angleIn,angleDel);
            % if isempty(angleN(angleDel >= 0))
            %     % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
            %     angleDel = angleDel + 2*pi;
            % end
            % indIn2 = find(angleN == min(angleN(angleDel >= 0)));
            % if isempty(angleN(angleDel < 0))
            %     angleDel = angleDel - 2*pi;
            % end
            % indIn3 = find(angleN == max(angleN(angleDel < 0)));
            [resTmp(1,indIn1),peakPtInTmp(:,indIn1),uLimTmp(:,indIn1)] = residual3D( ...
                toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
                toolData,toolRadius,uLimTmp(:,indIn1),indIn1,indIn2,indIn3);
    
            % debug
            % plot3(toolPathPtRes(1,ii),toolPathPtRes(2,ii),toolPathPtRes(3,ii), ...
            %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
            % toolSp1 = toolSp;
            % R1 = axesRot([0;0;1],[1;0;0],toolNormDirectRes(:,ii),toolCutDirectRes(:,ii),'zx');
            % toolSp1.coefs = R1*toolSp.coefs + toolPathPtRes(:,ii);
            % Q = fnval(toolSp1,uLimTmp(1,ii):0.01:uLimTmp(2,ii));
            % plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',1);
        end

        % if residual height does not satisfy the reqiurement
        maxRes = max(resTmp(1,loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp),[],"all");
        diffRes = maxRes - aimRes;
        if maxRes == Inf
            error('The radius %d is too large to fit the surface',r);
        end
%         if maxRes < aimRes
%             break;
%         else
%             fprintf('\tIter %d: The maximum residual height is %f %cm.\n',iter,maxRes,char([956]));
%             delta = delta/3;
%             r = r - delta;
%             iter = iter + 1;
%         end
    end

opt1 = optimset('Display','iter','MaxIter',50,'PlotFcns',{@optimplotx,@optimplotfval}, ...
    'TolFun',1,'TolX',1e-3);

opt2 = optimoptions('fsolve','Algorithm','levenberg-marquardt', ...
    'Display','iter-detailed','PlotFcn',{@optimplotx,@optimplotfval}, ...
    'MaxIterations',100,'FunctionTolerance',1e-6,'StepTolerance',1e-6);

while true
    tic
    angle = [];
    loopPtNumLast = [];
    loopPtNumTmp = [];
    toolPathPtRes = [];
    toolNormDirectRes = [];
    toolCutDirectRes = [];
    toolContactURes = [];
    surfPtTmp = [];
    surfNormTmp = [];
    surfDirectTmp = [];
    toolQuatTmp = [];
    toolVecTmp = [];
    toolContactUTmp = [];
    isCollisionTmp = [];
    toolPathPtTmp = [];
    toolNormDirectTmp = [];
    toolCutDirectTmp = [];
    peakPtInTmp = [];

%     [r1,diffRes1] = fzero(@iterfunc,r0,opt1);

    % [r2,diffRes2] = fminsearch(@iterfunc,r0,opt1);

    [r1,~] = fsolve(@iterfunc,r0,opt2);


    % outer side of each point in the tool path
    angleOut = angle(loopPtNumLast + 1:loopPtNumLast + loopPtNumTmp);
    parfor indOut1 = 1:loopPtNumLast
        angleDel = angleOut - angle(indOut1);
        [indOut2,indOut3] = getOuterLoopToolPathIndex(angleOut,angleDel,loopPtNumLast);
        % if isempty(angleN(angleDel >= 0))
        %     % to avoid that angle(indOut1) is cloesd to -pi, and smaller than each elements
        %     angleDel = angleDel + 2*pi;
        % end
        % indOut2 = loopPtNumLast + find(angleN == min(angleN(angleDel >= 0)));
        % if isempty(angleN(angleDel < 0))
        %     angleDel = angleDel - 2*pi;
        % end
        % indOut3 = loopPtNumLast + find(angleN == max(angleN(angleDel < 0)));
        [resTmp(2,indOut1),peakPtOutTmp(:,indOut1),uLimTmp(:,indOut1)] = residual3D( ...
            toolPathPtRes,toolNormDirectRes,toolCutDirectRes,toolContactURes, ...
            toolData,toolRadius,uLimTmp(:,indOut1),indOut1,indOut2,indOut3);
    end
    
    % then store the data of this loop
    loopPtNum = [loopPtNum,loopPtNumTmp];
    accumPtNum = [accumPtNum,accumPtNum(end) + loopPtNumTmp];
    toolREach = [toolREach,r1];
    toolRAccum = [toolRAccum,r1*ones(1,loopPtNumTmp)];
    toolPathAngle = [toolPathAngle,wrapTo2Pi(atan2(toolPathPtTmp(2,:),toolPathPtTmp(1,:)))];
    toolNAccum = [toolNAccum,(length(loopPtNum) - 1)*ones(1,loopPtNumTmp)];
    if toolPathAngle(accumPtNum(end - 1) + 1) > 6
        toolPathAngle(accumPtNum(end - 1) + 1) = toolPathAngle(accumPtNum(end - 1) + 1) - 2*pi;
    end
    surfPt = [surfPt,surfPtTmp];
    surfNorm = [surfNorm,surfNormTmp];
    surfDirect = [surfDirect,surfDirectTmp];
    toolQuat = [toolQuat;toolQuatTmp];
    toolVec = [toolVec,toolVecTmp];
    toolContactU = [toolContactU,toolContactUTmp];
    isCollision = [isCollision,isCollisionTmp];
    toolPathPt = [toolPathPt,toolPathPtTmp];
    toolNormDirect = [toolNormDirect,toolNormDirectTmp];
    toolCutDirect = [toolCutDirect,toolCutDirectTmp];

    uLim(:,end - loopPtNumLast + 1:end) = [];
    uLim = [uLim,uLimTmp]; % the interval of each toolpath
    peakPt(:,end - loopPtNumLast + 1:end) = [];
    peakPt = [peakPt,[peakPtInTmp;peakPtOutTmp]];
    res(:,end - loopPtNumLast + 1:end) = [];
    res = [res,resTmp];

    % debug
    % plot3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
    %     '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
    % quiver3(toolPathPtTmp(1,:),toolPathPtTmp(2,:),toolPathPtTmp(3,:), ...
    %     toolNormDirectTmp(1,:),toolNormDirectTmp(2,:),toolNormDirectTmp(3,:));
    % for jj = accumPtNum(end - 1) + 1:accumPtNum(end)
    %     toolSp1 = toolSp;
    %     toolSp1.coefs = quat2rotm(toolQuat(jj,:))*toolCoefs + toolVec(:,jj);
    %     Q = fnval(toolSp1,uLim(1,jj):0.01:uLim(2,jj));
    %     plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
    % end
   
    fprintf('r remains %f.\n\n',r1);
    fprintf('No.%d\tElapsed time is %f seconds.\n-----\n',length(loopPtNum),toc);
    if r1 > rMax, break; end
%     if toolREach(end) - toolREach(end-1) > 0
%         delta = toolREach(end) - toolREach(end-1);
%     else
%         delta = delta0;
%     end
    r0 = r1 + delta0;
%     iter = 1;
end

tRes = toc(tRes0);
fprintf('The time spent in the residual height calculating process is %fs.\n',tRes);

end

