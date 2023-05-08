function [tOut,fOut,uOut,varargout] = arclengthparam(arcInv,maxAng,t,f,u,var1,var2,options)
%ARCLENGTHPARAM Generate the function of arc length 
%
%
% 1st calculation method: numerical table query
%   The vector t remains the 

arguments
    arcInv {mustBeNumeric}
    maxAng {mustBeNumeric}
    t
    f 
    u (1,:) = []
    var1 = []
    var2 = []
    options.algorithm {mustBeMember(options.algorithm, ...
        {'analytical-function','list-interpolation','lq-fitting'})} ...
        = 'list-interpolation'
    options.interpType {mustBeMember(options.interpType, ...
        {'linear','nearest','next','previous','pchip','cubic,' ...
        'quadratic'})} = 'linear'
    options.uQTol {mustBePositive} = 1e-3
    options.cutDirection {mustBeMember(options.cutDirection, ...
        {'Edge to Center','Center to Edge'})} = 'Edge to Center'
end

% find the point needed to be recalculated (for constant arc length)
angDist = abs((t(2:end) - t(1:end - 1)) - maxAng);
angDistEqEps = 1e-6; % the epsilon to judge whether the angle difference is equal or not

switch options.cutDirection
    case 'Edge to Center'
        angularMaxInd = find(angDist < angDistEqEps,1,'first');
        indChange = 1:angularMaxInd; % constant-arc range
        indOri = (angularMaxInd + 1):length(t); % constant-angle range
    case 'Center to Edge'
        angularMaxInd = find(angDist > angDistEqEps,1,'first');
        indChange = (angularMaxInd + 1):length(t);
        indOri = 1:angularMaxInd;
end
tEq = t(indChange);
fEq = f(:,indChange);
uEq = u(indChange);

% t(1:angularMaxInd) = [];
% f(:,1:angularMaxInd) = [];

%% position interpolation
switch options.algorithm
    case 'analytical-function'
        % f remains the curve parametric function: [x,y,z]=f(t)
    case 'list-interpolation'
        % or f remains the scatters of the curve: (3,:) or (2,:)
        f_t = numericdiff(fEq,1,2);
        arcLen = cumtrapz(sqrt(sum(f_t.^2,1)));
        X = arcLen.'; % query arc length
        V = [tEq.',uEq.',fEq.']; % query point table
        Xq = linspace(arcLen(1),arcLen(end), ...
            ceil((arcLen(end) - arcLen(1))/arcInv));% equally spaced indices
        switch options.interpType
            case {'linear','nearest','next','previous','pchip','cubic'}
                Vq = interp1(X,V,Xq,options.interpType);
            case 'quadratic'
        end
        tOut = (Vq(:,1))';
        uOut = (Vq(:,2))';
        fOut = (Vq(:,3:end))';
        numOut = length(tOut);
    case 'lq-fitting'
        % the algorithm in "Tencent's Game Development"
        if length(tEq) ~= size(fEq,2)
            error('the parameters and point values do not match well.');
        end
        % Step one: trajactory build
        mu = 1;
        tang = mu*(fEq(:,3:end) - fEq(:,1:end-2));
        % Step two: bi-quadratic spline
        A = [zeros(3,3), zeros(3,3),     eye(3), zeros(3,3), zeros(3,3), zeros(3,3);
             zeros(3,3),     eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
                 eye(3),     eye(3),     eye(3), zeros(3,3), zeros(3,3),  -1*eye(3);
               2*eye(3),     eye(3), zeros(3,3), zeros(3,3),  -1*eye(3), zeros(3,3);
             zeros(3,3), zeros(3,3), zeros(3,3),     eye(3),     eye(3),     eye(3);
             zeros(3,3), zeros(3,3), zeros(3,3),   2*eye(3),     eye(3), zeros(3,3)];
        b = [fEq(:,1:end - 1);tang(:,1:end - 1);0;0;fEq(:,2:end);tang(:,2:end)];
        coef = A\b; % column: a1 b1 c1 a2 b2 c2; row: 1 ~ num - 1
        % Step three: constant-chord segmentation
        
        % tool contact u calculation
        spiralContactU = spiralcontactu(spiralPath,spiralQuat,toolSp,surfFunc);
end

%% orientation interpolation
% figure;
% plot3(fEq(1,1:angularMaxInd),fEq(2,1:angularMaxInd),fEq(3,1:angularMaxInd),'b.');
% hold on;

vecOut = {};
if isempty(var1) % input: (arcInv,maxAng,t,f,u,options)
    quatOut = [];
elseif size(var1,2) == 4 % input: (arcInv,maxAng,t,f,u,quat,vec,options)
    quat = var1;
    quatEq = quat(indChange,:);
    quat = quat(indOri,:); % quat that only have constant-angle
    vec = var2;
    vecLen = length(vec);
    for ii = 1:vecLen
        vecEq{ii} = vec{ii}(:,indChange);
    end
    quatOut = zeros(numOut,4);
    vecOut = cell(vecLen,1);
    % constant angle ones
%     quatOut(1:angularMaxInd,:) = quat(1:angularMaxInd,:);
%     for ii = 1:vecLen
%         vecOut{ii}(:,1:angularMaxInd) = vec{ii}(:,1:angularMaxInd);
%     end
    % constant arc length ones
    for ii = 1:numOut - 1

%         plot3(fOut(1,ii),fOut(2,ii),fOut(3,ii),'b.');

        ind1 = find(tOut(ii) >= tEq);
        ind1 = ind1(end); % the closest parameter smaller than tEq
        ind2 = find(tOut(ii) < tEq,1); % the closest parameter larger than tEq
        quat_t = (tOut(ii) - tEq(ind1))/(tEq(ind2) - tEq(ind1));
        quatOut(ii,:) = slerp(quatEq(ind1,:),quatEq(ind2,:),quat_t);
        for kk = 1:vecLen
            vecOut{kk}(:,ii) = quat2rotm(quatOut(ii,:))*vecEq{kk}(:,ind1);

%             q(kk) = quiver3(fOut(1,ii),fOut(2,ii),fOut(3,ii), ...
%                 vecOut{kk}(1,ii),vecOut{kk}(2,ii),vecOut{kk}(3,ii), ...
%                 50,'ShowArrowHead','on');
        end
%         delete(q(1));
%         delete(q(2));
    end

    % the last point
    quatOut(numOut,:) = quatEq(end,:);
    for kk = 1:vecLen
        vecOut{kk}(:,numOut) = vecEq{kk}(:,end);
    end

elseif iscell(var1) % input: (arcInv,maxAng,t,f,u,vec,toolData,options)
    vec = var1;
    vecLen = length(vec);
    for ii = 1:vecLen
        vecEq{ii} = vec{ii}(:,indChange);
    end
    quatOut = zeros(numOut,4);
    vecOut = cell(vecLen,1);
    % constant angle ones
    for ii = indOri
        quat(ii,:) = rotm2quat(axesRot(var2.toolEdgeNorm,var2.cutDirect, ...
            vec{1}(:,ii),vec{2}(:,ii),''));
    end

    % constant arc length ones
    for ii = 1:numOut - 1

%         plot3(fEq(1,ii),fEq(2,ii),fEq(3,ii),'b.');

        ind1 = find(tOut(ii) >= tEq);
        ind1 = ind1(end); % the closest parameter smaller than tEq
        ind2 = find(tOut(ii) < tEq,1); % the closest parameter larger than tEq
        quat_t = (tOut(ii) - tEq(ind1))/(tEq(ind2) - tEq(ind1));
        quat_12 = rotm2quat(axesRot(vecEq{1}(:,ind1),vecEq{2}(:,ind1), ...
            vecEq{1}(:,ind2),vecEq{2}(:,ind2),''));
        quat_1t = slerp([1,0,0,0],quat_12,quat_t);
        for kk = 1:vecLen
%             quat = vecQuat(vec{kk}(:,ind1),vec{kk}(:,ind2));
%             quatEq(ii,:) = slerp([1,0,0,0],quat,quat_t);
            vecOut{kk}(:,ii) = quat2rotm(quat_1t)*vecEq{kk}(:,ind1);

%             q(kk) = quiver3(fEq(1,ii),fEq(2,ii),fEq(3,ii), ...
%                 vecEq{kk}(1,ii),vecEq{kk}(2,ii),vecEq{kk}(3,ii), ...
%                 50,'ShowArrowHead','on');
        end
        quatOut(ii,:) = rotm2quat(axesRot(var2.toolEdgeNorm,var2.cutDirect, ...
            vecOut{1}(:,ii),vecOut{2}(:,ii),''));
%         delete(q(1));
%         delete(q(2));
    end

    % the last point
    quatOut(numOut,:) = rotm2quat(axesRot(var2.toolEdgeNorm,var2.cutDirect, ...
        vecEq{2}(:,end),vecEq{2}(:,end),''));
    for kk = 1:vecLen
        vecOut{kk}(:,numOut) = vecEq{kk}(:,end);
    end
end



%% output
switch options.cutDirection
    case 'Edge to Center'
        tOut = [tOut,t(indOri)];
        fOut = [fOut,f(:,indOri)];
        uOut = [uOut,u(indOri)];
        if exist('vecOut','var')
            for ii = 1:vecLen
                vecOut{ii} = [vecOut{ii},vec{ii}(:,indOri)];
            end
            if exist('quatOut','var')
                quatOut = [quatOut;quat(indOri,:)];
            end
        end
    case 'Center to Edge'
        tOut = [t(indOri),tOut];
        fOut = [f(:,indOri),fOut];
        uOut = [u(indOri),uOut];
        if exist('vecOut','var')
            for ii = 1:vecLen
                vecOut{ii} = [vec{ii}(:,indOri),vecOut{ii}];
            end
            if exist('quatOut','var')
                quatOut = [quat;quatOut];
            end
        end
end

switch nargout
    case 3
    case 4
        varargout{1} = quatOut;
    case 5
        varargout{1} = quatOut;
        varargout{2} = vecOut;
end
end