function [tEq,fEq,uEq,varargout] = arclengthparam(arcInv,maxAng,t,f,u,quat,vec,options)
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
    quat (:,4) = []
    vec cell = {}
    options.algorithm {mustBeMember(options.algorithm, ...
        {'analytical-function','list-interpolation','lq-fitting'})} ...
        = 'list-interpolation'
    options.interpType {mustBeMember(options.interpType, ...
        {'linear','nearest','next','previous','pchip','cubic,' ...
        'quadratic'})} = 'linear'
    options.uQTol {mustBePositive} = 1e-3
end

% find the point needed to be recalculated (for constant arc length)
angDist = abs((t(2:end) - t(1:end - 1)) - maxAng);
angDistEqEps = 1e-6; % the epsilon to judge whether the angle difference is equal or not
angularMaxInd = find(angDist > angDistEqEps,1);
tEq(1:angularMaxInd) = t(1:angularMaxInd);
fEq(:,1:angularMaxInd) = f(:,1:angularMaxInd);
uEq(1:angularMaxInd) = u(1:angularMaxInd);
% t(1:angularMaxInd) = [];
% f(:,1:angularMaxInd) = [];

%% position interpolation
switch options.algorithm
    case 'analytical-function'
        % f remains the curve parametric function: [x,y,z]=f(t)
    case 'list-interpolation'
        % or f remains the scatters of the curve: (3,:) or (2,:)
        f_t = numericdiff(f,1,2);
        arcLen = cumtrapz(sqrt(sum(f_t.^2,1)));
        X = arcLen.'; % query arc length
        V = [t.',u.',f.']; % query point table
        Xq = arcLen(angularMaxInd + 1):arcInv:arcLen(end);% equally spaced indices
        switch options.interpType
            case {'linear','nearest','next','previous','pchip','cubic'}
                Vq = interp1(X,V,Xq,options.interpType);
            case 'quadratic'
        end
        tEq = [tEq,(Vq(:,1))'];
        uEq = [uEq,(Vq(:,2))'];
        fEq = [fEq,(Vq(:,3:end))'];
        numEq = length(tEq);
    case 'lq-fitting'
        % the algorithm in "Tencent's Game Development"
        if length(t) ~= size(f,2)
            error('the parameters and point values do not match well.');
        end
        % Step one: trajactory build
        mu = 1;
        tang = mu*(f(:,3:end) - f(:,1:end-2));
        % Step two: bi-quadratic spline
        A = [zeros(3,3), zeros(3,3),     eye(3), zeros(3,3), zeros(3,3), zeros(3,3);
             zeros(3,3),     eye(3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
                 eye(3),     eye(3),     eye(3), zeros(3,3), zeros(3,3),  -1*eye(3);
               2*eye(3),     eye(3), zeros(3,3), zeros(3,3),  -1*eye(3), zeros(3,3);
             zeros(3,3), zeros(3,3), zeros(3,3),     eye(3),     eye(3),     eye(3);
             zeros(3,3), zeros(3,3), zeros(3,3),   2*eye(3),     eye(3), zeros(3,3)];
        b = [f(:,1:end - 1);tang(:,1:end - 1);0;0;f(:,2:end);tang(:,2:end)];
        coef = A\b; % column: a1 b1 c1 a2 b2 c2; row: 1 ~ num - 1
        % Step three: constant-chord segmentation
        
        % tool contact u calculation
        spiralContactU = spiralcontactu(spiralPath,spiralQuat,toolSp,surfFunc);
end

%% orientation interpolation
% figure;
% plot3(fEq(1,1:angularMaxInd),fEq(2,1:angularMaxInd),fEq(3,1:angularMaxInd),'b.');
% hold on;

vecEq = {};
if isempty(quat)
    quatEq = [];
else
    quatEq = zeros(numEq,4);
    vecLen = length(vec);
    vecEq = cell(vecLen,1);
    % constant angle ones
    quatEq(1:angularMaxInd,:) = quat(1:angularMaxInd,:);
    for ii = 1:vecLen
        vecEq{ii}(:,1:angularMaxInd) = vec{ii}(:,1:angularMaxInd);
    end
    % constant arc length ones
    for ii = angularMaxInd + 1:numEq

%         plot3(fEq(1,ii),fEq(2,ii),fEq(3,ii),'b.');

        ind1 = find(tEq(ii) >= t);
        ind1 = ind1(end); % the closest parameter smaller than tEq
        ind2 = find(tEq(ii) < t,1); % the closest parameter larger than tEq
        quat_t = (tEq(ii) - t(ind1))/(t(ind2) - t(ind1));
        quatEq(ii,:) = slerp(quat(ind1,:),quat(ind2,:),quat_t);
        for kk = 1:vecLen
            vecEq{kk}(:,ii) = quat2rotm(quatEq(ii,:))*vec{kk}(:,ind1);
%             q(kk) = quiver3(fEq(1,ii),fEq(2,ii),fEq(3,ii), ...
%                 vecEq{kk}(1,ii),vecEq{kk}(2,ii),vecEq{kk}(3,ii), ...
%                 50,'ShowArrowHead','on');
        end
%         delete(q(1));
%         delete(q(2));
    end
end



%% output
switch nargout
    case 3
    case 4
        varargout{1} = quatEq;
    case 5
        varargout{1} = quatEq;
        varargout{2} = vecEq;
end
end