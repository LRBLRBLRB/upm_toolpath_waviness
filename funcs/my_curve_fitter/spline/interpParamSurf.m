function [uQ,vQ] = interpParamSurf(Q,options)
% usage: 
% The function is used to process the parameterization during
% interpolation. 
% Inputs:
%   Q (m,n,3) matrix of data points
%   interpMethod interpolation method
% Outputs:
%   uQ (m,1) vQ (n,1) node parameters corresponding to Q

arguments
    Q {mustBeFinite}
    options.paramMethod {mustBeMember(options.paramMethod, ...
        ['uniform','concentric','chord',''])} = 'chord'
end

[m,n,~] = size(Q);

if strcmp(options.paramMethod,'uniform') % uniform parameterization
    uQ = transpose(linspace(0,1,m));
    vQ = transpose(linspace(0,1,n));
else
    if strcmp(options.paramMethod,'concentric') % concentric parameterization
        % sqrt of dist of adjacent Q, size(n-1,1)
        uDist = sum((Q(2:end,:,:)-Q(1:end-1,:,:)).^2,3).^(1/4); % (m-1,n)
        vDist = sum((Q(:,2:end,:)-Q(:,1:end-1,:)).^2,3).^(1/4); % (m,n-1)
    else % chord parameterization
        % dist of adjacent Q, size(n-1,1)
        uDist = sum((Q(2:end,:,:)-Q(1:end-1,:,:)).^2,3).^(1/2); % (m-1,n)
        vDist = sum((Q(:,2:end,:)-Q(:,1:end-1,:)).^2,3).^(1/2); % (m,n-1)
    end
    uSum = sum(uDist,1); % (1,n)
    vSum = sum(vDist,2); % (m,1)
    % 1st param of Q
    uQp = zeros(m,n);
    for ii = 2:m
        uQp(ii,:) = uQp(ii-1,:) + uDist(ii-1,:);
    end
    uQp = uQp./uSum;
    uQp(isnan(uQp)) = 0; 
    uQ = mean(uQp,2);
    uQ(end) = 1; % elliminate calcuation error
    % 2nd param of Q
    vQq = zeros(m,n);
    for jj = 2:n
        vQq(:,jj) = vQq(:,jj-1) + vDist(:,jj-1);
    end
    vQq = vQq./vSum;
    vQq(isnan(vQq)) = 0;
    vQ = transpose(mean(vQq,1));
    vQ(end) = 1;
end

end