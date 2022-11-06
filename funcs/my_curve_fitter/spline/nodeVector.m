function U = nodeVector(k,n,options)
% solve the node vector of the spline
% Input:
%   the order of the spline function: k
%   the number of control pts: n
%  options:
%   node generation method: options.nodeMethod
%   the coordinates of control pts (nodeMethod = 'Riesenfeld'/'Hartley-Judd'): options.cpts
%   the param vector of interpolation (nodeMethod = 'Interpolation'): options.cpts
% Output:
%   the node vector: U

arguments
    k (1,1) {mustBeInteger,mustBePositive}
    n (1,1) {mustBeInteger,mustBeGreaterThanOrEqual(n,k)}
    options.NodeMethod {mustBeMember(options.NodeMethod, ...
        {'Riesenfeld','Hartley-Judd','Interpolation','Approximation',''})} = 'Interpolation'
    options.Cpts {mustBeFloat} = []
    options.uQ (:,1) {mustBeInRange(options.uQ,0,1,'inclusive')} = []
end

if n-k-1<0
    error('control points are not enough, or the required order is too large');
end

switch options.NodeMethod
    case 'Riesenfeld'
        l = sum((options.Cpts(1:end-1,:)-options.Cpts(2:end,:)).^2,2).^(1/2);
        L = sum(l);
        UTemp = zeros(n-k-1,1);
        if mod(k,2) == 0
            tmp = sum(l(1:k/2));
            for ii = 1:n-k-1
                tmp = tmp + l((k+2*ii)/2);
                UTemp(ii) = tmp/L;
                tmp = tmp + l((k+2*ii+2)/2);
            end
        else
            tmp = sum(l(1:(k-1)/2));
            for ii = 1:n-k-1
                tmp = tmp + l((k+2*ii-1)/2);
                UTemp(ii) = tmp/L;
            end
        end
        U = [zeros(k+1,1); UTemp; ones(k+1,1)];
    case 'Hartley-Judd'
        l = sum((options.Cpts(1:end-1,:)-options.Cpts(2:end,:)).^2,2).^(1/2);
        UTemp = zeros(n-k-1,1);
        deltaU = zeros(n-k,1);
        for ii = 1:k
            deltaU = deltaU + l(ii:(n-k+ii-1));
        end
        deltaU = deltaU/sum(deltaU);
        tmp = 0;
        for ii = 1:n-k-1
            tmp = tmp + deltaU(ii);
            UTemp(ii) = tmp;
        end
        % normalized node vector
        U = [zeros(k+1,1); UTemp; ones(k+1,1)];
    case 'Interpolation' % interpolation node vector generation
        U = [zeros(n,1); ones(k+1,1)]; % normalized node vector
        for ii = k+2:n
            for j = 1:k
                U(ii) = U(ii) + options.uQ(ii-k+j-1)/k;
            end
        end
    case 'Approximation'
end

end