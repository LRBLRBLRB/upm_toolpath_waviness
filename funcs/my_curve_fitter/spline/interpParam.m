function uQ = interpParam(Q,interpMethod)
% usage: 
% The function is used to process the parameterization during
% interpolation. 
% Inputs:
%   Q (n,3) matrix of data points
%   interpMethod interpolation method
% Outputs:
%   uQ (n,1) node parameters corresponding to Q

arguments
    Q {mustBeFinite}
    interpMethod {mustBeMember(interpMethod, ...
        ['uniform','concentric','chord',''])} = 'chord'
end

n = size(Q,1);
if strcmp(interpMethod,'uniform') % uniform parameterization
    uQ = transpose(linspace(0,1,n));
else
    if strcmp(interpMethod,'concentric') % concentric parameterization
        l = sum((Q(1:end-1,:)-Q(2:end,:)).^2,2).^(1/4); % sqrt of dist of adjacent Q, size(n-1,1)
    else % chord parameterization
        l = sum((Q(1:end-1,:)-Q(2:end,:)).^2,2).^(1/2); % dist of adjacent Q, size(n-1,1)
    end
    L = sum(l);
    uQ = zeros(n,1);
    for ii = 2:n-1
        uQ(ii) = uQ(ii-1) + l(ii-1);
    end
    uQ = uQ/L; % control points are on the spline curve
    uQ(end) = 1;
end

end