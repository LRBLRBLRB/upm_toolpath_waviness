function [tEq,fEq] = arclengthparam(t,f,arcInr,options)
%ARCLENGTHPARAM Generate the function of arc length 
%
%
% 1st calculation method: numerical table query
%   The vector t remains the 

arguments
    t
    f (3,:)
    arcInr {mustBeNumeric}
    options.algorithm {mustBeMember(options.algorithm, ...
        {'analytical-function','list-interpolation','lq-fitting'})} ...
        = 'list-interpolation'
    options.interpType {mustBeMember(options.interpType, ...
        {'linear','nearest','next','previous','pchip','cubic,' ...
        'quadratic'})} = 'linear'
end

switch options.algorithm
    case 'analytical-function'
        % f remains the curve parametric function: [x,y,z]=f(t)
    case 'list-interpolation'
        % or f remains the scatters of the curve: (3,:) or (2,:)
        f_t = numericdiff(f,1,2);
        arcLen = cumtrapz(sqrt(sum(f_t.^2,1)));
        
        X = arcLen.'; % query arc length
        V = [t.',f.']; % query point table
        L = arcLen(end); % total arc length
        Xq = 0:arcInr:L;% equally spaced indices
        switch options.interpType
            case {'linear','nearest','next','previous','pchip','cubic'}
                Vq = interp1(X,V,Xq,options.interpType);
            case 'quadratic'
        end
        
        tEq = (Vq(:,1))';
        fEq = (Vq(:,2:end))';
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
        
end

end