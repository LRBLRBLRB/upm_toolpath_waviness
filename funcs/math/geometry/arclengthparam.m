function [tEq,fEq] = arclengthparam(t,f,arcInr,options)
%ARCLENGTHPARAM Generate the function of arc length 
%
%
% 1st calculation method: numerical table query
%   The vector t remains the 

arguments
    t
    f
    arcInr {mustBeNumeric}
    options.interpType {mustBeMember(options.interpType, ...
        {'linear','nearest','next','previous','pchip','cubic,' ...
        'quadratic'})} = 'linear'
end


if isa(f,'function_handle')
    % f remains the curve parametric function: [x,y,z]=f(t)
else
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
end

end