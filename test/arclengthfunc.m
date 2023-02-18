function arclen = arclengthfunc(f,t)
%ARCLENGTHFUNC Summary of this function goes here
%
%   t remains the paramter of the parametric function f
%   f can be 'pp' form or function handle
%
% the arclen will be the sym variable

switch class(f)
    case {'function_handle','sym'}
        % the input is an analytical function
        f_t = diff(f,t);
        arclen = sqrt(sum(f_t.^2));
    case 'struct'
        % the input is a structural function
        if isfield(f,'form')
            switch f.form
                case 'pp'
                    % piecewise polynomial
                    
                case 'B-'
            end
        end
    case 'double'
        % the input is an array representing a polynomial function
        error('This method has not been realized.');
end


end