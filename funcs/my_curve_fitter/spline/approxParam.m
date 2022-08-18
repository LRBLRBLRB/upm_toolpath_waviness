function approxParam(Q,options)
% usage: 
% The function is used to process the parameterization during
% approximation. 
% Inputs:
%   Q (n,3) matrix of data points
%   paramMethod {'uniform','concentric','chord'} parameterization method
% Outputs:
%   uQ (n,1) node parameters corresponding to Q

arguments
    Q {mustBeFinite}
    options.paramMethod {mustBeMember(options.paramMethod, ...
        [''])} = ''
end

n = size(Q,1);



end