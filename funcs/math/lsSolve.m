function X = lsSolve(A,B)
% usage: X = LSC(A,B)
% to solve the least square solution of over-determined equations
% Input:
%   A (m,n) 
%   B (m,1)
% Output:
%   X (n,1)

arguments
    A (:,:)
    B (:,1)
end

[m,~] = size(A);
if size(B,1) ~= m
    error('Error: the size of B does not match A');
end

% X = zeros(n,1);
X = (A'*A)\(A'*B);

end
