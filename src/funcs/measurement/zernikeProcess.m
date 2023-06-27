function [reconstructed,a] = zernikeProcess(data,n)
%ZERNIKEPROCESS 此处显示有关此函数的摘要
%   此处显示详细说明

% % normalization
% dataR = sqrt(data(:,:,1).^2 + data(:,:,2).^2);
% maxR = max(dataR,[],'all');
% data(:,:,1:2) = data(:,:,1:2)./maxR;

L = size(data,1);
X = -1:2/(L-1):1;
[x,y] = meshgrid(X);
x = x(:); y = y(:);
[theta,r] = cart2pol(x,y); 

N = []; M = [];
for i = 0:n
    N = [N i*ones(1,i+1)];
    M = [M -i:2:i];
end

isInCircle = (r <= 1);
isNotNan = ~isnan(data(:));
isFit = isNotNan & isInCircle;
Z = zernfun(N,M,r(isFit),theta(isFit));
% (here a householder transform can be used to avoid a pathological equation)
a = Z\data(isFit); % coeffs of Zernike 

reconstructed = NaN(size(data));
reconstructed(isFit) = Z*a; 

end