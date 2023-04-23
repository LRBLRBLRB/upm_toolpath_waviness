function pts = bSplineSurfPts(cpts,k,u,U,l,v,V)
% usage: pt = bSplineSurfacePt(cpts,k,u,U,l,v,V)
% solve a series of points of B-spline surface with param vector u,v
% Input：
%   control points: cpts (m+1,n+1,dimention)
%   order of B-spline: k, l
%   node parameter vector: u (:,1), v (:,1)
%   node vector: U (n+k+1,1), V (m+l+1,1)
% Output:
%   B-spline surface points: pts

dim = size(cpts,3);
uNum = length(u);
vNum = length(v);
pts = zeros(uNum,vNum,dim);

%% method one: reverse all the pts in vector u & v
% ptTmp = zeros(dim,uNum*vNum);
for ii = 1:uNum
    for jj = 1:vNum
        ptTmp = bSplineSurfPt(cpts,k,u(ii),U,l,v(jj),V);
        pts(ii,jj,:) = reshape(ptTmp,[1 1 dim]);
        % ptTmp(:,(ii-1)*vNum+jj) = bSplineSurfPt(cpts,k,u(ii),U,l,v(jj),V);
    end
end

% for kk = 1:dim
%     pts(:,:,kk) = reshape(ptTmp(kk,:),[uNum,vNum]);
% end
% clear ptTmp;

end
