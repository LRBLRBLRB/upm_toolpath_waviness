function pt = bSplineSurfPt(cpts,k,u,U,l,v,V)
% usage: pt = bSplineSurfacePt(cpts,k,u,U,l,v,V)
% solve the point of B-spline surface with parameter u,v
% Input：
%   control points: cpts (m+1,n+1,dimention)
%   order of B-spline: k, l
%   node parameter: u, v
%   node vector: U (n+k+1,1), V (m+l+1)
% Output:
%   B-spline surface point: pt

[m,n,dim] = size(cpts);
ui = findSpan(m,k,u,U);
Nu = bBasisFuns(ui,k,u,U);
vj = findSpan(n,l,v,V);
Nv = bBasisFuns(vj,l,v,V);

pt = zeros(dim,1);
for ii = 1:dim
    pt(ii) = transpose(Nu)*cpts((ui-k):ui,(vj-l):vj,ii)*Nv;
end

% matrix calculation, which can be replaced by one line
% pt = 0;
% uind = ui - k;
% temp = zeros(l+1,1);
% for q = 0:l
%     vind = vj-l+q;
%     for p = 0:k
%         temp(q+1) = temp(q+1) + Nu(p+1)*cpt(uind+p+1)(vind+1);
%     end
%     pt = pt + Nv(q+1)*temp(q+1);
% end

end
