function pts = bSplinePts(cpt,k,u)
% DeBoor-Coox to solve all the pts on B-spline curve with param u
% Inputs:
%   control points: cpt (n,dim)
%   order pf B-spline: k
%   parameter vector: u (:,1)
% Outputs:
%   curve points: pts (:,dim)

% --------------STATE: still need to debug--------------

for i = 1:length(u)
    j = round(u(i));
    for ii = j-k:j
        pts = cpt;
    end
    for r = 1:k
        for ii = j:j-k+r
            pts(ii,:) = (t-ii)/(k-r+1)*pts(ii,:) + (ii+k-r+1-t)/(k-r+1)*pts(ii-1,:);
        end
    end
end
end
