function [res,peakPt,uLim1,uLim2] = residual2D_multi(sp1,sp2,eps,uLim2)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,peakPt,ind1,ind2] = residual2D_numeric(s1,s2,eps,p1,p2,method)
%   s1 & s2 struct  the B-form struct of the two tool edge
%   eps (1,1)  the discretization of the parameter u
%   p1 & p2 () 
%   res (1,1)  the residual within the two position
%   
%
% [res,peakPt,ind1,ind2] = residual2D_numeric(sp1,sp2,...)
%   all the same except that the sp1 and sp2 remain the B-form spline
%   struct of the tool edge

%% solve the intersection point of the two spline curve
% switch method
%     case 'pt2surf'
%         if dim == 2
%             dist1 = dist2curve(sp1,surfFunc);
%             dist2 = dist2curve(sp2,surfFunc);
%         else
%             dist1 = dist2surf(sp1,surfFunc);
%             dist2 = dist2surf(sp2,surfFunc);
%         end
%         u1 = dist1 <= dist2; % the closer one will be selected
%         u2 = dist1 > dist2;
%         cmp = dist1 <= dist2;
%     case 'DSearchn'

% find the intersection range
eps0 = 1e-3;
if eps < eps0
    u = 0:eps0:1;
    s1 = fnval(sp1,u);
    s2 = fnval(sp2,u);
    % index0: the point s1(:,index0(ii)) is the closest point from s2(:,ii)
    % dist0(ii): the closest distance between s2(:,ii) and s1(:,index0(ii))
    [index0,dist0] = dsearchn(s1',s2');
    range0 = mean(vecnorm(diff(s1,1,2),2,1));
    ind2(1) = find(dist0 < 1.5*range0,1,"first");
    ind2(2) = find(dist0 < 1.5*range0,1,"last");
    ind1(1) = index0(ind2(1));
    ind1(2) = index0(ind2(2));
    u1 = u(ind1(1)):eps:u(ind1(2)); % the intersection range in sp1 
    u2 = u(ind2(1)):eps:u(ind2(2)); % the intersection range in sp2
else
    u1 = 0:eps:1;
    u2 = 0:eps:1;
end
% find the intersection points
s1 = fnval(sp1,u1);
s2 = fnval(sp2,u2);
% same as the above dsearchn function
[index,dist] = dsearchn(s1',s2');
% get the indices of the local minimum points
cmp1 = find(diff(sign(diff(dist))) > 0) + 1;
range = mean(vecnorm(diff(s1,1,2),2,1));
% get the indices of the intersection points, 
% i.e., the local minimals which is close enough to both sp1 and sp2
cmp2 = cmp1(dist(cmp1) < range);
% exclude the tagent case
diff1 = diff(s1,1,2);
der1 = diff1(3,:)./diff1(1,:);
diff2 = diff(s2,1,2);
der2 = diff2(3,:)./diff2(1,:);
for ii = 1:length(cmp2)
    % 比较每个交点前后10个点处两条tool tip的切线斜率之差
    abs(der1(index(cmp2(ii))) - der2(cmp2(ii)))
end
cmp2(diff(cmp2) <= 10) = []; % edxclude the case that two local minimals are the same

uInter2 = u2(cmp2);
uInter1 = u1(index(cmp2));
res = dist(cmp2);
peakPt = 0.5*(s1(:,index(cmp2)) + s2(:,cmp2));

fig = figure('WindowState','maximized');
tiledlayout(2,2);
nexttile(1,[1 2]);
plot(s1(1,:),s1(3,:),'LineWidth',1.5,'Color',[0,0.4470,0.7410]);
hold on;
plot(s2(1,:),s2(3,:),'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]);
plot(peakPt(1,:),peakPt(3,:),'.','MarkerSize',36,'Color',[0.8500 0.3250 0.0980]);
title('tool tip edge');
xlabel('x'); zlabel('z');
legend('edge 1','edge 2','peak pt','Location','best');
set(gca,'FontSize',14,'FontName','Times New Roman');
nexttile(3);
plot(1:length(dist),dist);
title('dist-No.');
set(gca,'FontSize',14,'FontName','Times New Roman');
nexttile(4);
plot(u2,dist);
title('dist-u2');
set(gca,'FontSize',14,'FontName','Times New Roman');
% pause(1);
delete(fig);

if ~mod(length(cmp2),2)
    warning('Something wrong in the intersection point calculation process.\n');
    pause();
    return;
end

uLim2(2,end) = uInter2(1);
if length(uInter2) > 2
    uLim1 = [uInter1(1:2:end);uInter1(2:2:end),1];
    uLim2 = [uLim2,[uInter2(2:2:end);uInter2(3:2:end)]];
else
    uLim1 = [uInter1;1];
end

end