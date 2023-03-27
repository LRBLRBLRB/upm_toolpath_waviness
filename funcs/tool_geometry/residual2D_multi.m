function [res,peakPt,interPt,uLim1,uLim2] = residual2D_multi(sp1,sp2,eps,pt1,pt2,uLim2)
% Solve the residual height between two adjacent points on the tool path in
% 2-dimension plane, supposing that the residual height is the distance
% between the intersection point of tool edges.
% 
% Usage:
% [res,peakPt,ind1,ind2] = intersect2D_multi(s1,s2,eps,p1,p2,method)
%   s1 & s2 struct  the B-form struct of the two tool edge
%   eps (1,1)  the discretization of the parameter u
%   p1 & p2 () 
%   res (1,1)  the residual within the two position
%   
%
% [res,peakPt,ind1,ind2] = intersect2D_multi(sp1,sp2,...)
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
u = 0:eps:1;
s1 = fnval(sp1,u);
s2 = fnval(sp2,u);
% eliminate the intersection points that are beyond the curve area
%   the projection of vec{pt1 -> interPt} to vec{pt1 -> pt2}
sRange1 = dot(ndgrid(pt2 - pt1,1:size(s1,2)),s1 - pt1,1)/norm(pt2 - pt1);
isRange1 = all([sRange1 >= 0; sRange1 <= norm(pt2 - pt1)],1);
sRange2 = dot(ndgrid(pt2 - pt1,1:size(s2,2)),s2 - pt1,1)/norm(pt2 - pt1);
isRange2 = all([sRange2 >= 0; sRange2 <= norm(pt2 - pt1)],1);
% % index0: the point s1(:,index0(ii)) is the closest point from s2(:,ii)
% % dist0(ii): the closest distance between s2(:,ii) and s1(:,index0(ii))
% [index0,dist0] = dsearchn(s1',s2');
% range0 = mean(vecnorm(diff(s1,1,2),2,1));
% ind2(1) = find(dist0 < 1.5*range0,1,"first");
% ind2(2) = find(dist0 < 1.5*range0,1,"last");
% ind1(1) = index0(ind2(1));
% ind1(2) = index0(ind2(2));
u1 = u(isRange1); % the intersection range in sp1 
u2 = u(isRange2); % the intersection range in sp2
s1 = s1(:,isRange1);
s2 = s2(:,isRange2);

% same as the above dsearchn function
[index,dist] = dsearchn(s1',s2');
% get the indices of the local minimum points
cmp = find(diff(sign(diff(dist))) > 0) + 1; % length = size(s2,2)

range = mean(vecnorm(diff(s1,1,2),2,1));
% cmp2: the indices of the intersection points, 
%   i.e., the local minimals which is close enough to both sp1 and sp2
cmp = cmp(dist(cmp) < range);

%% eliminate the redundant cases
diff1 = diff(s1,1,2);
der1 = diff1(3,:)./diff1(1,:);
diff2 = diff(s2,1,2);
der2 = diff2(3,:)./diff2(1,:);
% edxclude the case that two local minimals are the same
MultiInter = find(diff(cmp) <= 10);
[num,ari] = arithmeticsequences(MultiInter);
cmp2clear = [];
for jj = 1:num
    % mid: the middle point among the same intersection point set
    midInter = round((MultiInter(ari(1,jj)) + MultiInter(ari(2,jj)) + 1)/2); 
    same = MultiInter(ari(1,jj)):MultiInter(ari(2,jj)) + 1;
    [minInter,sameInd] = min(abs(der1(index(cmp(same))) - der2(cmp(same))));
    minInd = same(sameInd);
    maxInter = max(abs(der1(index(cmp(same))) - der2(cmp(same))));
    if maxInter/minInter > 5 % abs(der1(index(cmp2(midInter))) - der2(cmp2(midInter))) < 1e-3
        % the repeated points are tagent points
        cmp2clear = [cmp2clear,MultiInter(ari(1,jj)):MultiInter(ari(2,jj)) + 1];
    else
        % the repeated points are not tagent points
        cmp2clear = [cmp2clear,MultiInter(ari(1,jj)):minInd - 1,minInd + 1:MultiInter(ari(2,jj)) + 1];
    end
end

% for ii = 1:length(cmp2)
%     % 比较每个交点前后10个点处两条tool tip的切线斜率之差
%     abs(der1(index(cmp2(ii))) - der2(cmp2(ii)))
% end
cmp(cmp2clear) = [];

% eliminate the tagent case
cmptan = abs(der1(index(cmp)) - der2(cmp));
meantan = mean(cmptan);
istan0 = find(meantan./cmptan > 5); % test the differentiate of each intersection point
istan = zeros(1,length(istan0));
for ii = 1:length(istan0)
    vec0 = s1(:,index(cmp(ii))) - s2(:,cmp(ii));
    vec1 = s1(:,index(cmp(ii) + round(1e-3/eps))) - s2(:,cmp(ii) + round(1e-3/eps));
    vec2 = s1(:,index(cmp(ii) - round(1e-3/eps))) - s2(:,cmp(ii) - round(1e-3/eps));
    istan(ii) = dot(vec0,vec1) * dot(vec0,vec2) > 0;
end
cmp(istan0) = [];

% calculate the intersection points
uInter2 = u2(cmp);
uInter1 = u1(index(cmp));
interPt = 0.5*(s1(:,index(cmp)) + s2(:,cmp));

% % eliminate the intersection points that are beyond the curve area
% %   the projection of vec{pt1 -> interPt} to vec{pt1 -> pt2}
% interRange = dot(ndgrid(pt2 - pt1,1:length(cmp)),interPt - pt1,1)/norm(pt2 - pt1);
% % if the projection is beyond the interval [0,norm(pt2 - pt1)], then eliminate
% isRange = any([interRange < 0;interRange > norm(pt2 - pt1)],1);
% cmp(isRange) = [];
% uInter2(isRange) = [];
% uInter1(isRange) = [];
% interPt(:,isRange) = [];

%% u range
if ~mod(length(cmp),2)
    warning('Something wrong in the intersection point calculation process.\n');
    res = 0;
    peakPt = [];
    interPt = [];
    uLim1 = [];
    uLim2 = [];
    return;
end

uLim2(2,end) = uInter2(1);
if length(uInter2) > 2
    uLim1 = [uInter1(1:2:end);uInter1(2:2:end),1];
    uLim2 = [uLim2,[uInter2(2:2:end);uInter2(3:2:end)]];
else
    uLim1 = [uInter1;1];
end

%% calculate the residual height and the coordinate of the peak point
% calculate all the curve pt of sp1 within uInter1 and sp2 wthin uInter2,
%   then pick out the one with the largest PV
peakPtTmp = [];
peakUTmp = [];
for ii = 1:size(uLim1,2)
    ind1 = u1 >= uLim1(1,ii) & u1 <= uLim1(2,ii);
    peakUTmp = [peakUTmp,u1(ind1)];
    peakPtTmp = [peakPtTmp,s1(:,ind1)];
end
kk = length(peakUTmp);
for ii = 1:size(uLim2,2)
    ind2 = u2 >= uLim2(1,ii) & u2 <= uLim2(2,ii);
    peakUTmp = [peakUTmp,u2(ind2)];
    peakPtTmp = [peakPtTmp,s2(:,ind2)];
end
[res,resInd] = max(pt2Line(peakPtTmp,pt1,pt2));
peakPt = peakPtTmp(:,resInd);
peakPt(4) = peakUTmp(resInd);
peakPt(5) = -1*(resInd > kk);

%% plot the intersection results
fig = figure('WindowState','maximized');
tiledlayout(2,2);
nexttile([1 2]);
plot(s1(1,:),s1(3,:),'--','LineWidth',1.5,'Color',[0,0.4470,0.7410]); % tool edge 1
hold on;
plot(s2(1,:),s2(3,:),'--','LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]); % tool edge 2
axis equal;
ax = get(gca,'XLim');
ay = get(gca,'YLim');
curveSlope = (pt2(3) - pt1(3))/(pt2(1) - pt1(1));
plot([pt1(1),pt1(1) - curveSlope*10],[pt1(3),pt1(3) + 10], ...
    '--','Color',[0.4660 0.6740 0.1880]); % the intersection range
plot([pt2(1),pt2(1) - curveSlope*10],[pt2(3),pt2(3) + 10], ...
    '--','Color',[0.4660 0.6740 0.1880]); % the intersection range
plot([pt1(1),pt2(1)],[pt1(3),pt2(3)],'LineWidth',1,'Color',[0.4660 0.6740 0.1880], ...
    'Marker','.','MarkerSize',18); % theoretical curve edge in line
plot(peakPtTmp(1,1:kk),peakPtTmp(3,1:kk),'.', ...
    'LineWidth',1.5,'Color',[0,0.4470,0.7410]); % actual curve edge
plot(peakPtTmp(1,kk + 1:end),peakPtTmp(3,kk + 1:end),'.', ...
    'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250]); % actual curve edge
plot(interPt(1,:),interPt(3,:),'.','MarkerSize',36, ...
    'Color',[0.8500 0.3250 0.0980]); % intersection points
plot(peakPt(1),peakPt(3),'.','MarkerSize',36, ...
    'Color',[0.4940 0.1840 0.5560]); % peak point
set(gca,'XLim',ax,'YLim',ay);
title('tool tip edge');
xlabel('x'); zlabel('z');
legend('edge 1','edge 2','intersection range','','theoretical surface', ...
    'actual surface','intersection pt','peak pt','Location','bestoutside');
set(gca,'FontSize',14,'FontName','Times New Roman');
nexttile(3);
plot(1:length(dist),dist);
hold on;
line(1:length(dist),zeros(1,length(dist)),'LineStyle','--','LineWidth',0.3, ...
    'Color',[0.8500 0.3250 0.0980]);
title('dist-No.');
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'XDir','reverse');
nexttile(4);
plot(u2,dist);
hold on;
line(u2,zeros(1,length(u2)),'LineStyle','--','LineWidth',0.3, ...
    'Color',[0.8500 0.3250 0.0980]);
title('dist-u2');
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'XDir','reverse');
% pause(1);
delete(fig);

end