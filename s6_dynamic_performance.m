% test whether the dynamic performances is good or not

%% concentric acceleration check
acce = diff(diff(curvePathPt,1,2),1,2);
maxAcce = max(acce);
figure;
plot(curvePathPt(1,2:end - 1),acce(3,:)./acce(1,:));
xlabel(['r (',unit,')']);
ylabel(['a (',unit,'^{-1})']);

%% 