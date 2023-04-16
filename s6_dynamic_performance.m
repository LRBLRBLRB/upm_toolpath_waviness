% test whether the dynamic performances is good or not

%% concentric acceleration check
acce = diff(diff(curvePathPt,1,2),1,2);
maxAcce = max(acce);
figure;
plot(curvePathPt(1,2:end - 1),acce(3,:)./acce(1,:));
xlabel(['r (',unit,')']);
ylabel(['a (',unit,'^{-1})']);

%% spiral acceration check
spiralZ = spiralPath(3,:);
spiralR = vecnorm(spiralPath(1:2,:),2,1);

ztt = diff(diff(spiralZ));
rtt = diff(diff(spiralR));

figure;
tiledlayout(2,1);
nexttile;
plot(ztt);
hold on;
ylabel(['a_{z} (',unit,'/s^2)']);
nexttile;
plot(rtt);
ylabel(['a_{r} (',unit,'/s^2)']);