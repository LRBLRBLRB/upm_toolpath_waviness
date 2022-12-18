addpath(genpath('funcs'));
fitLineFcn = @(pts) arcFit3D(pts','displayType','off');  % fit function

%% test whetger the functions above is true
cx0 = 1*1000; % unit:mu m
cy0 = 2*1000; % unit:mu m
cz0 = 3*1000; % unit:mu m
r0 = 0.1*1000; % unit:mu m
openAng = pi/3; % unit: rad
edgePV = 200; % low-frequency error
k = -edgePV/openAng;
noise = r0*5e-3; % mid-frequency error
zNoise = r0*0.5; % data pre-processing error
theta = transpose(linspace(0,openAng,300));
r = r0 + edgePV/2 + k*theta + (noise*rand(length(theta),1) - 0.5*noise);
testOri(:,1) = cx0 + r.*cos(theta);
testOri(:,2) = cy0 + r.*sin(theta);
testOri(:,3) = cz0 + (zNoise*rand(length(theta),1) - 0.5*zNoise);
testOri = testOri*(rotz(pi/3))'*(roty(pi/6))';
fitCirc = fitLineFcn(testOri);

%% plot the data
figure('Name','Function Testification');
plot3(testOri(:,1),testOri(:,2),testOri(:,3),'.','Color',[0,0.45,0.74]); hold on;

% plot the fitting center of the circle
scatter3(fitCirc{1}(1),fitCirc{1}(2),fitCirc{1}(3),36,[0.6350,0.0780,0.1840],'filled');
quiver3(fitCirc{1}(1),fitCirc{1}(2),fitCirc{1}(3),fitCirc{4}(1),fitCirc{4}(2),fitCirc{4}(3), ...
    0.6*fitCirc{2},'filled','Color',[0.6350,0.0780,0.1840]);

startV = fitCirc{4};

% plot the plane of the fitting circle
plotOpts.FaceColor = [0.9290,0.6940,0.1250];
plotOpts.FaceAlpha = 0.1;
plotOpts.EdgeColor = 'none';
drawCirclePlane(fitCirc{1},fitCirc{2},fitCirc{5},plotOpts);

% plot the fitting circle
R = vecRot([0;0;1],fitCirc{5});
scaThe = linspace(0,2*pi);
scat(1,:) = fitCirc{2}*cos(scaThe);
scat(2,:) = fitCirc{2}*sin(scaThe);
scat(3,:) = zeros(1,length(scaThe));
circFit = R*scat + fitCirc{1};
plot3(circFit(1,:),circFit(2,:),circFit(3,:),'k--','LineWidth',1);

% plot the fitting arc
R = axesRot((rotz(0.5*fitCirc{3}))*[1;0;0],[0;0;1],fitCirc{4},fitCirc{5},'');
scaThe = linspace(0,fitCirc{3});
scat(1,:) = fitCirc{2}*cos(scaThe);
scat(2,:) = fitCirc{2}*sin(scaThe);
scat(3,:) = zeros(1,length(scaThe));
circFit = R*scat + fitCirc{1};
plot3(circFit(1,:),circFit(2,:),circFit(3,:), ...
    'Color',[0.8500,0.3250,0.0980],'LineWidth',3);
hold off;
grid on;

% set the x & y axis equal
xtick = get(gca,'XTick');
ytick = get(gca,'YTick');
ztick = get(gca,'ZTick');
% ytick间距，并将xtick间距设为与y相同
N = (max(xtick) - min(xtick))/(ztick(2)-ztick(1));
N_ = (max(ytick) - min(ytick)) / (ztick(2) - ztick(1));
xtick = xtick(1) + (0:(N  - 1))*(ztick(2)-ztick(1));
ytick = ytick(1) + (0:(N_ - 1)) * (ztick(2) - ztick(1));
set(gca,'XTick',xtick');%此时横轴和纵轴坐标刻度一致
set(gca,'YTick',ytick');

