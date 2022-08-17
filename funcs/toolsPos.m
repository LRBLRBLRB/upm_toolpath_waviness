%% 相邻刀位点的残高
x1 = [0;0]; vec1 = [0;1];
x2 = [r0;r0/30]; vec2 = [cosd(90+8);sind(90+8)];
R = rotz(5,'deg'); % rotation matrix about z axis
tool1 = toolPt(1:2:end,:) + x1;
tool2 = R(1:2,1:2)*toolPt(1:2:end,:) + x2;

[res,interPt] = residualHigh(x1,vec1,tool1,x2,vec2,tool2);

figure('Name','Residual of the adjacent tool');
plot(tool1(1,:),tool1(2,:),'Color',[0,0.45,0.74]); hold on;
plot(tool2(1,:),tool2(2,:),'Color',[0.85,0.33,0.10]);
plot(interPt(1),interPt(2),'*','Color',[0.49,0.18,0.56]);
plot(x1(1)+radius*vec1(1),x1(2)+radius*vec1(2),'o','Color',[0,0.45,0.74],'MarkerFaceColor',[0,0.45,0.74]);
quiver(x1(1),x1(2),radius*vec1(1),radius*vec1(2),'LineWidth',2,'Color',[0,0.45,0.74],'AutoScale','off');
plot(x2(1)+radius*vec2(1),x2(2)+radius*vec2(2),'o','Color',[0.85,0.33,0.10],'MarkerFaceColor',[0.85,0.33,0.10]);
quiver(x2(1),x2(2),radius*vec2(1),radius*vec2(2),'LineWidth',2,'Color',[0.85,0.33,0.10],'AutoScale','off');
axis equal