% continue to extract the residual error from the whole error in both 2D
% and 3D cases

zernikeZ = zernikeProcess(deltaZ(:,:,3),3);

% rot = rotz(lineAng);
% lineData

figure('Name','2D error extraction');
tiledlayout(2,1);
nexttile(1);
surf(deltaZ(:,:,1),deltaZ(:,:,2),deltaZ(:,:,3) - zernikeZ,'EdgeColor','none');
axis equal;
hold('on');
colormap(turbo(256));
colorbar('eastoutside');
% clim([min(deltaZ(:,:,3),[],'all'),max(deltaZ(:,:,3),[],'all')]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['\Deltaz (',unit,')']);
% 
% nexttile(2);
% 
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% grid on;
% xlabel(['r (',unit,')']);
% ylabel(['\Deltaz (',unit,')']);

%% 3D cases
% deltaZ