function viewError(deltaZ,textFontSize,textFontType,unit)
%VIEWERROR 此处显示有关此函数的摘要
%   此处显示详细说明

fig = uifigure('Name','Check the useless extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [900,500];
figGridLayout = uigridlayout(fig,[3,3]);
figGridLayout.RowHeight = {'fit','fit','1x'};
figGridLayout.ColumnWidth = {'2x','fit','1x'};

fig3DAxes = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig3DAxes,'Surface Error');
fig3DAxes.Layout.Row = [1,3];
fig3DAxes.Layout.Column = 1;

surf(fig3DAxes,deltaZ(:,:,1),deltaZ(:,:,2),deltaZ(:,:,3),'EdgeColor','none');
colormap(fig3DAxes,turbo(256));
colorbar(fig3DAxes,'eastoutside');
clim(fig3DAxes,[min(deltaZ(:,:,3),[],'all'),max(deltaZ(:,:,3),[],'all')]);
xlabel(fig3DAxes,['x (',unit,')']);
ylabel(fig3DAxes,['y (',unit,')']);
zlabel(fig3DAxes,['\Deltaz (',unit,')']);
grid(fig3DAxes,'on');
% axis(fig3DAxes,'equal');

figLabel = uilabel(figGridLayout,'Text','Cross-Plane Angle', ...
    'FontName',textFontType,'FontSize',textFontSize);
figLabel.Layout.Row = 1;
figLabel.Layout.Column = [2,3];

figEdit = uieditfield(figGridLayout,'numeric','Value',0,'Limits',[0,180], ...
    'FontName',textFontType,'FontSize',textFontSize);
figEdit.Layout.Row = 2;
figEdit.Layout.Column = 2;
figEdit.ValueChangedFcn = @(sld,event) editfieldUpdate(sld);

figSlider = uislider(figGridLayout,'Value',0,'Limits',[0,180], ...
    'MajorTicks',[0,30,60,90,120,150,180], ...
    'MajorTickLabels',{'0','30','60','90','120','150','°'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figSlider.Layout.Row = 2;
figSlider.Layout.Column = 3;
figSlider.ValueChangedFcn = @(sld,event) sliderUpdate(sld);

fig2DAxes = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig2DAxes,'Surface Error in a Cross-Plane');
fig2DAxes.Layout.Row = 3;
fig2DAxes.Layout.Column = [2,3];

n = size(deltaZ);
minX = min(deltaZ(:,:,1),[],'all');
maxX = max(deltaZ(:,:,1),[],'all');

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        delete(app);
    end

    function sliderUpdate(sld)
        ang = sld.Value;
        figEdit.Value = sld.Value;
        updateCrossplane(ang);
    end

    function editfieldUpdate(sld)
        ang = sld.Value;
        figSlider.Value = sld.Value;
        updateCrossplane(ang);
    end

    % Create ValueChangedFcn callback
    function updateCrossplane(ang)
        rot = rotz(ang);
        surfData = reshape(deltaZ,[],3);
        surfData = (rot*surfData')';
        surfData = reshape(surfData,n);
        [Xq,Yq] = meshgrid(0.5*minX:1:0.5*maxX,-100:1:100);
        [resUnique,peakPtUnique] = groupsummary(surfData,surfData(1:2,:)',@max);
        lineData = interp2(surfData(:,:,1),surfData(:,:,2),surfData(:,:,3), ...
            Xq,Yq);
        %bsplineSurfPts_spapi(surfData,3,3,u,v)
        hold(fig2DAxes,'off');
        plot(fig2DAxes,lineData(:,1),lineData(:,3));
        xlabel(fig2DAxes,['x (',unit,')']);
        ylabel(fig2DAxes,['y (',unit,')']);
    end
end