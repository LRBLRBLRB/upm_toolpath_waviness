function [lineData,ang] = viewError(deltaZ,textFontSize,textFontType,unit)
%VIEWERROR 此处显示有关此函数的摘要
%   此处显示详细说明

% n = size(deltaZ);
ang = 0;
plotNum = 1000;
lineData = [];

fig = uifigure('Name','Check the useless extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [900,450];
figGridLayout = uigridlayout(fig,[3,6]);
figGridLayout.RowHeight = {'fit','fit','1x'};
figGridLayout.ColumnWidth = {'fit','4x','fit','fit','fit','4x'};

figNumLabel = uilabel(figGridLayout,'Text','Re-Plot Points No.', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figNumLabel.Layout.Row = 1;
figNumLabel.Layout.Column = [1,2];

figNumEdit = uieditfield(figGridLayout,'numeric','Limits',[50,2000], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'RoundFractionalValues','on','ValueDisplayFormat','%i');
figNumEdit.Value = plotNum;
figNumEdit.Layout.Row = 2;
figNumEdit.Layout.Column = 1;
figNumEdit.ValueChangedFcn = @(sld,event) figNumEditfieldUpdate(sld);

figNumSlider = uislider(figGridLayout,'Limits',[50,2000], ...
    'MajorTicks',[50,200,400,500,600,800,1000,1200,1400,1500,1600,1800,2000], ...
    'MajorTickLabels',{'50','','','500','','','1000','','','1500','','','2000'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figNumSlider.Value = plotNum;
figNumSlider.Layout.Row = 2;
figNumSlider.Layout.Column = 2;
figNumSlider.ValueChangedFcn = @(sld,event) figNumSliderUpdate(sld);

figRotLabel = uilabel(figGridLayout,'Text','Cross-Plane Angle', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figRotLabel.Layout.Row = 1;
figRotLabel.Layout.Column = [5,6];

figRotEdit = uieditfield(figGridLayout,'numeric','Limits',[0,360], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'ValueDisplayFormat','%.2f °');
figRotEdit.Value = ang;
figRotEdit.Layout.Row = 2;
figRotEdit.Layout.Column = 5;
figRotEdit.ValueChangedFcn = @(sld,event) figRotEditfieldUpdate(sld);

figRotSlider = uislider(figGridLayout,'Limits',[0,360], ...
    'MajorTicks',[0,30,60,90,120,150,180,210,240,270,300,330,360], ...
    'MajorTickLabels',{'0','','60','','120','','180','','240','','300','','°'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figRotSlider.Value = ang;
figRotSlider.Layout.Row = 2;
figRotSlider.Layout.Column = 6;
figRotSlider.ValueChangedFcn = @(sld,event) figRotSliderUpdate(sld);

figEnterButton = uibutton(figGridLayout,'push','Text','Enter', ...
    'FontName',textFontType,'FontSize',textFontSize);
figEnterButton.Layout.Row = 1;
figEnterButton.Layout.Column = [3,4];
figEnterButton.ButtonPushedFcn = @(sld,event) figEnterButtonPushed(fig);

figExportButton = uibutton(figGridLayout,'push','Text','Export', ...
    'FontName',textFontType,'FontSize',textFontSize);
figExportButton.Layout.Row = 2;
figExportButton.Layout.Column = [3,4];
figExportButton.ButtonPushedFcn = @(sld,event) figExportButtonPushed(fig);

fig3DAxes = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig3DAxes,'Surface Error');
fig3DAxes.Layout.Row = 3;
fig3DAxes.Layout.Column = [1,3];

surf(fig3DAxes,deltaZ(:,:,1),deltaZ(:,:,2),deltaZ(:,:,3),'EdgeColor','none');
colormap(fig3DAxes,turbo(256));
colorbar(fig3DAxes,'eastoutside');
clim(fig3DAxes,[min(deltaZ(:,:,3),[],'all'),max(deltaZ(:,:,3),[],'all')]);
xlabel(fig3DAxes,['x (',unit,')']);
ylabel(fig3DAxes,['y (',unit,')']);
zlabel(fig3DAxes,['\Deltaz (',unit,')']);
grid(fig3DAxes,'on');
hold(fig3DAxes,'off');
% axis(fig3DAxes,'equal');

fig2DAxes = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig2DAxes,'Surface Error in a Cross-Plane');
fig2DAxes.Layout.Row = 3;
fig2DAxes.Layout.Column = [4,6];

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        delete(app);
    end

    function figNumSliderUpdate(sld)
        plotNum = round(sld.Value);
        sld.Value = plotNum;
        figNumEdit.Value = plotNum;
    end

    function figNumEditfieldUpdate(sld)
        plotNum = sld.Value;
        figNumSlider.Value = sld.Value;
    end

    function figRotSliderUpdate(sld)
        ang = sld.Value;
        figRotEdit.Value = sld.Value;
    end

    function figRotEditfieldUpdate(sld)
        ang = sld.Value;
        figRotSlider.Value = sld.Value;
    end

    % Create ValueChangedFcn callback
    function figEnterButtonPushed(app)
        figWaitbar = uiprogressdlg(app,'Title','Line Data Extraction', ...
            'Cancelable','on','Icon','success','Indeterminate','on');
        figWaitbar.Message = 'Data Extracting ... ';

        % rotation transform
        rot = rotz(ang);
        surfData = reshape(deltaZ,[],3);
        surfData = (rot*surfData')';
        minX = min(surfData(:,1));
        maxX = max(surfData(:,1));
        minY = min(surfData(:,2));
        maxY = max(surfData(:,2));
        minZ = min(surfData(:,3));
        maxZ = max(surfData(:,3));
        if figWaitbar.CancelRequested, return; end

        % same-XY elimination
        surfData = surfData(all(~isnan(surfData),2),:);
        [surfDataUniZ,surfDataUniXY] = groupsummary( ...
            surfData(:,3),surfData(:,1:2),@mean);
        if figWaitbar.CancelRequested, return; end
        xPlot = linspace(minX,maxX,plotNum);
        yPlot = linspace(minY,maxY,plotNum);
        [deltaZUni(:,:,1),deltaZUni(:,:,2)] = meshgrid(xPlot,yPlot);
        deltaZUni(:,:,3) = griddata(surfDataUniXY{1},surfDataUniXY{2}, ...
            surfDataUniZ,deltaZUni(:,:,1),deltaZUni(:,:,2));
        if figWaitbar.CancelRequested, return; end

        % line data
        clear('lineData');
        [lineData(:,1),Yq(:,1)] = meshgrid((minX:0.1:maxX)',0);
        lineData(:,2) = interp2(deltaZUni(:,:,1),deltaZUni(:,:,2), ...
            deltaZUni(:,:,3),lineData(:,1),Yq);
        %bsplineSurfPts_spapi(surfData,3,3,u,v)
        hold(fig2DAxes,'off');
        plot(fig2DAxes,lineData(:,1),lineData(:,2),'.-');
        xlabel(fig2DAxes,['x (',unit,')']);
        ylabel(fig2DAxes,['y (',unit,')']);
        grid(fig2DAxes,'on');

        % surf data
        surf(fig3DAxes,deltaZUni(:,:,1),deltaZUni(:,:,2),deltaZUni(:,:,3), ...
            'EdgeColor','none');
        hold(fig3DAxes,'on');
        fill3(fig3DAxes,[1.2*minX;1.2*maxX;1.2*maxX;1.2*minX],[0;0;0;0], ...
            [1.2*minZ;1.2*minZ;1.2*maxZ;1.2*maxZ],[0.6350 0.0780 0.1840], ...
            'EdgeColor','none','FaceAlpha',0.3);
        line(fig3DAxes,[1.2*minX;1.2*maxX],[0;0],[1.2*maxZ;1.2*maxZ], ...
            'Color',[0.6350 0.0780 0.1840],'LineWidth',3, ...
            'Marker','.','MarkerSize',18);
        hold(fig3DAxes,'off');

        close(figWaitbar);
    end

    function figExportButtonPushed(app)
        uiresume(app);
    end
end