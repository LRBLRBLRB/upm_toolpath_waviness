function [ang,varargout] = viewError(deltaZ,textFontSize,textFontType,unit)
%VIEWERROR 此处显示有关此函数的摘要
%   此处显示详细说明

% n = size(deltaZ);
ang = 0;
plotNum = 1000;
lineData = [];
zernike2 = [];
zernike3 = [];
zOrder = 4;
dataName = [];

maxR = max(sqrt(deltaZ(:,:,1).^2 + deltaZ(:,:,2).^2),[],'all');
minZ = min(deltaZ(:,:,3),[],'all');
maxZ = max(deltaZ(:,:,3),[],'all');

fig = uifigure('Name','Check the useless extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(1:2) = fig.Position(1:2) - [900 - 560,600 - 420];
fig.Position(3:4) = [900,600];
figGridLayout = uigridlayout(fig,[4,6]);
figGridLayout.RowHeight = {'fit','fit','2x','1x'};
figGridLayout.ColumnWidth = {'fit','fit','1x','fit','fit','1x'};

figNumLabel = uilabel(figGridLayout,'Text','Re-Plot Points No.', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figNumLabel.Layout.Row = 1;
figNumLabel.Layout.Column = 1;

figNumEdit = uieditfield(figGridLayout,'numeric','Limits',[50,2000], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'RoundFractionalValues','on','ValueDisplayFormat','%i');
figNumEdit.Value = plotNum;
figNumEdit.Layout.Row = 1;
figNumEdit.Layout.Column = 2;
figNumEdit.ValueChangedFcn = @(sld,event) figNumEditfieldUpdate(sld);

figNumSlider = uislider(figGridLayout,'Limits',[50,2000], ...
    'MajorTicks',[50,200,400,500,600,800,1000,1200,1400,1500,1600,1800,2000], ...
    'MajorTickLabels',{'50','','','500','','','1000','','','1500','','','2000'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figNumSlider.Value = plotNum;
figNumSlider.Layout.Row = 1;
figNumSlider.Layout.Column = 3;
figNumSlider.ValueChangedFcn = @(sld,event) figNumSliderUpdate(sld);

figRotLabel = uilabel(figGridLayout,'Text','Cross-Plane Angle', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figRotLabel.Layout.Row = 2;
figRotLabel.Layout.Column = 1;

figRotEdit = uieditfield(figGridLayout,'numeric','Limits',[0,360], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'ValueDisplayFormat','%.2f °');
figRotEdit.Value = ang;
figRotEdit.Layout.Row = 2;
figRotEdit.Layout.Column = 2;
figRotEdit.ValueChangedFcn = @(sld,event) figRotEditfieldUpdate(sld);

figRotSlider = uislider(figGridLayout,'Limits',[0,360], ...
    'MajorTicks',[0,30,60,90,120,150,180,210,240,270,300,330,360], ...
    'MajorTickLabels',{'0','','60','','120','','180','','240','','300','','°'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figRotSlider.Value = ang;
figRotSlider.Layout.Row = 2;
figRotSlider.Layout.Column = 3;
figRotSlider.ValueChangedFcn = @(sld,event) figRotSliderUpdate(sld);

figOrderLabel = uilabel(figGridLayout,'Text','Zernike Order', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figOrderLabel.Layout.Row = 1;
figOrderLabel.Layout.Column = 4;

figOrderEdit = uieditfield(figGridLayout,'numeric','Limits',[4,80], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'RoundFractionalValues','on','ValueDisplayFormat','%i','Editable','off');
figOrderEdit.Value = zOrder;
figOrderEdit.Layout.Row = 1;
figOrderEdit.Layout.Column = 5;
figOrderEdit.ValueChangedFcn = @(sld,event) figOrderEditfieldUpdate(sld);

figOrderSlider = uislider(figGridLayout,'Limits',[4,80], ...
    'MajorTicks',[4,10,20,30,40,50,60,70,80], ...
    'MajorTickLabels',{'4','','20','','','50','','','80'}, ...
    'FontName',textFontType,'FontSize',textFontSize,'Enable','off');
figOrderSlider.Value = zOrder;
figOrderSlider.Layout.Row = 1;
figOrderSlider.Layout.Column = 6;
figOrderSlider.ValueChangedFcn = @(sld,event) figOrderSliderUpdate(sld);

figExtractOriButton = uibutton(figGridLayout,'push','Text','Extract Original', ...
    'FontName',textFontType,'FontSize',textFontSize);
figExtractOriButton.Layout.Row = 2;
figExtractOriButton.Layout.Column = 4;
figExtractOriButton.ButtonPushedFcn = @(sld,event) figExtractOriButtonPushed(fig);

figExtractZernikeButton = uibutton(figGridLayout,'push','Text','Extract Zernike', ...
    'FontName',textFontType,'FontSize',textFontSize);
figExtractZernikeButton.Layout.Row = 2;
figExtractZernikeButton.Layout.Column = 5;
figExtractZernikeButton.ButtonPushedFcn = @(sld,event) figExtractZernikeButtonPushed(fig);

figExportButton = uibutton(figGridLayout,'push','Text','Export', ...
    'FontName',textFontType,'FontSize',textFontSize);
figExportButton.Layout.Row = 2;
figExportButton.Layout.Column = 6;
figExportButton.ButtonPushedFcn = @(sld,event) figExportButtonPushed(fig);

fig3DAxes = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig3DAxes,'Surface Error');
fig3DAxes.Layout.Row = 3;
fig3DAxes.Layout.Column = [1,3];
surf(fig3DAxes,deltaZ(:,:,1),deltaZ(:,:,2),deltaZ(:,:,3),'EdgeColor','none');
hold(fig3DAxes,'on');
planeFill = [];
planeLine = [];
[~,planeFill,planeLine] = surfRot(fig3DAxes,ang,planeFill,planeLine);
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
fig2DAxes.Layout.Row = 4;
fig2DAxes.Layout.Column = [1,3];

fig3DAxesZernike = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig3DAxesZernike,'Surface Error');
fig3DAxesZernike.Layout.Row = 3;
fig3DAxesZernike.Layout.Column = [4,6];
planeFill2 = [];
planeLine2 = [];

fig2DAxesZernike = uiaxes(figGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig2DAxesZernike,'Surface Error in a Cross-Plane');
fig2DAxesZernike.Layout.Row = 4;
fig2DAxesZernike.Layout.Column = [4,6];

fig.Visible = 'on';

uiwait(fig);

switch nargout
    case 2
        varargout{1} = lineData;
        if strcmp(dataName,filesep)
            save(dataName,"lineData");
        end
    case 4
        varargout{1} = lineData;
        varargout{2} = zernike3;
        varargout{3} = zernike2;
        if strcmp(dataName,filesep)
            save(dataName,"lineData","zernike2","zernike3");
        end
end

    function [rot,planeFill,planeLine] = surfRot(fig,ang,planeFill,planeLine)
        % rotation of the plane
        rot = rotz(ang);
        planePt = [-1.2*maxR,1.2*maxR,1.2*maxR,-1.2*maxR; ...
            0,0,0,0; ...
            1.2*minZ,1.2*minZ,1.2*maxZ,1.2*maxZ];
        planePt = rot*planePt;
        hold(fig,'on');
        if exist('planeFill','var')
            delete(planeFill);
            delete(planeLine);
        end
        planeFill = fill3(fig,planePt(1,:),planePt(2,:),planePt(3,:), ...
            [0.6350 0.0780 0.1840],'EdgeColor','none','FaceAlpha',0.3);
        planeLine = line(fig,planePt(1,3:4),planePt(2,3:4),planePt(3,3:4), ...
            'Color',[0.6350 0.0780 0.1840],'LineWidth',3, ...
            'Marker','.','MarkerSize',18);
        hold(fig,'off');
    end

    function extractPlot(surfData,fig,figWaitbar)
        minX = min(surfData(:,1));
        maxX = max(surfData(:,1));
        minY = min(surfData(:,2));
        maxY = max(surfData(:,2));

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
        hold(fig,'off');
        plot(fig,lineData(:,1),lineData(:,2),'.-');
        xlabel(fig,['x (',unit,')']);
        ylabel(fig,['y (',unit,')']);
        grid(fig,'on');
    end

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
        [rot,planeFill,planeLine] = surfRot(fig3DAxes,ang,planeFill,planeLine);
    end

    function figRotEditfieldUpdate(sld)
        ang = sld.Value;
        figRotSlider.Value = sld.Value;
        [rot,planeFill,planeLine] = surfRot(fig3DAxes,ang,planeFill,planeLine);
    end

    function figOrderEditfieldUpdate(sld)
        zOrder = sld.Value;
        figOrderSlider.Value = sld.Value;
    end

    function figOrderSliderUpdate(sld)
        zOrder = sld.Value;
        figOrderEdit.Value = sld.Value;
        close(figWaitbar);
    end

    % Create ValueChangedFcn callback
    function figExtractOriButtonPushed(app)
        figWaitbar = uiprogressdlg(app,'Title','Line Data Extraction', ...
            'Cancelable','on','Icon','success','Indeterminate','on');
        figWaitbar.Message = 'Data Extracting ... ';

        % rotation transform
        [rot,planeFill,planeLine] = surfRot(fig3DAxes,ang,planeFill,planeLine);
        surfData = reshape(deltaZ,[],3);
        surfData = (rot*surfData')';
        if figWaitbar.CancelRequested, return; end

        extractPlot(surfData,fig2DAxes,figWaitbar);
        % if figWaitbar.CancelRequested, return; end

        close(figWaitbar);
    end

    function figExtractZernikeButtonPushed(app)
        figOrderSlider.Enable = 'on';
        figOrderEdit.Editable = 'on';
        figWaitbar = uiprogressdlg(app,'Title','Line Data Extraction', ...
            'Cancelable','on','Icon','success','Indeterminate','on');
        figWaitbar.Message = 'Data Extracting ... ';

        % rotation transform
        [rot,planeFill2,planeLine2] = surfRot(fig3DAxesZernike,ang,planeFill2,planeLine2);
        if figWaitbar.CancelRequested, return; end
        surfData = reshape(deltaZ,[],3);
        % surfData = (rot*surfData')';

        zernikeZ = zernikeProcess(deltaZ(:,:,3),zOrder);
        if figWaitbar.CancelRequested, return; end
        surf(fig3DAxesZernike,deltaZ(:,:,1),deltaZ(:,:,2),zernikeZ,'EdgeColor','none');
        hold(fig3DAxesZernike,'on');
        colormap(fig3DAxesZernike,turbo(256));
        colorbar(fig3DAxesZernike,'eastoutside');
        clim(fig3DAxesZernike,[min(deltaZ(:,:,3),[],'all'),max(deltaZ(:,:,3),[],'all')]);
        xlabel(fig3DAxesZernike,['x (',unit,')']);
        ylabel(fig3DAxesZernike,['y (',unit,')']);
        zlabel(fig3DAxesZernike,['\Deltaz (',unit,')']);
        grid(fig3DAxesZernike,'on');

        if figWaitbar.CancelRequested, return; end

        extractPlot(surfData,fig2DAxesZernike,figWaitbar);
        % if figWaitbar.CancelRequested, return; end

        close(figWaitbar);
    end

    function figExportButtonPushed(app)
        [dataFileName,dataDirName] = uiputfile({ ...
            '*.mat','MAT-file(*.mat)'; ...
            '*.*','all file(*.*)';...
            }, ...
            'Select a file to save the data processing results');
        dataName = fullfile(dataDirName,dataFileName);
        uiresume(app);
    end
end