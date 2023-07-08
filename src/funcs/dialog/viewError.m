function [ang,varargout] = viewError(dataOri,textFontSize,textFontType,unit)
%VIEWERROR an app to view the error and extract the high-order error

ang = 0;
plotNum = 1000;
zOrder = 4;
dataName = [];
dataLineOri = [];
dataZernike = dataOri;
dataLineZernike = [];
dataError = dataOri;
dataLineError = [];

maxR = max(sqrt(dataOri(:,:,1).^2 + dataOri(:,:,2).^2),[],'all');
minZ = min(dataOri(:,:,3),[],'all');
maxZ = max(dataOri(:,:,3),[],'all');

fig = uifigure('Name','Check the useless extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(1:2) = fig.Position(1:2) - [1000 - fig.Position(3),900 - fig.Position(4)];
fig.Position(3:4) = [1000,900];
figGridLayout = uigridlayout(fig,[1,2]);
figGridLayout.RowHeight = {'2x','1x','1x','2x'};
figGridLayout.ColumnWidth = {'1x','1x'};

figParamPanel = uipanel(figGridLayout,'Title','Parameters','TitlePosition','lefttop', ...
    'FontName',textFontType,'FontSize',textFontSize);
figParamPanel.Layout.Row = 1;
figParamPanel.Layout.Column = 1;
figParamPanelGridLayout = uigridlayout(figParamPanel,[4,3]);
figParamPanelGridLayout.RowHeight = {'fit','fit','fit','1x'};
figParamPanelGridLayout.ColumnWidth = {'fit','fit','1x'};

figNumLabel = uilabel(figParamPanelGridLayout,'Text','Re-Plot Points No.', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figNumLabel.Layout.Row = 1;
figNumLabel.Layout.Column = 1;

figNumEdit = uieditfield(figParamPanelGridLayout,'numeric','Limits',[50,2000], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'RoundFractionalValues','on','ValueDisplayFormat','%i');
figNumEdit.Value = plotNum;
figNumEdit.Layout.Row = 1;
figNumEdit.Layout.Column = 2;
figNumEdit.ValueChangedFcn = @(sld,event) figNumEditfieldUpdate(sld);

figNumSlider = uislider(figParamPanelGridLayout,'Limits',[50,2000], ...
    'MajorTicks',[50,200,400,500,600,800,1000,1200,1400,1500,1600,1800,2000], ...
    'MajorTickLabels',{'50','','','500','','','1000','','','1500','','','2000'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figNumSlider.Value = plotNum;
figNumSlider.Layout.Row = 1;
figNumSlider.Layout.Column = 3;
figNumSlider.ValueChangedFcn = @(sld,event) figNumSliderUpdate(sld);

figRotLabel = uilabel(figParamPanelGridLayout,'Text','Cross-Plane Angle', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figRotLabel.Layout.Row = 2;
figRotLabel.Layout.Column = 1;

figRotEdit = uieditfield(figParamPanelGridLayout,'numeric','Limits',[0,360], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'ValueDisplayFormat','%.2f °');
figRotEdit.Value = ang;
figRotEdit.Layout.Row = 2;
figRotEdit.Layout.Column = 2;
figRotEdit.ValueChangedFcn = @(sld,event) figRotEditfieldUpdate(sld);

figRotSlider = uislider(figParamPanelGridLayout,'Limits',[0,360], ...
    'MajorTicks',[0,30,60,90,120,150,180,210,240,270,300,330,360], ...
    'MajorTickLabels',{'0','','60','','120','','180','','240','','300','','°'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
figRotSlider.Value = ang;
figRotSlider.Layout.Row = 2;
figRotSlider.Layout.Column = 3;
figRotSlider.ValueChangedFcn = @(sld,event) figRotSliderUpdate(sld);
% figRotSlider.ValueChangingFcn = @(sld,event) figRotSliderUpdate(sld);

figOrderLabel = uilabel(figParamPanelGridLayout,'Text','Zernike Order', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'HorizontalAlignment','center');
figOrderLabel.Layout.Row = 3;
figOrderLabel.Layout.Column = 1;

figOrderEdit = uieditfield(figParamPanelGridLayout,'numeric','Limits',[0,20], ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'RoundFractionalValues','on','ValueDisplayFormat','%i','Editable','off');
figOrderEdit.Value = zOrder;
figOrderEdit.Layout.Row = 3;
figOrderEdit.Layout.Column = 2;
figOrderEdit.ValueChangedFcn = @(sld,event) figOrderEditfieldUpdate(sld);

figOrderSlider = uislider(figParamPanelGridLayout,'Limits',[0,20], ...
    'MajorTicks',[0,5,10,15,20], ...
    'MajorTickLabels',{'0','5','10','15','20'}, ...
    'FontName',textFontType,'FontSize',textFontSize,'Enable','off');
figOrderSlider.Value = zOrder;
figOrderSlider.Layout.Row = 3;
figOrderSlider.Layout.Column = 3;
figOrderSlider.ValueChangedFcn = @(sld,event) figOrderSliderUpdate(sld);

figExtract2DButton = uibutton(figParamPanelGridLayout,'push','Text','2D Data Extraction', ...
    'FontName',textFontType,'FontSize',textFontSize,'WordWrap','on');
figExtract2DButton.Layout.Row = 4;
figExtract2DButton.Layout.Column = 2;
figExtract2DButton.ButtonPushedFcn = @(sld,event) figExtract2DButtonPushed(fig);

figFitZernikeButton = uibutton(figParamPanelGridLayout,'push','Text','Zernike Fitting', ...
    'FontName',textFontType,'FontSize',textFontSize,'WordWrap','on');
figFitZernikeButton.Layout.Row = 4;
figFitZernikeButton.Layout.Column = 1;
figFitZernikeButton.ButtonPushedFcn = @(sld,event) figFitZernikeButtonPushed(fig);

figExportButton = uibutton(figParamPanelGridLayout,'push','Text','Data Exporting', ...
    'FontName',textFontType,'FontSize',textFontSize);
figExportButton.Layout.Row = 4;
figExportButton.Layout.Column = 3;
figExportButton.ButtonPushedFcn = @(sld,event) figExportButtonPushed(fig);

figOriPanel = uipanel(figGridLayout,'Title','Origin','TitlePosition','lefttop', ...
    'FontName',textFontType,'FontSize',textFontSize);
figOriPanel.Layout.Row = [2,4];
figOriPanel.Layout.Column = 1;
figOriPanelGridLayout = uigridlayout(figOriPanel,[2,1], ...
    'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',0);
figOriPanelGridLayout.RowHeight = {'2x','1x'};

fig3DAxesOri = uiaxes(figOriPanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
fig3DAxesOri.Layout.Row = 1;
fig3DAxesOri.Layout.Column = 1;
title(fig3DAxesOri,'Surface Error');
xlim(fig3DAxesOri,[-1.2*maxR,1.2*maxR]);
ylim(fig3DAxesOri,[-1.2*maxR,1.2*maxR]);
zlim(fig3DAxesOri,[1.5*minZ,1.5*maxZ]);
surfPlot(fig3DAxesOri,dataOri);
hold(fig3DAxesOri,'on');
[planeFillOri,planeLineOri] = planeRot(fig3DAxesOri,ang,[],[]);
hold(fig3DAxesOri,'off');

fig2DAxesOri = uiaxes(figOriPanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
title(fig2DAxesOri,'Surface Error in a Cross-Plane');
fig2DAxesOri.Layout.Row = 2;
fig2DAxesOri.Layout.Column = 1;

figZernikePanel = uipanel(figGridLayout,'Title','Zernike', ...
    'TitlePosition','righttop', ...
    'FontName',textFontType,'FontSize',textFontSize);
figZernikePanel.Layout.Row = [1,2];
figZernikePanel.Layout.Column = 2;
figZernikePanelGridLayout = uigridlayout(figZernikePanel,[2,1], ...
    'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',0);
figZernikePanelGridLayout.RowHeight = {'2x','1x'};

fig3DAxesZernike = uiaxes(figZernikePanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
fig3DAxesZernike.Layout.Row = 1;
fig3DAxesZernike.Layout.Column = 1;
xlim(fig3DAxesZernike,[-1.2*maxR,1.2*maxR]);
ylim(fig3DAxesZernike,[-1.2*maxR,1.2*maxR]);
zlim(fig3DAxesZernike,[1.5*minZ,1.5*maxZ]);
planeFillZernike = [];
planeLineZernike = [];

fig2DAxesZernike = uiaxes(figZernikePanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
fig2DAxesZernike.Layout.Row = 2;
fig2DAxesZernike.Layout.Column = 1;

figErrorPanel = uipanel(figGridLayout,'Title','Error', ...
    'TitlePosition','righttop', ...
    'FontName',textFontType,'FontSize',textFontSize);
figErrorPanel.Layout.Row = [3,4];
figErrorPanel.Layout.Column = 2;
figErrorPanelGridLayout = uigridlayout(figErrorPanel,[2,1], ...
    'Padding',[0 0 0 0],'ColumnSpacing',0,'RowSpacing',0);
figErrorPanelGridLayout.RowHeight = {'2x','1x'};

fig3DAxesError = uiaxes(figErrorPanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
fig3DAxesError.Layout.Row = 1;
fig3DAxesError.Layout.Column = 1;
xlim(fig3DAxesError,[-1.2*maxR,1.2*maxR]);
ylim(fig3DAxesError,[-1.2*maxR,1.2*maxR]);
zlim(fig3DAxesError,[1.5*minZ,1.5*maxZ]);
planeFillError = [];
planeLineError = [];

fig2DAxesError = uiaxes(figErrorPanelGridLayout, ...
    'FontName',textFontType,'FontSize',textFontSize);
fig2DAxesError.Layout.Row = 2;
fig2DAxesError.Layout.Column = 1;

fig.Visible = 'on';

uiwait(fig);
fig.WindowStyle = 'normal';

switch nargout
    case 2
        varargout{1} = dataLineOri;
        if strcmp(dataName,filesep)
            save(dataName,"dataLineOri");
        end
    case 4
        varargout{1} = dataLineOri;
        varargout{2} = dataZernike;
        varargout{3} = dataLineZernike;
        if strcmp(dataName,filesep)
            save(dataName,"dataLineOri","dataZernike","dataLineZernike");
        end
    case 6
        varargout{1} = dataLineOri;
        varargout{2} = dataZernike;
        varargout{3} = dataLineZernike;
        varargout{4} = dataError;
        varargout{5} = dataLineError;
        if strcmp(dataName,filesep)
            save(dataName,"dataLineOri","dataZernike","dataLineZernike", ...
                "dataError","dataLineError");
        end
end

    function surfPlot(fig,data)
        cla(fig);
        surf(fig,data(:,:,1),data(:,:,2),data(:,:,3), ...
            'EdgeColor','none');
        hold(fig,'on');
        colormap(fig,turbo(256));
        colorbar(fig,'eastoutside');
        minData = min(data(:,:,3),[],'all');
        maxData = max(data(:,:,3),[],'all');
        if minData ~= maxData
            clim(fig,[minData,maxData]);
        end
        xlabel(fig,['x (',unit,')']);
        ylabel(fig,['y (',unit,')']);
        zlabel(fig,['\Deltaz (',unit,')']);
        grid(fig,'on');
    end

    function [planeFill,planeLine] = planeRot(fig,ang,planeFill,planeLine)
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

    function [lineData,figWaitbar] = extractPlot(surfData,fig,figWaitbar)
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
        [lineData(:,1),Yq(:,1)] = meshgrid((minX:0.1:maxX)',0);
        lineData(:,2) = interp2(deltaZUni(:,:,1),deltaZUni(:,:,2), ...
            deltaZUni(:,:,3),lineData(:,1),Yq);
        %bsplineSurfPts_spapi(surfData,3,3,u,v)
        hold(fig,'off');
        plot(fig,lineData(:,1),lineData(:,2),'-','LineWidth',0.5);
        xlabel(fig,['x (',unit,')']);
        ylabel(fig,['y (',unit,')']);
        grid(fig,'on');
    end

    function fig_close_req(sld)
        delete(sld);
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
        [planeFillOri,planeLineOri] = planeRot( ...
            fig3DAxesOri,ang,planeFillOri,planeLineOri);
        cla(fig2DAxesOri);
        [planeFillZernike,planeLineZernike] = planeRot( ...
            fig3DAxesZernike,ang,planeFillZernike,planeLineZernike);
        cla(fig2DAxesZernike);
        [planeFillError,planeLineError] = planeRot( ...
            fig3DAxesError,ang,planeFillError,planeLineError);
        cla(fig2DAxesError);
    end

    function figRotEditfieldUpdate(sld)
        ang = sld.Value;
        figRotSlider.Value = sld.Value;
        [planeFillOri,planeLineOri] = planeRot( ...
            fig3DAxesOri,ang,planeFillOri,planeLineOri);
    end

    function figOrderEditfieldUpdate(sld)
        zOrder = sld.Value;
        figOrderSlider.Value = sld.Value;
    end

    function figOrderSliderUpdate(sld)
        zOrder = round(sld.Value);
        figOrderEdit.Value = sld.Value;
    end

    function figFitZernikeButtonPushed(app)
        figOrderSlider.Enable = 'on';
        figOrderEdit.Editable = 'on';
        figWaitbar = uiprogressdlg(app,'Title','Zernike Fitting', ...
            'Cancelable','on','Icon','success','Indeterminate','on');
        figWaitbar.Message = 'Data Fitting ... ';

        % plot the Zernike result
        [planeFillOri,planeLineOri] = planeRot( ...
            fig3DAxesOri,ang,planeFillOri,planeLineOri);

        % plot the Zernike result
        dataZernike(:,:,3) = zernikeProcess(dataOri(:,:,3),zOrder);
        surfPlot(fig3DAxesZernike,dataZernike);
        if figWaitbar.CancelRequested, return; end

        [planeFillZernike,planeLineZernike] = planeRot( ...
            fig3DAxesZernike,ang,planeFillZernike,planeLineZernike);
        if figWaitbar.CancelRequested, return; end

        % plot the error result
        dataError(:,:,3) = dataOri(:,:,3) - dataZernike(:,:,3);
        surfPlot(fig3DAxesError,dataError);
        if figWaitbar.CancelRequested, return; end

        [planeFillError,planeLineError] = planeRot( ...
            fig3DAxesError,ang,planeFillError,planeLineError);
        if figWaitbar.CancelRequested, return; end

        close(figWaitbar);
    end

    function figExtract2DButtonPushed(app)
        figWaitbar = uiprogressdlg(app,'Title','Line Data Extraction', ...
            'Cancelable','on','Icon','success','Indeterminate','on');
        figWaitbar.Message = 'Data Extracting ... ';

        % rotation transform
        rot = rotz(ang);
        if figWaitbar.CancelRequested, return; end

        % original data
        surfData = reshape(dataOri,[],3);
        surfData = (rot*surfData')';
        [dataLineOri,figWaitbar] = extractPlot(surfData,fig2DAxesOri,figWaitbar);
        if figWaitbar.CancelRequested, return; end

        % Zernike data
        surfData = reshape(dataZernike,[],3);
        surfData = (rot*surfData')';
        [dataLineZernike,figWaitbar] = extractPlot(surfData,fig2DAxesZernike,figWaitbar);
        if figWaitbar.CancelRequested, return; end

        % error data
        surfData = reshape(dataError,[],3);
        surfData = (rot*surfData')';
        [dataLineError,figWaitbar] = extractPlot(surfData,fig2DAxesError,figWaitbar);
        if figWaitbar.CancelRequested, return; end

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