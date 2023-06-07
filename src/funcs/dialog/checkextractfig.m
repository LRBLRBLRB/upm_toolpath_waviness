function [isContinue,lineFitMaxDist] = checkextractfig(lineFitMaxDist,figShow)

fig = uifigure('Name','Check the line extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);

fig.Position(3:4) = [300,150];
figGridLayout = uigridlayout(fig,[3,2]);
figGridLayout.RowHeight = {'1x','fit','fit'};
figGridLayout.ColumnWidth = {'1x','1x'};

figLabel = uilabel(figGridLayout,'Text','Line extraction is finished.', ...
    'HorizontalAlignment','center','FontSize',14);
figLabel.Layout.Row = 1;
figLabel.Layout.Column = [1,2];

figLineMaxDistLabel = uilabel(figGridLayout,'Text','Ransac max distance:');
figLineMaxDistLabel.Layout.Row = 2;
figLineMaxDistLabel.Layout.Column = 1;

figLineMaxDistEditfield = uieditfield(figGridLayout,'numeric','Limits',[0 inf], ...
    'HorizontalAlignment','center','ValueDisplayFormat','%d', ...
    'Value',lineFitMaxDist);
figLineMaxDistEditfield.Layout.Row = 2;
figLineMaxDistEditfield.Layout.Column = 2;

figEnterButton = uibutton(figGridLayout,'push','Text','Enter','Visible','on');
figEnterButton.Layout.Row = 3;
figEnterButton.Layout.Column = 1;
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figReFitButton = uibutton(figGridLayout,'push','Text','Re-fit','Visible','on');
figReFitButton.Layout.Row = 3;
figReFitButton.Layout.Column = 2;
figReFitButton.ButtonPushedFcn = @(app,event)fig_refit_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        isContinue = 0;
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        fig_close_req(fig);
    end

    function fig_refit_button_pushed
        isContinue = 0;
        lineFitMaxDist = figLineMaxDistEditfield.Value;
        delete(figShow);
        fig_close_req(fig);
    end
end

