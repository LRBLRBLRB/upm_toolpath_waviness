function isContinue = checkextractcir(hPlot)

fig = uifigure('Name','Check the useless extraction', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);

fig.Position(3:4) = [300,150];
figGridLayout = uigridlayout(fig,[2,2]);
figGridLayout.RowHeight = {'1x','fit'};
figGridLayout.ColumnWidth = {'1x','1x'};

figLabel = uilabel(figGridLayout,'Text','Useless extraction is finished.', ...
    'HorizontalAlignment','center','FontSize',14);
figLabel.Layout.Row = 1;
figLabel.Layout.Column = [1,2];

figEnterButton = uibutton(figGridLayout,'push','Text','Enter','Visible','on');
figEnterButton.Layout.Row = 3;
figEnterButton.Layout.Column = 1;
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figReFitButton = uibutton(figGridLayout,'push','Text','More','Visible','on');
figReFitButton.Layout.Row = 3;
figReFitButton.Layout.Column = 2;
figReFitButton.ButtonPushedFcn = @(app,event)fig_refit_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        fig_close_req(fig);
    end

    function fig_refit_button_pushed
        isContinue = 0;
        delete(hPlot);
        fig_close_req(fig);
    end
end

