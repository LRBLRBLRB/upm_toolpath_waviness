function isContinue = checksavedfig
fig = uifigure('Name','Select the saved figure', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);

fig.Position(3:4) = [300,150];
figGridLayout = uigridlayout(fig,[2,2]);
figGridLayout.RowHeight = {'1x','fit'};
figGridLayout.ColumnWidth = {'1x','1x'};

figLabel = uilabel(figGridLayout,'Text', ...
    {['Any of the following operation can be used ', ...
    '  to select the saved figure:'], ...
    '',...
    '  - Shift-click the left mouse button.',...
    '  - Click the middle mouse button.',...
    '  - Click both left and right mouse buttons.'}, ...
    'HorizontalAlignment','left','VerticalAlignment','center', ...
    'FontSize',14,'WordWrap','on');
figLabel.Layout.Row = 1;
figLabel.Layout.Column = [1,2];

figEnterButton = uibutton(figGridLayout,'push','Text','Enter');
figEnterButton.Layout.Row = 2;
figEnterButton.Layout.Column = 1;
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figCancelButton = uibutton(figGridLayout,'push','Text','Cancel');
figCancelButton.Layout.Row = 2;
figCancelButton.Layout.Column = 2;
figCancelButton.ButtonPushedFcn = @(app,event)fig_cancel_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function fig_close_req(app)
        isContinue = 0;
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        delete(fig);
    end

    function fig_cancel_button_pushed
        isContinue = 0;
        delete(figShow);
        delete(fig);
    end

end