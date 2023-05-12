function [isContinue,approxOut] = selectfr(textFontType,textFontSize)
%SELECTFR select feedrate smoothing (approximation) method
%   

fig = uifigure('Name','Feed rate approximation selection', ...
                'WindowStyle','alwaysontop','Visible','off');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [400,300];
figGridLayout = uigridlayout(fig,[3,2]);
figGridLayout.RowHeight = {'1x','fit'};
figGridLayout.ColumnWidth = {'1x','1x'};
color = fig.Color;

figParamPanel = uipanel('Title','Parameters Selection','TitlePosition','centertop', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    "Scrollable","on",'BorderType','etchedin');
figParamPanel.Layout.Row = 1;
figParamPanel.Layout.Column = [1,2];

figParamGridLayout = uigridlayout(figParamPanel,[2,2]);
figParamGridLayout.RowHeight = {'1x','1x'};
figParamGridLayout.ColumnWidth = {'fit','1x'};

frMethodLabel = uilabel(figParamGridLayout,'Text','Fr Method');
frMethodLabel.RowHeight = 1;
frMethodLabel.ColumnWidth = 1;
frMethodDropdown = uidropdown(figParamGridLayout,'Items',{'Cubic Spline','Basic Spline'}, ...
    'Value','Cubic Spline','Editable','off');
frMethodDropdown.RowHeight = 1;
frMethodDropdown.ColumnWidth = 2;

frParamLabel = uilabel(figParamGridLayout,'Text','Fr Parameter');
frParamLabel.RowHeight = 2;
frParamLabel.ColumnWidth = 1;
frParamEditfield = uieditfield(figParamGridLayout,'numeric','Value',1);
frParamEditfield.RowHeight = 2;
frParamEditfield.ColumnWidth = 2;

figEnterButton = uibutton(figGridLayout,'push','Text','OK & Continue','Visible','on');
figEnterButton.Layout.Row = 2;
figEnterButton.Layout.Column = 1;
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figCancelButton = uibutton(figGridLayout,'push','Text','Cancel','Visible','on');
figCancelButton.Layout.Row = 2;
figCancelButton.Layout.Column = 2;
figCancelButton.ButtonPushedFcn = @(app,event)fig_cancel_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function valuechanged_fcn(app)
        approxOut.approxMethod = frMethodDropdown.Value;
        approxOut.approxParam = frParamEditfield.Value;
    end

    function fig_close_req(app)
        isContinue = 0;
        valuechanged_fcn(app);
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        valuechanged_fcn(app)
        fig_close_req(fig);
    end

    function fig_cancel_button_pushed
        isContinue = 0;
        valuechanged_fcn(app)
        fig_close_req(fig);
    end


end

