function [approxOut,isContinue] = selectfr(textFontType,textFontSize)
%SELECTFR select feedrate smoothing (approximation) method
%   

fig = uifigure('Name','Feed rate approximation selection', ...
                'Visible','off','Scrollable','on', ...
                'WindowStyle','alwaysontop');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [300,200];
figGridLayout = uigridlayout(fig,[3,4]);
figGridLayout.RowHeight = {'1x','1x','fit'};
figGridLayout.ColumnWidth = {'1x','1x','1x','1x'};
% color = fig.Color;

% Eq = texlabel('');
% 
% frEquationTextarea = uitextarea(figGridLayout,'Value',Eq,'Editable','off', ...
%     'FontName','Symath','FontSize',textFontSize);
% frEquationTextarea.Layout.Row = 1;
% frEquationTextarea.Layout.Column = [1,4];

frMethodLabel = uilabel(figGridLayout,'Text',{'Approximation', 'Method'});
frMethodLabel.Layout.Row = 1;
frMethodLabel.Layout.Column = 1;
frMethodDropdown = uidropdown(figGridLayout,'Items',{'Cubic Spline','Basic Spline'}, ...
    'Value','Cubic Spline','Editable','off');
frMethodDropdown.Layout.Row = 1;
frMethodDropdown.Layout.Column = [2,3];

frParamLabel = uilabel(figGridLayout,'Text',{'Approximation', 'Parameter'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
frParamLabel.Layout.Row = 2;
frParamLabel.Layout.Column = 1;
frParamEditfield = uieditfield(figGridLayout,'numeric','Value',-1);
frParamEditfield.Layout.Row = 2;
frParamEditfield.Layout.Column = [2,3];

figEnterButton = uibutton(figGridLayout,'push','Text','OK & Continue','Visible','on');
figEnterButton.Layout.Row = 3;
figEnterButton.Layout.Column = [1,2];
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figCancelButton = uibutton(figGridLayout,'push','Text','Cancel','Visible','on');
figCancelButton.Layout.Row = 3;
figCancelButton.Layout.Column = [3,4];
figCancelButton.ButtonPushedFcn = @(app,event)fig_cancel_button_pushed;

fig.Visible = 'on';

uiwait(fig);

    function valuechanged_fcn
        switch frMethodDropdown.Value
            case 'Cubic Spline'
                approxOut.approxMethod = 'csaps';
            case 'Basic Spline'
                approxOut.approxMethod = 'spap2';
        end
        approxOut.approxParam = frParamEditfield.Value;
    end

    function fig_close_req(app)
        isContinue = 0;
        valuechanged_fcn;
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        valuechanged_fcn;
        fig_close_req(fig);
    end

    function fig_cancel_button_pushed
        isContinue = 0;
        valuechanged_fcn;
        fig_close_req(fig);
    end


end

