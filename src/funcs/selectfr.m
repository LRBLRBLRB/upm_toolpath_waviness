function [approxOut,isContinue] = selectfr(dataX,dataY,textFontType,textFontSize)
%SELECTFR select feedrate smoothing (approximation) method
%   

dlgOpts.Interpreter = 'tex';
dlgOpts.WindowStyle = 'modal';

fig = uifigure('Name','Feed rate approximation selection', ...
                'Visible','off','Scrollable','on');
fig.CloseRequestFcn = @(app,event)fig_close_req(app);
fig.Position(3:4) = [400,250];
figGridLayout = uigridlayout(fig,[4,4]);
figGridLayout.RowHeight = {'fit','1x','fit','fit'};
figGridLayout.ColumnWidth = {'fit','1x','1x','fit'};
% color = fig.Color;

% Eq = texlabel('');
% 
% frEquationTextarea = uitextarea(figGridLayout,'Value',Eq,'Editable','off', ...
%     'FontName','Symath','FontSize',textFontSize);
% frEquationTextarea.Layout.Row = 1;
% frEquationTextarea.Layout.Column = [1,4];

%% Categories dropdown
methodCateDropdown = uidropdown(figGridLayout,'Editable','off', ...
    'Items',{'','CurveFitter','polyfit','csape','csaps','spaps','Basic Spline'}, ...
    'Value','','Placeholder','Fitters', ...
    'FontName',textFontType,'FontSize',textFontSize);
methodCateDropdown.Layout.Row = 1;
methodCateDropdown.Layout.Column = 1;
methodCateDropdown.ValueChangedFcn = @(app,event)fitter_value_changed(app,event);

%% Basic spline Parameters
paramPanel = uipanel(figGridLayout,'Title','Parameters', ...
    'TitlePosition','centertop','BorderType','line', ...
    'FontName',textFontType,'FontSize',textFontSize);
paramPanel.Layout.Row = [1,3];
paramPanel.Layout.Column = [2,4];
paramGridLayout = uigridlayout(paramPanel,[2,2]);
paramGridLayout.RowHeight = {'fit','fit'};
paramGridLayout.ColumnWidth = {'fit','1x'};

blankLabel = uilabel(paramGridLayout,'Text', ...
    {'No parameters to set.','Push the Launch button.'}, ...
    'FontName',textFontType,'FontSize',textFontSize);
blankLabel.Layout.Row = 1;
blankLabel.Layout.Column = 1;

appLabel = uilabel(paramGridLayout,'Text','curveFitter file', ...
    'FontName',textFontType,'FontSize',textFontSize);
appLabel.Visible = 'off';
appLabel.Layout.Row = 1;
appLabel.Layout.Column = 1;
appButton = uibutton(paramGridLayout,'push','Text','Choose', ...
    'FontName',textFontType,'FontSize',textFontSize);
appButton.Visible = 'off';
appButton.Layout.Row = 1;
appButton.Layout.Column = 2;
appButton.ButtonPushedFcn = @(app,event)appButton_button_pushed;
appEditfield = uieditfield(paramGridLayout,'text');
appEditfield.Visible = 'off';
appEditfield.Layout.Row = 2;
appEditfield.Layout.Column = [1,2];
appEditfield.ValueChangedFcn = @(app,event)appEditfield_value_changed(event);

nurbsLabel = uilabel(paramGridLayout,'Text','Method');
nurbsLabel.Layout.Row = 1;
nurbsLabel.Layout.Column = 1;
nurbsLabel.Visible = 'off';
nurbsDropdown = uidropdown(paramGridLayout,'Items',{'Basic Spline','NURBS'}, ...
    'Value','Basic Spline','Editable','off');
nurbsDropdown.Layout.Row = 1;
nurbsDropdown.Layout.Column = 2;
nurbsDropdown.Visible = 'off';

nurbsParamLabel = uilabel(paramGridLayout,'Text','Parameter', ...
    'FontName',textFontType,'FontSize',textFontSize);
nurbsParamLabel.Layout.Row = 2;
nurbsParamLabel.Layout.Column = 1;
nurbsParamLabel.Visible = 'off';
nurbsParamEditfield = uieditfield(paramGridLayout,'numeric','Value',-1);
nurbsParamEditfield.Layout.Row = 2;
nurbsParamEditfield.Layout.Column = 2;
nurbsParamEditfield.Visible = 'off';

LaunchButton = uibutton(figGridLayout,'push','Text','Launch the fitter', ...
    'FontName',textFontType,'FontSize',textFontSize, ...
    'ButtonPushedFcn',@(app,event)launch_button_pushed);
LaunchButton.Layout.Row = 3;
LaunchButton.Layout.Column = 1;

%% buttons
figEnterButton = uibutton(figGridLayout,'push','Text','OK & Continue','Visible','on');
figEnterButton.Layout.Row = 4;
figEnterButton.Layout.Column = [1,2];
figEnterButton.ButtonPushedFcn = @(app,event)fig_enter_button_pushed;

figCancelButton = uibutton(figGridLayout,'push','Text','Cancel','Visible','on');
figCancelButton.Layout.Row = 4;
figCancelButton.Layout.Column = [3,4];
figCancelButton.ButtonPushedFcn = @(app,event)fig_cancel_button_pushed;

approxOut = struct('fittedmodel',[],'param',0,'goodness',[],'output',[]);
isContinue = 0;
curveFitterFile = [];

fig.Visible = 'on';
methodCateDropdown.Visible = 'on';

uiwait(fig);

    function fitter_value_changed(src,evt)
        switch evt.Value
            case ''
                set_all_invisible;
                blankLabel.Visible = 'on';
                warndlg({sprintf(['\\fontname{%s}\\fontSize{%d} ', ...
                    'No fitter has been selected.\n'], ...
                    textFontType,textFontSize), ...
                    '\bf Curve fitter will not begin'},'Warning',dlgOpts);
            case 'CurveFitter'
                set_all_invisible;
                appLabel.Visible = 'on';
                appButton.Visible = 'on';
                appEditfield.Visible = 'on';
            case 'polyfit'
                set_all_invisible;
            case 'csape'
                set_all_invisible;
            case 'csaps'
                set_all_invisible;
            case 'spaps'
                set_all_invisible;
            case 'Basic Spline'
                set_all_invisible;
                nurbsLabel.Visible = 'on';
                nurbsDropdown.Visible = 'on';
                nurbsParamLabel.Visible = 'on';
                nurbsParamEditfield.Visible = 'on';
        end
%         fittype
        disp(['Previous: ' evt.PreviousValue]);
        disp(['Current: ' evt.Value]);
        disp('------------------');
    end

    % Code that workspace directory editfield changed
    function appEditfield_value_changed(event)
        curveFitterFile = event.Value;
        mkdir(curveFitterFile);
    end
    
    % Code that choose the workspace directory
    function appButton_button_pushed
        curveFitterFile = uigetdir('Select the Workspace Directory');
        if isempty(curveFitterFile)
            uialert(fig,{'Invalid workspace directory:', ...
                'workspace directory should not be empty'},'Alert Message');
            return;
        end
        appEditfield.Value = curveFitterFile;
    end

    function launch_button_pushed
        switch methodCateDropdown.Value
            case '' 
                warndlg({sprintf(['\\fontname{%s}\\fontsize{%d} ' ...
                    'No fitter has been selected.\n'], ...
                    textFontType,textFontSize), ...
                    '\bf Curve fitter cannot starts.'},'Warning',dlgOpts);
            case 'CurveFitter'
                curveFitter(dataX,dataY);
                approxOut.approxParam = nurbsParamEditfield.Value;
            case 'polyfit'
                fitOrder = 5;
                fitType = sprintf('poly%d',fitOrder);
                fit(dataX,dataY,'poly1')
            case 'csape'
            case 'csaps'
            case 'spaps'
            case 'Basic Spline'
                approxOut.approxResult = spap2(3,3,dataX,dataY,approxOut.approxParam);
        end
    end

    function set_all_invisible
        blankLabel.Visible = 'off';
        appLabel.Visible = 'off';
        appButton.Visible = 'off';
        appEditfield.Visible = 'off';
        nurbsLabel.Visible = 'off';
        nurbsDropdown.Visible = 'off';
        nurbsParamLabel.Visible = 'off';
        nurbsParamEditfield.Visible = 'off';
    end

    function fig_close_req(app)
        isContinue = 0;
        delete(app);
    end

    function fig_enter_button_pushed
        isContinue = 1;
        fig_close_req(fig);
    end

    function fig_cancel_button_pushed
        isContinue = 0;
        fig_close_req(fig);
    end


end

