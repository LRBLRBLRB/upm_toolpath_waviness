classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure         matlab.ui.Figure
        GridLayout       matlab.ui.container.GridLayout
        TextArea         matlab.ui.control.TextArea
        TextAreaLabel    matlab.ui.control.Label
        Label            matlab.ui.control.Label
        Panel            matlab.ui.container.Panel
        GridLayout2      matlab.ui.container.GridLayout
        Button           matlab.ui.control.Button
        EditField        matlab.ui.control.EditField
        EditFieldLabel   matlab.ui.control.Label
        GridLayout3      matlab.ui.container.GridLayout
        Button2          matlab.ui.control.Button
        EditField2       matlab.ui.control.EditField
        EditField2Label  matlab.ui.control.Label
        UIAxes           matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            delete(app)
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 638 482];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.WindowStyle = 'alwaysontop';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'3x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Layout.Row = 4;
            app.UIAxes.Layout.Column = [2 4];

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout);
            app.GridLayout3.ColumnWidth = {'fit', '1x', 'fit'};
            app.GridLayout3.RowHeight = {'fit'};
            app.GridLayout3.Layout.Row = 3;
            app.GridLayout3.Layout.Column = [1 4];

            % Create EditField2Label
            app.EditField2Label = uilabel(app.GridLayout3);
            app.EditField2Label.HorizontalAlignment = 'right';
            app.EditField2Label.Layout.Row = 1;
            app.EditField2Label.Layout.Column = 1;
            app.EditField2Label.Text = 'Edit Field2';

            % Create EditField2
            app.EditField2 = uieditfield(app.GridLayout3, 'text');
            app.EditField2.Layout.Row = 1;
            app.EditField2.Layout.Column = 2;

            % Create Button2
            app.Button2 = uibutton(app.GridLayout3, 'push');
            app.Button2.Layout.Row = 1;
            app.Button2.Layout.Column = 3;
            app.Button2.Text = 'Button2';

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.GridLayout);
            app.GridLayout2.ColumnWidth = {'fit', '1x', 'fit'};
            app.GridLayout2.RowHeight = {'fit'};
            app.GridLayout2.Layout.Row = 2;
            app.GridLayout2.Layout.Column = [1 4];

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.GridLayout2);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Layout.Row = 1;
            app.EditFieldLabel.Layout.Column = 1;
            app.EditFieldLabel.Text = 'Edit Field';

            % Create EditField
            app.EditField = uieditfield(app.GridLayout2, 'text');
            app.EditField.Layout.Row = 1;
            app.EditField.Layout.Column = 2;

            % Create Button
            app.Button = uibutton(app.GridLayout2, 'push');
            app.Button.Layout.Row = 1;
            app.Button.Layout.Column = 3;

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Title = 'Panel';
            app.Panel.Layout.Row = [4 5];
            app.Panel.Layout.Column = 1;

            % Create Label
            app.Label = uilabel(app.GridLayout);
            app.Label.Layout.Row = 1;
            app.Label.Layout.Column = [1 4];

            % Create TextAreaLabel
            app.TextAreaLabel = uilabel(app.GridLayout);
            app.TextAreaLabel.HorizontalAlignment = 'right';
            app.TextAreaLabel.Layout.Row = 6;
            app.TextAreaLabel.Layout.Column = 1;
            app.TextAreaLabel.Text = 'Text Area';

            % Create TextArea
            app.TextArea = uitextarea(app.GridLayout);
            app.TextArea.Layout.Row = 6;
            app.TextArea.Layout.Column = [2 4];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end