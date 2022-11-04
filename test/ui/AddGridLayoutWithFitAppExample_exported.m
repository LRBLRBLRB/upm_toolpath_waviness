classdef AddGridLayoutWithFitAppExample_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        UITable                      matlab.ui.control.Table
        SelectFontSizeDropDown       matlab.ui.control.DropDown
        SelectFontSizeDropDownLabel  matlab.ui.control.Label
        SelectFontNameDropDown       matlab.ui.control.DropDown
        SelectFontNameDropDownLabel  matlab.ui.control.Label
        UIAxes                       matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.UITable.Data = magic(5);
            plot(app.UIAxes,app.UITable.Data)
        end

        % Value changed function: SelectFontNameDropDown
        function SelectFontNameDropDownValueChanged(app, event)
            fontName = app.SelectFontNameDropDown.Value;
            
            % Get all children that have a 'FontName' property
            children = allchild(app.UIFigure);
            childrenWithFontNameProperty = findall(children,'flat','-property','FontName');
            
            if fontName == "Helvetica"
                set(childrenWithFontNameProperty,'FontName','Helvetica')
            else
                set(childrenWithFontNameProperty,'FontName','Lucida Console')
            end
        end

        % Value changed function: SelectFontSizeDropDown
        function SelectFontSizeDropDownValueChanged(app, event)
            fontSize = app.SelectFontSizeDropDown.Value;
            
            % Get all children that have a 'FontSize' property
            children = allchild(app.UIFigure);
            childrenWithFontSizeProperty = findall(children,'flat','-property','FontSize');
           
            if fontSize == "12"
                set(childrenWithFontSizeProperty,'FontSize',12)
            else
                set(childrenWithFontSizeProperty,'FontSize',24)
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 742 484];
            app.UIFigure.Name = 'UI Figure';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.XTickLabelRotation = 0;
            app.UIAxes.YTickLabelRotation = 0;
            app.UIAxes.ZTickLabelRotation = 0;
            app.UIAxes.Position = [16 53 354 396];

            % Create SelectFontNameDropDownLabel
            app.SelectFontNameDropDownLabel = uilabel(app.UIFigure);
            app.SelectFontNameDropDownLabel.HorizontalAlignment = 'right';
            app.SelectFontNameDropDownLabel.Position = [451 430 102 22];
            app.SelectFontNameDropDownLabel.Text = 'Select Font Name';

            % Create SelectFontNameDropDown
            app.SelectFontNameDropDown = uidropdown(app.UIFigure);
            app.SelectFontNameDropDown.Items = {'Helvetica', 'Lucida Console', ''};
            app.SelectFontNameDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectFontNameDropDownValueChanged, true);
            app.SelectFontNameDropDown.Position = [568 430 100 22];
            app.SelectFontNameDropDown.Value = 'Helvetica';

            % Create SelectFontSizeDropDownLabel
            app.SelectFontSizeDropDownLabel = uilabel(app.UIFigure);
            app.SelectFontSizeDropDownLabel.HorizontalAlignment = 'right';
            app.SelectFontSizeDropDownLabel.Position = [460 379 93 22];
            app.SelectFontSizeDropDownLabel.Text = 'Select Font Size';

            % Create SelectFontSizeDropDown
            app.SelectFontSizeDropDown = uidropdown(app.UIFigure);
            app.SelectFontSizeDropDown.Items = {'12', '24'};
            app.SelectFontSizeDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectFontSizeDropDownValueChanged, true);
            app.SelectFontSizeDropDown.Position = [568 379 100 22];
            app.SelectFontSizeDropDown.Value = '12';

            % Create UITable
            app.UITable = uitable(app.UIFigure);
            app.UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'; 'Column 5'};
            app.UITable.RowName = {};
            app.UITable.Position = [423 53 274 280];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AddGridLayoutWithFitAppExample_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

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