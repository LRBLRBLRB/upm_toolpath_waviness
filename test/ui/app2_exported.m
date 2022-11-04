classdef app2_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure    matlab.ui.Figure
        Toolbar     matlab.ui.container.Toolbar
        PushTool    matlab.ui.container.toolbar.PushTool
        ToggleTool  matlab.ui.container.toolbar.ToggleTool
        TabGroup    matlab.ui.container.TabGroup
        Tab         matlab.ui.container.Tab
        Tab2        matlab.ui.container.Tab
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create Toolbar
            app.Toolbar = uitoolbar(app.UIFigure);

            % Create PushTool
            app.PushTool = uipushtool(app.Toolbar);
            app.PushTool.Icon = fullfile(pathToMLAPP, '1.png');

            % Create ToggleTool
            app.ToggleTool = uitoggletool(app.Toolbar);
            app.ToggleTool.Icon = fullfile(pathToMLAPP, '2.png');

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [1 0 640 481];

            % Create Tab
            app.Tab = uitab(app.TabGroup);
            app.Tab.Tooltip = {'ded'};
            app.Tab.Title = 'Tab';

            % Create Tab2
            app.Tab2 = uitab(app.TabGroup);
            app.Tab2.Title = 'Tab2';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app2_exported

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