classdef post_process < matlab.apps.AppBase

    properties (Constant)
        workpieceCoordDefault
        toolCoordDefault
    end

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        UIGridLayout    matlab.ui.container.GridLayout
    end

    properties (Access = public)
        workpieceCoord  string
        toolCoord       string
    end

    methods (Access = private)
        function UIFigureCloseReq(app)
            selection = uiconfirm(app.UIFigure, 'Close the figure window?', ...
                'Confirmation');
            switch selection
                case 'OK'
                    delete(app.UIFigure);
                case 'Cancel'
                    return
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible','off','WindowStyle','alwaysontop');
            app.UIFigure.Position(3:4) = [300, 150];
            app.UIFigure.Name = 'Select the post process parameters';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app,@UIFigureCloseReq,true);

            app.UIGridLayout = uigridlayout(app.UIFigure, [2 2]);
            app.UIGridLayout.ColumnWidth = {'1x', '1x'};
            app.UIGridLayout.RowHeight = {'1x', '1x'};

            % 

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end

        % Code that executes after component creation
        function startupFcn(app)
            % default set
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = post_process
            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app,@startupFcn);
            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure);
        end
    end
end
