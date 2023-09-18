classdef upm_toolpath_gui < matlab.apps.AppBase
    % to create a ui figure for tool data & parameters input

    % Default properties
    properties (Constant)
        workspaceDirDefault     = 'D:\Code\2021-11_ToolWaviness\upm_toolpath_waviness\workspace';
        unitDefault             = '\mum'
        toolUnitDefault         = 'mm'
        fontNameDefault         = 'Times New Roman'
        fontSizeDefault         = 12
        toolFitTypeDefault      = 'lineArc'
        arcFitMethodDefault     = 'levenberg-marquardt'
        arcRansacMaxDistDefault = 1e-2
        lineFitMethodDefault    = 'polyfit'
        lineFitMaxDistDefault   = 1
        paramMethodDefault      = 'centripetal'
        radius0Default          = 192
        ImportCell              = {'Point Cloud','Mesh'}
        Geometry2DCell          = {'Rotating Paraboloid','Aspheric'}
        Geometry3DCell          = {'Ellipsoid','Function-Based'}
        cutDirectionDefault     = 'Edge to Center'
        startDirectionDefault   = 'X Minus'
        angularDiscreteDefault  = 'Constant Arc'
        radialIncrementDefault  = 'On-Axis'
        spiralMethodDefault     = 'Radius-Number'
        zAllowanceDefault        = 1.2
        aimResDefault           = 0.0005
        rStepDefault            = 0.2
        maxIterDefault          = 50
        arcLengthDefault        = 0.03
        maxAngPtDistDefault     = 6*pi/180
        angularLengthDefault    = 6*pi/180
        surfPlotSparDefault     = 501
        
    end

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        FigureMenuFiles         matlab.ui.container.Menu
        FigureMenuFilesNew      matlab.ui.container.Menu
        FigureMenuFilesPrint     matlab.ui.container.Menu
        FigureMenuHelp          matlab.ui.container.Menu
        FigureMenuHelpDoc       matlab.ui.container.Menu
        FigureMenu              matlab.ui.container.Menu
        FigureToolbar           matlab.ui.container.Toolbar
        AxisEqualPushtool       matlab.ui.container.toolbar.PushTool
        CloseAllFigurePushtool  matlab.ui.container.toolbar.PushTool
        BoldInfoToggletool      matlab.ui.container.toolbar.ToggleTool
        TopToggletool           matlab.ui.container.toolbar.ToggleTool
        FigureGl                matlab.ui.container.GridLayout
        FigureTbGp              matlab.ui.container.TabGroup
        ToolTb                  matlab.ui.container.Tab
        WorkspaceDirEf          matlab.ui.control.EditField
        WorkspaceDirBtn         matlab.ui.control.Button
        ParamTabToolImport      matlab.ui.container.Tab
        UnitDd                  matlab.ui.control.DropDown
        FontNameDd              matlab.ui.control.DropDown
        FontSizeEf              matlab.ui.control.NumericEditField
        CommonResetBtn          matlab.ui.control.Button
        ToolFileEf              matlab.ui.control.EditField
        ToolFileBtn             matlab.ui.control.Button
        ToolUnitDd              matlab.ui.control.DropDown
        ToolimportResetBtn      matlab.ui.control.Button
        ToolimportUpdateBtn     matlab.ui.control.Button
        ToolimportPlotBtn       matlab.ui.control.Button
        ParamTabToolfit         matlab.ui.container.Tab
        ToolFitTypeDd           matlab.ui.control.DropDown
        LineFitMethodDd         matlab.ui.control.DropDown
        LineFitMaxDistEf        matlab.ui.control.NumericEditField
        ArcFitMethodDd          matlab.ui.control.DropDown
        ArcRansacMaxDistEf      matlab.ui.control.NumericEditField
        Radius0Ef               matlab.ui.control.NumericEditField
        S1ToolExtract2DLineBtn  matlab.ui.control.Button
        S1ToolExtract3DLineBtn  matlab.ui.control.Button
        S1ToolExtractSurfBtn    matlab.ui.control.Button
        ToolfitResetBtn         matlab.ui.control.Button
        ToolfitUpdateBtn        matlab.ui.control.Button
        ToolfitPlotBtn          matlab.ui.control.Button
        ParamTabToolmod         matlab.ui.container.Tab
        S1Tool2DBtn             matlab.ui.control.Button
        S1Tool3DBtn             matlab.ui.control.Button
        ParamTabToolinterp      matlab.ui.container.Tab
        S1ToolModelBtn          matlab.ui.control.Button
        ParamMethodDd           matlab.ui.control.DropDown
        ToolinterpResetBtn      matlab.ui.control.Button
        ToolinterpUpdateBtn     matlab.ui.control.Button
        ToolinterpPlotBtn       matlab.ui.control.Button
        ToolDataAxes            matlab.ui.control.UIAxes
        ToolDataAxesClearBtn    matlab.ui.control.Button
        SurfaceTb               matlab.ui.container.Tab
        AddSurfaceBtn           matlab.ui.control.Button
        SurfacePlotSparSpin     matlab.ui.control.Spinner
        SurfaceUnitDd           matlab.ui.control.DropDown
        SurfaceDetailTa         matlab.ui.control.TextArea
        SurfaceDataAxes         matlab.ui.control.UIAxes
        Surface2DPlotBtn        matlab.ui.control.Button
        Surface3DPlotBtn        matlab.ui.control.Button
        SurfaceSavedBtn         matlab.ui.control.Button
        SurfaceCancelBtn        matlab.ui.control.Button
        ProgramTb               matlab.ui.container.Tab
        OptimTb                 matlab.ui.container.Tab
        CheckToolLamp           matlab.ui.control.Lamp
        CheckSurfLamp           matlab.ui.control.Lamp
        OptimParamPathTab       matlab.ui.container.Tab
        CutDirectionDd          matlab.ui.control.DropDown
        StartDirectionDd        matlab.ui.control.DropDown
        RadialIncrementDd       matlab.ui.control.DropDown
        AngularDiscreteDd       matlab.ui.control.DropDown
        AimResEf                matlab.ui.control.NumericEditField
        MaxIterSpin             matlab.ui.control.Spinner
        RStepEf                 matlab.ui.control.NumericEditField
        ArcLengthLb             matlab.ui.control.Label
        ArcLengthEf             matlab.ui.control.NumericEditField
        ArcLengthUnitLb         matlab.ui.control.Label
        MaxAngPtDistLb          matlab.ui.control.Label
        MaxAngPtDistEf          matlab.ui.control.NumericEditField
        MaxAngPtDistUnitLb      matlab.ui.control.Label
        AngularLengthLb         matlab.ui.control.Label
        AngularLengthEf         matlab.ui.control.NumericEditField
        AngularLengthUnitLb     matlab.ui.control.Label
        ZAllowanceEf            matlab.ui.control.NumericEditField
        SpiralMethodDd          matlab.ui.control.DropDown
        OptimParamFeedTab       matlab.ui.container.Tab
        OptimUpdateBtn          matlab.ui.control.Button
        OptimResetBtn           matlab.ui.control.Button
        Optim2DSingleBtn        matlab.ui.control.Button
        Optim2DSingleIterBtn    matlab.ui.control.Button
        Optim2DMultiBtn         matlab.ui.control.Button
        Optim2DMultiIterBtn     matlab.ui.control.Button
        InfoTa                  matlab.ui.control.TextArea
        MsgState                logical
        Msg                     char
        MsgNum                  int16
        S2DesignSimulAsphericConcentricBtn  matlab.ui.control.Button
        S2DesignSimulFreeformBtn            matlab.ui.control.Button
    end

    % properties for the opening app
    properties (Access = private)
        AddSurfaceApp
    end
    
    % properties that should be used in the .m program
    properties (Access = public)
        saveFig
        workspaceDir                string
        toolPathName                string
        unit                        char
        fontName                    string
        fontSize            (1,1)   double
        toolUnit                    char
        toolFitType                 string
        arcFitMethod                string
        arcRansacMaxDist    (1,1)   double
        lineFitMethod       (1,1)   string
        lineFitMaxDist      (1,1)   double
        paramMethod                 string
        radius0             (1,1)   double
        toolOri
        % post process between tool fit and interpolation
        toolFit
        radius              (1,1)   double
        openAngle
        toolDataFile                char
        toolData                    struct

        surfaceUnit                 char
        surfType                    char
        surfDomain          (:,2)   double
        surfPathName                char
        surfPt              (3,:)   double
        surfFuncs                   function_handle
        surfFx                      sym
        surfFy                      sym
        surfPlotSpar        (1,1)   double
        surfMesh            (:,:,3) double

        cutDirection                string
        startDirection              string
        angularDiscrete             string
        arcLength           (1,1)   double
        maxAngPtDist        (1,1)   double
        angularLength       (1,1)   double
        radialIncrement             string
        aimRes              (1,1)   double
        rStep               (1,1)   double
        maxIter             (1,1)   double
        rMax                (1,1)   double
        spiralMethod                string
        zAllowance           (1,1)   double
    end

    % Functions that update the properties when multi-apps are used
    methods (Access = public)
        function updateSurface(app,surfType,surfString,surfDomain)
            % Stores the inputs as properties
            app.surfType = surfType;
            app.surfDomain = surfDomain;
            % app.SurfaceDetailEf.Value = surfDetail;

            % Update surface detail
            app.SurfaceDetailTa.Value = {['- Surface type: ',app.surfType]};
            switch app.surfType
                case 'Point Cloud'
                    app.SurfaceDetailTa.Value{2} = '- Surface path: ';
                    app.SurfaceDetailTa.Value{3} = char(surfString);
                    app.SurfaceDetailTa.Value{4} = '- Surface Domain: ';
                case 'Function-Based'
                    app.SurfaceDetailTa.Value{2} = '- Surface function: ';
                    app.SurfaceDetailTa.Value{3} = char(surfString);
                    app.SurfaceDetailTa.Value{4} = '- Surface Domain: ';
                otherwise
                    app.SurfaceDetailTa.Value{2} = '- Surface Function: ';
                    % app.SurfaceDetailTa.Value{3} = ...
                    %     '$Z = \sqrt{C^{2}(D-\frac{x^{2}}{A^{2}}-\frac{y^{2}}{B^{2}})}$';
                    app.SurfaceDetailTa.Value{3} = char(surfString);
                    app.SurfaceDetailTa.Value{4} = '- Surface Domain: ';
                    app.surfFuncs = matlabFunction(surfString);
            end

            switch size(app.surfDomain,1)
                case 2
                    app.SurfaceDetailTa.Value{5} = ['    X: ',num2str(app.surfDomain(1,:)),' (',app.unit,')'];
                    app.SurfaceDetailTa.Value{6} = ['    Y: ',num2str(app.surfDomain(2,:)),' (',app.unit,')'];
                case 1
                    app.SurfaceDetailTa.Value{5} = ['    R: ',num2str(app.surfDomain),' (',app.unit,')'];
            end

            % Re-enable the AddSurface buttion
            app.AddSurfaceBtn.Enable = 'on';
        end
    end
    
    % Functions that is used in the callbacks
    methods (Access = private)
        resetCommonParams(app);

        resetToolfitParams(app);

        resetSurfaceParams(app);

        resetOptimParams(app);

        toolFitCheck(app);

        surfaceMeshGen(app);

        surfacePlot(app,method);

        function [figName,figFormat] = saveAxes2Fig(app)
            [figFileName,figDirName,figIndex] = uiputfile({ ...
                '*.fig','MATLAB Figure(*.fig)'; ...
                '*.png','Portable Network Graphics file(*.png)';...
                '*.jpg','JPEG image(*.jpg)';...
                '*.bmp','Bitmap file(*.bmp)';...
                '*.tif','TIFF image(*.tif)';...
                '*.tif','TIFF no compression image(*.tif)';...
                '*.svg','Scalable Vector Graphics fig(*.svg)';...
                '*.pdf','Full page Portable Document Format(*.pdf)'...
                }, ...
                'Select a file to save the figure');
            switch figIndex
                case 1
                    figFormat = 'fig';
                case 2
                    figFormat = 'png';
                case 3
                    figFormat = 'jpeg';
                case 4
                    figFormat = 'bmp';
                case 5
                    figFormat = 'tif';
                case 6
                    figFormat = 'tiffn';
                case 7
                    figFormat = 'svg';
                case 8
                    figFormat = 'pdf';
            end
            figName = fullfile(figDirName,figFileName);
        end
    end

    % Callbacks that handle component events
    methods (Access = private)
        % Code that executes after component creation
        function startupFcn(app)
            app.MsgNum = 1;
            app.Msg = 'APP is successfully initialized!';
            InfoTaValueChanged(app,true);
            resetCommonParams(app);
            resetToolfitParams(app);
            resetSurfaceParams(app);
            resetOptimParams(app);
            app.unit = app.UnitDd.Value;
            % updateSurface(app,app.surfType,app.surfString);
            addpath(genpath('funcs'));
            app.Msg = 'All the parameters have been reset.';
            InfoTaValueChanged(app,true);
            app.Msg = 'Please choose a directory for the workspace.';
            InfoTaValueChanged(app,true);
        end

        % Menu selection function: new file
        function FigureMenuFilesNewMenuSelect(app,event)
            ToolDataAxesClearBtnPushed(app,event);
            SurfaceCancelBtnPushed(app,event);
            startupFcn(app);
            close all;
        end

        % Menu selection function: save figures
        function FigureMenuFilesPrintMenuSelect(app,event)
            [figFileName,figDirName] = uiputfile({ ...
                '*.png','Portable Network Graphics file(*.png)';...
                '*.jpg','JPEG image(*.jpg)';...
                '*.tif','TIFF image(*.tif)';...
                '*.pdf','Full page Portable Document Format(*.pdf)'...
                }, ...
                'Select a file to save the figure');
            exportapp(app.UIFigure,fullfile(figFileName,figDirName));
        end

        function FigureMenuHelpFileMenuSelect(app,event)
            winopen('@upm_toolpath_gui\Documentation\index.html');
        end

        % Clicked function: set the current figure equal
        function AxisEqualPushtoolClicked(app,event)
            % app.UIFigure.CurrentAxes
            axis(app.UIFigure.CurrentAxes,'equal');
        end

        % Clicked toolbar pushtool: close all the figures
        function CloseAllFigurePushtoolClicked(app,event)
            close all;
        end

        % CLicked toolbar toggletool: bold / normalize the infomation font
        function BoldInfoToggletoolClicked(app,event)
            switch app.InfoTa.FontWeight
                case 'normal'
                    app.InfoTa.FontWeight = 'bold';
                case 'bold'
                    app.InfoTa.FontWeight = 'normal';
            end
            app.Msg = ['The font weight of the window is set as ',app.InfoTa.FontWeight,'.'];
            InfoTaValueChanged(app,true);
        end

        % CLicked toolbar toggletool: set the window style
        function TopToggletoolClicked(app,event)
            switch app.UIFigure.WindowStyle
                case 'normal'
                    app.UIFigure.WindowStyle = 'alwaysontop';
                    app.Msg = 'The window is set to be always on top.';
                    InfoTaValueChanged(app,true);
                case 'alwaysontop'
                    app.UIFigure.WindowStyle = 'normal';
                    app.Msg = 'The window is set to be normal.';
                    InfoTaValueChanged(app,true);
            end
        end

        % Button down: shift to the toolbar tab, and refresh the message
        function ToolTbButtonDown(app,event)
            app.Msg = ['Switch to the tool tab. ', ...
                'Please choose a directory for the workspace.'];
            InfoTaValueChanged(app,true);
        end

        % Code that workspace directory editfield changed
        function WorkspaceDirEfValueChanged(app,event)
            app.workspaceDir = app.WorkspaceDirEf.Value;
            [app.MsgState,app.Msg] = mkdir(app.workspaceDir);
            if ~app.MsgState
                InfoTaValueChanged(app,true)
            end
            app.Msg = 'Please choose a tool file.';
            InfoTaValueChanged(app,true);
        end
        
        % Code that choose the workspace directory
        function WorkspaceDirBtnPushed(app,event)
            app.workspaceDir = uigetdir(app.workspaceDirDefault, ...
                'Select the Workspace Directory');
            if isempty(app.workspaceDir)
                uialert(app.UIFigure,{'Invalid workspace directory:', ...
                    'workspace directory should not be empty'},'Alert Message');
                return;
            end
            app.WorkspaceDirEf.Value = app.workspaceDir;
            app.Msg = 'Please choose a tool file.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: DdFontName i.e., font type selection
        function FontNameDdValueChanged(app,event)
            app.fontName = app.FontNameDd.Value;
            set(app.ToolDataAxes,'FontName',app.fontName);

            % Report the infomation
            app.Msg = ['The font name is set as ',app.fontName,'.'];
            InfoTaValueChanged(app,true);
        end
        
        % Value changed function: font size selection
        function FontSizeEfValueChanged(app,event)
            app.fontSize = app.FontSizeEf.Value;
            set(app.ToolDataAxes,'FontSize',app.fontSize);

            % Report the infomation
            app.Msg = ['The font size is set as ',app.fontSize,'.'];
            InfoTaValueChanged(app,true);
        end

        % Value changed function: unit
        function UnitDdValueChanged(app,event)
            app.unit = app.UnitDd.Value;

            % Report the infomation
            app.Msg = ['The unit is set as ',app.unit,'.'];
            InfoTaValueChanged(app,true);
        end

        % Button pushed function: reset the common parameters
        function CommonResetButtonPushed(app,event)
            resetCommonParams(app);

            % Report the infomation
            app.Msg = 'All the parameters in Common tab are reset.';
            InfoTaValueChanged(app,true);
        end

        % Code that tool file path editfield changed
        function ToolFileEfValueChanged(app,event)
            if isempty(app.workspaceDir)
                uialert(app.UIFigure, ...
                    'Workspace directory doesn''t exist, please choose one first!', ...
                    'Alert Message','CloseFcn',createCallbackFcn(app,@WorkspaceDirBtnPushed,true));
                return
            end
            app.toolPathName = app.ToolFileEf.Value;
            app.Msg = 'Please select the parameters and click Update.';
            InfoTaValueChanged(app,true);
        end
        
        % Value changed function: to choose the tool file path
        ToolFileBtnPushed(app,event);
        
        % Button Value changed function: tool import reset
        function ToolImportResetBtnPushed(app,event)
            app.ToolUnitDd.Value = app.toolUnitDefault;
            app.toolUnit = [];
            app.ToolFileEf.Value = '';
            app.toolPathName = app.ToolFileEf.Value;
            app.toolData = [];

%             app.ParamTabToolfit.HandleVisibility = 'on';
%             app.ParamTabToolmod.HandleVisibility = 'on';
%             app.ParamTabToolinterp.HandleVisibility = 'on';

            app.Msg = 'Tool importing process has been reset.';
            InfoTaValueChanged(app,true);
        end

        % Button value changed function: tool import update
        function ToolImportUpdateBtnPushed(app,event)
            if ~strcmp(app.toolUnit,app.ToolUnitDd.Value)
                app.toolUnit = app.ToolUnitDd.Value;
    
                % unit conversion (now: app.toolUnit, aim: app.unit)
                unitList = {'m','mm','\mum','nm'};
                presUnit = find(strcmp(unitList,app.toolUnit),1);
                aimUnit = find(strcmp(unitList,app.unit),1);
                app.toolOri = 1000^(aimUnit - presUnit)*app.toolOri;
            end
        end

        % Button value changed function: tool import plot
        function ToolImportPlotBtnPushed(app,event)
            ToolDataAxesClearBtnPushed(app,event);
            % plot the importing results
            plot(app.ToolDataAxes,app.toolOri(1,:),app.toolOri(2,:), ...
                '.','MarkerSize',2);
            hold(app.ToolDataAxes,'on');
            grid(app.ToolDataAxes,'on');
            xlabel(app.ToolDataAxes,['x (',app.unit,')']);
            ylabel(app.ToolDataAxes,['y (',app.unit,')']);
            title(app.ToolDataAxes,'Tool Original Data');
            app.Msg = 'Tool data is successfully loaded.';
            InfoTaValueChanged(app,true);
        end

        % Dropdown opening function: tool fit type
        function ToolFitTypeDdOpening(app,event)
            app.toolFitType = app.ToolFitTypeDd.Value;
            % enable the parameter selection based on the input tool type
            switch app.toolFitType
                case 'onlyArc'
                    app.ArcRansacMaxDistEf.Enable = "off";
                    app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMethodDd.Enable = "off";
                    app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMaxDistEf.Enable = "off";
                    app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                case 'arcRansac'
                    app.ArcRansacMaxDistEf.Enable = "on";
                    app.ArcRansacMaxDistEf.BackgroundColor = [1 1 1];
                    app.LineFitMethodDd.Enable = "off";
                    app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMaxDistEf.Enable = "off";
                    app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                case 'lineArc'
                    app.ArcRansacMaxDistEf.Enable = "off";
                    app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMethodDd.Enable = "on";
                    app.LineFitMethodDd.BackgroundColor = [1 1 1];
                    app.LineFitMaxDistEf.Enable = "on";
                    app.LineFitMaxDistEf.BackgroundColor = [1 1 1];
            end
        end

        % Value changed function: tool fit type
        function ToolFitTypeDdValueChanged(app,event)
            app.toolFitType = app.ToolFitTypeDd.Value;
            % enable the parameter selection based on the input tool type
            switch app.toolFitType
                case 'onlyArc'
                    app.ArcRansacMaxDistEf.Enable = "off";
                    app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMethodDd.Enable = "off";
                    app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMaxDistEf.Enable = "off";
                    app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                case 'arcRansac'
                    app.ArcRansacMaxDistEf.Enable = "on";
                    app.ArcRansacMaxDistEf.BackgroundColor = [1 1 1];
                    app.LineFitMethodDd.Enable = "off";
                    app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMaxDistEf.Enable = "off";
                    app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                case 'lineArc'
                    app.ArcRansacMaxDistEf.Enable = "off";
                    app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
                    app.LineFitMethodDd.Enable = "on";
                    app.LineFitMethodDd.BackgroundColor = [1 1 1];
                    app.LineFitMaxDistEf.Enable = "on";
                    app.LineFitMaxDistEf.BackgroundColor = [1 1 1];
            end
        end

        % Value changed function: reset the parameters
        function ToolFitResetBtnPushed(app,event)
            % Reset all the values to the default
            resetToolfitParams(app);

            % Report the infomation
            app.Msg = ['All the parameters in Tool File Processing are reset.', ... 
                'Modify them and click ''Update'' to set.'];
            InfoTaValueChanged(app,true);
        end

        % Value changed function: update the selection selection
        function ToolFitUpdateBtnPushed(app,event)
            app.toolFitType = app.ToolFitTypeDd.Value;
            app.arcFitMethod = app.ArcFitMethodDd.Value;
            app.arcRansacMaxDist = app.ArcRansacMaxDistEf.Value;
            app.lineFitMethod = app.LineFitMethodDd.Value;
            app.lineFitMaxDist = app.LineFitMaxDistEf.Value;
            app.radius0 = app.Radius0Ef.Value;

            app.Msg = 'All the parameters are set. Click ''Plot'' to plot the data.';
            InfoTaValueChanged(app,true);
        end

        function S1ToolExtract2DLineBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Continue to fit the tool tip of 2D line format?', ...
                'Confirmation');
            switch selection
                case 'OK'
                    s1_toolExtract_2Dline;
                    app.S1ToolExtract2DLineBtn.Enable = 'off';
                    app.S1ToolExtract3DLineBtn.Enable = 'off';
                    app.S1ToolExtractSurfBtn.Enable = 'off';
                    app.S1ToolModelBtn.Enable = 'on';
                case 'Cancel'
                    return
            end
        end

        function S1ToolExtract3DLineBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Continue to fit the tool tip of 3D line format?', ...
                'Confirmation');
            switch selection
                case 'OK'
                    s1_toolExtract_3Dline;
                    app.S1ToolExtract2DLineBtn.Enable = 'off';
                    app.S1ToolExtract3DLineBtn.Enable = 'off';
                    app.S1ToolExtractSurfBtn.Enable = 'off';
                    app.S1ToolModelBtn.Enable = 'on';
                case 'Cancel'
                    return;
            end
        end

        function S1ToolExtractSurfBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Continue to fit the tool tip of surface format?', ...
                'Confirmation');
            switch selection
                case 'OK'
                    s1_toolExtract_surf;
                    app.S1ToolExtract2DLineBtn.Enable = 'off';
                    app.S1ToolExtract3DLineBtn.Enable = 'off';
                    app.S1ToolExtractSurfBtn.Enable = 'off';
                    app.S1ToolModelBtn.Enable = 'on';
                case 'Cancel'
                    return;
            end
        end

        % Value changed function: plot the tool file data
        function ToolfitPlotBtnPushed(app,event)
            ToolDataAxesClearBtnPushed(app,event);
            % plot the importing results
            plot(app.ToolDataAxes,app.toolFit(2,:),app.toolFit(3,:), ...
                'Color',[0,0.45,0.74],'LineWidth',0.75); % tool edge scatters
            hold(app.ToolDataAxes,'on');
            theta = (pi/2 - app.openAngle/2):0.01:(pi/2 + app.openAngle/2);
            xtmp = app.radius*cos(theta);
            ytmp = app.radius*sin(theta);
            plot(app.ToolDataAxes,xtmp,ytmp,'Color',[0.85,0.33,0.10], ...
                'LineWidth',1,'LineStyle','--'); % tool edge circle
            scatter(app.ToolDataAxes,0,0, ...
                'MarkerFaceColor',[0.85,0.33,0.10], ...
                'MarkerEdgeColor',[0.85,0.33,0.10]); % tool edge center
            grid(app.ToolDataAxes,'on');
            xlabel(app.ToolDataAxes,['x (',app.unit,')']);
            ylabel(app.ToolDataAxes,['y (',app.unit,')']);
            title(app.ToolDataAxes,'Tool Fit Results');
            app.Msg = 'Tool data is successfully loaded.';
            InfoTaValueChanged(app,true);
            clear theta xtmp ytmp;
        end


        function S1Tool2DBtnPushed(app,event)
            s1_tool2D;
        end

        function S1Tool3DBtnPushed(app,event)
            s1_tool3D;
        end

        function ToolinterpResetBtnPushed(app,event)
            app.ParamMethodDd.Value = app.paramMethodDefault;
            app.paramMethod = app.ParamMethodDd.Value;
        end

        function ToolinterpUpdateBtnPushed(app,event)
            app.paramMethod = app.ParamMethodDd.Value;
        end

        function S1ToolModelBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Continue to fit the tool tip of 2D line format?', ...
                'Confirmation');
            switch selection
                case 'OK'
                    s1_toolModel;
                    toolFitCheck(app);
                case 'Cancel'
                    return
            end
        end

        function ToolinterpPlotBtnPushed(app,event)
            ToolDataAxesClearBtnPushed(app,event);
            % plot the importing results
            plot(app.ToolDataAxes,app.toolFit(2,:),app.toolFit(3,:), ...
                '--.','MarkerSize',8,'Color',[0,0.447,0.741]);
            hold(app.ToolDataAxes,'on');
            plot(app.ToolDataAxes,app.toolData.toolCpts(2,:),app.toolData.toolCpts(3,:), ...
                'x','Color',[0.32,0.55,0.19],'MarkerSize',5);
            plot(app.ToolDataAxes,app.toolData.toolEdgePt(2,:),app.toolData.toolEdgePt(3,:), ...
                'Color',[0.635,0.078,0.184]);
            legend(app.ToolDataAxes,'Measured Pts','Control Pts', ...
                'Fitting Pts','Location','best');
            axis(app.ToolDataAxes,'on');
            % grid(app.ToolDataAxes,'on');
            xlabel(app.ToolDataAxes,['x (',app.unit,')']);
            ylabel(app.ToolDataAxes,['y (',app.unit,')']);
            title(app.ToolDataAxes,'Tool Interpolation Results');
            app.Msg = 'Tool data is successfully loaded.';
            InfoTaValueChanged(app,true);
        end

        function ToolDataAxesSaveMenuSelected(app,event)
            app.saveFig = app.ToolDataAxes;
            [figName,figFormat] = saveAxes2Fig(app);
            saveas(app.saveFig,figName,figFormat);
        end

        % Value changed function: cancel the process with nothing to be saved
        function ToolDataAxesClearBtnPushed(app,event)
            cla(app.ToolDataAxes,'reset');
            title(app.ToolDataAxes,'');
            app.CheckToolLamp.Color = 'r';
        end

        % Button down function: shift to the surface tab, and refresh the message
        function SurfaceTbButtonDown(app,event)
            app.Msg = ['Switch to the surface tab. ', ...
                'Press the corresponding button to finish the programming process.'];
            InfoTaValueChanged(app,true);
        end

        % Button down function: open the surface adding page
        function AddSurfaceBtnPushed(app,event)
            % Disable Plot Options button while dialog is open
            app.AddSurfaceBtn.Enable = 'off';

            % Open the add surface dialog
            app.surfType = '3D Geometry';
            app.AddSurfaceApp = add_surface(app,app.surfType,app.unit);
            app.surfPlotSpar = app.SurfacePlotSparSpin.Value;
            app.surfMesh = zeros(app.surfPlotSpar,app.surfPlotSpar,3);
        end

        function SurfacePlotSparSpinValueChanged(app,event)
            app.surfPlotSpar = app.SurfacePlotSparSpin.Value;
            app.surfMesh = zeros(app.surfPlotSpar,app.surfPlotSpar,3);

            app.Msg = ['The discretization of the surface plotting process', ...
                'is changed to ',num2str(app.surfPlotSpar),'.'];
            InfoTaValueChanged(app,true);
        end

        function SurfaceUnitDdValueChanged(app,event)
            app.surfaceUnit = app.SurfaceUnitDd.Value;

            % unit conversion (now: app.toolUnit, aim: app.unit)
            unitList = {'m','mm','\mum','nm'};
            presUnit = find(strcmp(unitList,app.surfaceUnit),1);
            aimUnit = find(strcmp(unitList,app.unit),1);
            % app.toolOri = 1000^(aimUnit - presUnit)*app.toolOri;
        end

        % Value changed function: plot the tool file data
        function Surface2DPlotBtnPushed(app,event)
            surfacePlot(app,'2D');
            app.Msg = 'Surface 2D image is successfully shown.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: plot the tool file data
        function Surface3DPlotBtnPushed(app,event)
            surfacePlot(app,'3D');
            app.Msg = 'Surface 3D image is successfully shown.';
            InfoTaValueChanged(app,true);
        end

        function SurfaceDataAxesSaveMenuSelected(app,event)
            app.saveFig = app.SurfaceDataAxes;
            [figName,figFormat] = saveAxes2Fig(app);
            saveas(app.saveFig,figName,figFormat);
        end

        % Value changed function: save the surface data
        function SurfaceSavedBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Save the surface data?','Comfirmation');
            switch selection
                case 'OK'
                    if ~isempty(app.surfFuncs)
                        if isempty(app.surfFx)
                            surfaceMeshGen(app);
                        end
                        app.Msg = 'Surface is corrected.';
                        InfoTaValueChanged(app,true);
                        app.CheckSurfLamp.Color = 'g';
                    else
                        uialert(app.UIFigure,'Error occurs in surface loading process!!!','Error load');
                        app.Msg = 'Error occurs in surface loading process!!!';
                        InfoTaValueChanged(app,true);
                        app.CheckSurfLamp.Color = 'r';
                    end
                    app.SurfaceSavedBtn.Enable = 'off';
                case 'Cancel'
                    return
            end
        end

        % Value changed function: cancel the process with nothing to be saved
        function SurfaceCancelBtnPushed(app,event)
            resetSurfaceParams(app);
            cla(app.SurfaceDataAxes,'reset');
            title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);

            app.Msg = 'Add-surface tab is reset.';
            InfoTaValueChanged(app,true);

            app.CheckSurfLamp.Color = 'r';
        end

        % Button down: shift to the program tab, and refresh the message
        function ProgramTbButtonDown(app,event)
            app.Msg = ['Switch to the program tab. ', ...
                'Press the corresponding button to finish the programming process.'];
            InfoTaValueChanged(app,true);
        end

        % Button down: shift to the optimization tab, and refresh the message
        function OptimTbButtonDown(app,event)
            app.Msg = ['Switch to the optimization tab. ', ...
                'Press the corresponding button to finish the programming process.'];
            InfoTaValueChanged(app,true);
        end

        % Dropdown : change the parameter list
        function AngularDiscreteDdValueChanged(app,event)
            switch app.AngularDiscreteDd.Value
                case 'Constant Angle'
                    app.ArcLengthLb.Visible = 'off';
                    app.ArcLengthEf.Visible = 'off';
                    app.ArcLengthUnitLb.Visible = 'off';
                    app.MaxAngPtDistLb.Visible = 'off';
                    app.MaxAngPtDistEf.Visible = 'off';
                    app.MaxAngPtDistUnitLb.Visible = 'off';

                    app.AngularLengthLb.Visible = 'on';
                    app.AngularLengthEf.Visible = 'on';
                    app.AngularLengthUnitLb.Visible = 'on';
                case 'Constant Arc'
                    app.ArcLengthLb.Visible = 'on';
                    app.ArcLengthEf.Visible = 'on';
                    app.ArcLengthUnitLb.Visible = 'on';
                    app.MaxAngPtDistLb.Visible = 'on';
                    app.MaxAngPtDistEf.Visible = 'on';
                    app.MaxAngPtDistUnitLb.Visible = 'on';

                    app.AngularLengthLb.Visible = 'off';
                    app.AngularLengthEf.Visible = 'off';
                    app.AngularLengthUnitLb.Visible = 'off';
            end
        end

        % Value changed function: update the selection selection
        function OptimUpdateBtnPushed(app,event)
            app.cutDirection = app.CutDirectionDd.Value;
            app.startDirection = app.StartDirectionDd.Value;
            app.angularDiscrete = app.AngularDiscreteDd.Value;
            app.radialIncrement = app.RadialIncrementDd.Value;
            app.aimRes = app.AimResEf.Value;
            app.maxIter = app.MaxIterSpin.Value;
            app.rStep = app.RStepEf.Value;
            app.arcLength = app.ArcLengthEf.Value;
            app.maxAngPtDist = app.MaxAngPtDistEf.Value;
            app.angularLength = app.AngularLengthEf.Value;
            app.spiralMethod = app.SpiralMethodDd.Value;
            app.zAllowance = app.ZAllowanceEf.Value;

            app.Msg = 'All the parameters are set.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: reset the parameters
        function OptimResetBtnPushed(app,event)
            % Reset all the values to the default
            resetOptimParams(app);

            % Report the infomation
            app.Msg = ['All the parameters in Optimization tab are reset.', ... 
                'Modify them and click ''Update'' to set.'];
            InfoTaValueChanged(app,true);
        end
        
        function Optim2DSingleBtnPushed(app,event)
            s4_optim_2D_aspheric_actual_tip_single;
            app.Msg = 'Toolpath optimization is successfully finished.';
            InfoTaValueChanged(app,true);
        end
        
        function Optim2DSingleIterBtnPushed(app,event)
            s4_optim_2D_aspheric_actual_tip_single_iter;
            app.Msg = 'Toolpath optimization is successfully finished.';
            InfoTaValueChanged(app,true);
        end
        
        function Optim2DMultiBtnPushed(app,event)
            s4_optim_2D_aspheric_actual_tip_multi;
            app.Msg = 'Toolpath optimization is successfully finished.';
            InfoTaValueChanged(app,true);
        end
        
        function Optim2DMultiIterBtnPushed(app,event)
            s4_optim_2D_aspheric_actual_tip_multi_iter;
            app.Msg = 'Toolpath optimization is successfully finished.';
            InfoTaValueChanged(app,true);
        end

        % --------------------------Function Execution--------------------------
        function S2DesignSimulAsphericConcentricBtnPushed(app,event)
            s2_design_simul_aspheric_concentric
        end

        function S2DesignSimulFreeformBtnPushed(app,event)
            s2_design_simul_freeform
        end

        % Code that update the infomation of the EfInfo window
        function InfoTaValueChanged(app,event)
            app.MsgNum = app.MsgNum + 1;
            app.InfoTa.Value{app.MsgNum} = char([num2str(app.MsgNum),' ',app.Msg]);
            scroll(app.InfoTa,"bottom");
        end

        % Cross pushed function: execute when pushing the cross of the UIFigure
        function UIFigureCloseReq(app,event)
            selection = uiconfirm(app.UIFigure,'Close the figure window?',...
                'Confirmation');
            switch selection
                case 'OK'
                    delete(app.UIFigure)
                case 'Cancel'
                    return
            end
        end
    end

    % Component initialization
    methods (Access = private)
        % Create UI figure and all the components
        function createComponents(app)
            % Create the figure panel and hide until all components are created
            app.UIFigure = uifigure('Name','Tool Data & Parameters Input', ...
                'WindowStyle','alwaysontop','WindowState','normal','Visible','off','SelectionType','extend');
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app,@UIFigureCloseReq,true);
            app.UIFigure.Resize = "on";
            app.UIFigure.Position = [500,500,800,700];
            app.UIFigure.Scrollable = "on";

            % Manage menu
            app.FigureMenuFiles = uimenu(app.UIFigure,'Text','Files', ...
                'Accelerator','F','Tooltip','File Operation');
            app.FigureMenuFilesNew = uimenu(app.FigureMenuFiles,'Text','New Window', ...
                'Accelerator','N','Tooltip','New File', ...
                'MenuSelectedFcn',createCallbackFcn(app,@FigureMenuFilesNewMenuSelect,true));
            app.FigureMenuFilesPrint = uimenu(app.FigureMenuFiles,'Text','Print', ...
                'Accelerator','S','Tooltip','Save Plot', ...
                'MenuSelectedFcn',createCallbackFcn(app,@FigureMenuFilesPrintMenuSelect,true));

            app.FigureMenuHelp = uimenu(app.UIFigure,'Text','Help', ...
                'Accelerator','H','Tooltip','Help for Software');
            app.FigureMenuHelpDoc = uimenu(app.FigureMenuHelp,'Text','Documentation', ...
                'Accelerator','F','Tooltip','Help Documentation', ...
                'MenuSelectedFcn',createCallbackFcn(app,@FigureMenuHelpFileMenuSelect,true));

            % Manage toolbar
            app.FigureToolbar = uitoolbar(app.UIFigure);

            app.AxisEqualPushtool = uipushtool(app.FigureToolbar, ...
                'Icon','@upm_toolpath_gui/Icon/AxisEqualPushtool.svg');
            app.AxisEqualPushtool.ClickedCallback = createCallbackFcn( ...
                app,@AxisEqualPushtoolClicked,true);
            app.AxisEqualPushtool.Tooltip = 'Set the axis of the current figure equal';

            app.CloseAllFigurePushtool = uipushtool(app.FigureToolbar, ...
                'Icon','@upm_toolpath_gui/Icon/CloseAllFigure.svg');
            app.CloseAllFigurePushtool.ClickedCallback = createCallbackFcn( ...
                app,@CloseAllFigurePushtoolClicked,true);
            app.CloseAllFigurePushtool.Tooltip = 'Close all the figures';

            app.BoldInfoToggletool = uitoggletool(app.FigureToolbar, ...
                'Icon','@upm_toolpath_gui/Icon/Bold.svg','Separator','on');
            app.BoldInfoToggletool.ClickedCallback = createCallbackFcn( ...
                app,@BoldInfoToggletoolClicked,true);
            app.BoldInfoToggletool.Tooltip = 'Bold/Normalize the infomation.';

            app.TopToggletool = uitoggletool(app.FigureToolbar, ...
                'Icon','@upm_toolpath_gui/Icon/top.svg');
            app.TopToggletool.ClickedCallback = createCallbackFcn( ...
                app,@TopToggletoolClicked,true);
            app.CloseAllFigurePushtool.Tooltip = 'Put the APP on top';

            % Manage tab groups
            app.FigureGl = uigridlayout(app.UIFigure,[3,1],'Padding',[5,5,5,5]);
            app.FigureGl.RowHeight = {'fit','1x','fit'};
            app.FigureGl.ColumnWidth = {'1x'};


            % ------------------------Workspace directory------------------------
            WorkspaceDirGl = uigridlayout(app.FigureGl,[2,7]);
            WorkspaceDirGl.Layout.Row = 1;
            WorkspaceDirGl.RowHeight = {'fit','fit'};
            WorkspaceDirGl.ColumnWidth = {'fit','2x','fit','fit','1x','fit','1x'};
            WorkspaceDirGl.Padding = [0,0,0,0];
            
            WorkspaceDirLb = uilabel(WorkspaceDirGl,'Text','Workspace directory:');
            WorkspaceDirLb.Layout.Row = 1;
            WorkspaceDirLb.Layout.Column = 1;
            
            app.WorkspaceDirEf = uieditfield(WorkspaceDirGl,'text');
            app.WorkspaceDirEf.Layout.Row = 2;
            app.WorkspaceDirEf.Layout.Column = [1,3];
            app.WorkspaceDirEf.ValueChangedFcn = createCallbackFcn(app,@WorkspaceDirEfValueChanged,true);
            
            app.WorkspaceDirBtn = uibutton(WorkspaceDirGl,'push','Text','Choose');
            app.WorkspaceDirBtn.Layout.Row = 1;
            app.WorkspaceDirBtn.Layout.Column = 3;
            app.WorkspaceDirBtn.ButtonPushedFcn = createCallbackFcn(app,@WorkspaceDirBtnPushed,true);
            
            %%%%%%%%% Create the font name selection label and dropdown box
            FontNameLb = uilabel(WorkspaceDirGl,'Text','Ploting font type');
            FontNameLb.Layout.Row = 1;
            FontNameLb.Layout.Column = 4;
            app.FontNameDd = uidropdown(WorkspaceDirGl,'Items',listfonts,'Value',app.fontNameDefault, ...
                'BackgroundColor',[1,1,1],'Visible','on');
            app.FontNameDd.Layout.Row = 1;
            app.FontNameDd.Layout.Column = 5;
            app.FontNameDd.ValueChangedFcn = createCallbackFcn(app,@FontNameDdValueChanged,true);
            
            %%%%%%%%% Create the font size selection label and editfield box
            FontSizeLb = uilabel(WorkspaceDirGl,'Text','Ploting font size');
            FontSizeLb.Layout.Row = 2;
            FontSizeLb.Layout.Column = 4;
            app.FontSizeEf = uieditfield(WorkspaceDirGl,'numeric','Limits',[6 72], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.fontSizeDefault);
            app.FontSizeEf.Layout.Row = 2;
            app.FontSizeEf.Layout.Column = 5;
            app.FontSizeEf.ValueChangedFcn = createCallbackFcn(app,@FontSizeEfValueChanged,true);
            
            %%%%%%%%% Create the unit selection label and dropdown box
            UnitLB = uilabel(WorkspaceDirGl,'Text','Project unit');
            UnitLB.Layout.Row = 1;
            UnitLB.Layout.Column = 6;
            app.UnitDd = uidropdown(WorkspaceDirGl,'Items',{'m','mm','\mum','nm'}, ...
                'Value',app.unitDefault,'BackgroundColor',[1,1,1]);
            app.UnitDd.Layout.Row = 1;
            app.UnitDd.Layout.Column = 7;
            app.UnitDd.ValueChangedFcn = createCallbackFcn(app,@UnitDdValueChanged,true);
            
            %%%%%%%%% ---Create the reset button---
            app.CommonResetBtn = uibutton(WorkspaceDirGl,'push','Text','Reset','Visible','on');
            app.CommonResetBtn.Layout.Row = 2;
            app.CommonResetBtn.Layout.Column = [6,7];
            app.CommonResetBtn.ButtonPushedFcn = createCallbackFcn(app,@CommonResetButtonPushed,true);

            % figure tab group
            app.FigureTbGp = uitabgroup(app.FigureGl,'SelectedTab',app.ToolTb);
            app.FigureTbGp.Layout.Row = 2;
            app.FigureTbGp.Layout.Column = 1;
            app.FigureTbGp.ContextMenu
            
            % ------------------------------------------------------------------------
            % --------------------------Tool File Processing--------------------------
            % ------------------------------------------------------------------------

            %%%%%%%%% tool bar object
            app.ToolTb = uitab(app.FigureTbGp,'Title','Tool Processing','Scrollable','on');
            app.ToolTb.ButtonDownFcn = createCallbackFcn(app,@ToolTbButtonDown,true);

            % Manage tool processing layout
            ToolTbGl = uigridlayout(app.ToolTb,[2,4]);
            ToolTbGl.RowHeight = {'1x','fit'};
            ToolTbGl.ColumnWidth = {'2x','1x','1x','1x'};
            % the elements of the cell array can be 'fit', fixed pixels, or '1x' '2x'
            %   'fit': to adjust the size to show the whole text, or basedon the default size
            %   fixed pixels: the size is fixed at the number of pixels
            %   variables ('2x'): to fill the remaining space, and the number is a
            %   weight for dividing up the remaining space among all variables

            % ------------------------Parameter editing tabgroup------------------------
            % Create the parameter selection tab group
            ParamTabgroup = uitabgroup(ToolTbGl,'SelectedTab',app.ParamTabToolImport);
            ParamTabgroup.Layout.Row = [1,2];
            ParamTabgroup.Layout.Column = 1;

            % ---------Data Import parameter selection---------
            app.ParamTabToolImport = uitab(ParamTabgroup,'Title','Data Import','Scrollable','on');
            ParamTabToolImportGl = uigridlayout(app.ParamTabToolImport,[5,3]);
            ParamTabToolImportGl.RowHeight = {'fit','2x','fit','fit','1x'};
            ParamTabToolImportGl.ColumnWidth = {'1x','1x','1x'};

            %%%%%%%%% Create the tool path selection 
            ToolFileLb = uilabel(ParamTabToolImportGl,'Text','Tool data path:');
            ToolFileLb.Layout.Row = 1;
            ToolFileLb.Layout.Column = 1;
            app.ToolFileEf = uieditfield(ParamTabToolImportGl,'text');
            app.ToolFileEf.Layout.Row = 2;
            app.ToolFileEf.Layout.Column = [1,3];
            app.ToolFileEf.ValueChangedFcn = createCallbackFcn(app,@ToolFileEfValueChanged,true);
            app.ToolFileBtn = uibutton(ParamTabToolImportGl,'push','Text','Choose');
            app.ToolFileBtn.Layout.Row = 3;
            app.ToolFileBtn.Layout.Column = 1;
            app.ToolFileBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolFileBtnPushed,true);

            % Create the tool data unit selection label and dropdown box
            ToolUnitLB = uilabel(ParamTabToolImportGl,'Text','Tool Data unit');
            ToolUnitLB.Layout.Row = 4;
            ToolUnitLB.Layout.Column = 1;
            app.ToolUnitDd = uidropdown(ParamTabToolImportGl,'Items',{'m','mm','\mum','nm'}, ...
                'Value',app.toolUnitDefault,'BackgroundColor',[1,1,1]);
            app.ToolUnitDd.Layout.Row = 4;
            app.ToolUnitDd.Layout.Column = [2,3];

            %%%%%%%%% Create the reset button for tool importing process
            app.ToolimportResetBtn = uibutton(ParamTabToolImportGl,'push','Text','Reset','Visible','on');
            app.ToolimportResetBtn.Layout.Row = 5;
            app.ToolimportResetBtn.Layout.Column = 3;
            app.ToolimportResetBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolImportResetBtnPushed,true);

            %%%%%%%%% Create the update button for tool importing process
            app.ToolimportUpdateBtn = uibutton(ParamTabToolImportGl,'push','Text','Update','Visible','on');
            app.ToolimportUpdateBtn.Layout.Row = 5;
            app.ToolimportUpdateBtn.Layout.Column = 1;
            app.ToolimportUpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolImportUpdateBtnPushed,true);

            %%%%%%%%% Create the plot button for tool importing process
            app.ToolimportPlotBtn = uibutton(ParamTabToolImportGl,'push','Text',{'Plot ','(Orig)'});
            app.ToolimportPlotBtn.Layout.Row = 5;
            app.ToolimportPlotBtn.Layout.Column = 2;
            app.ToolimportPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolImportPlotBtnPushed,true);

            % Create the cancle and clear button for the imported tool ploting
%             app.ToolImportCancelBtn = uibutton(ParamTabToolImportGl,'push','Text','Plot (Orig)');
%             app.ToolImportCancelBtn.Layout.Row = 3;
%             app.ToolImportCancelBtn.Layout.Column = 2;
%             app.ToolImportCancelBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolImportCancelBtnPushed,true);

            % ---------Tool fitting parameter selection---------
            app.ParamTabToolfit = uitab(ParamTabgroup,'Title','Tool fit','Scrollable','on');
            ParamTabToolfitGl = uigridlayout(app.ParamTabToolfit,[3,1]);
            ParamTabToolfitGl.RowHeight = {'1x','fit','fit'};
            ParamTabToolfitGl.ColumnWidth = {'1x'};

            ParamTabToolfitParamGl = uigridlayout(ParamTabToolfitGl,[7,2]);
            ParamTabToolfitParamGl.Layout.Row = 1;
            ParamTabToolfitParamGl.RowHeight = {'fit','fit','fit','fit','fit','fit','fit'};
            ParamTabToolfitParamGl.ColumnWidth = {'fit','1x'};
            ParamTabToolfitParamGl.Scrollable = 'on';

            %%%%%%%%% Create the tool fitting type label and dropdown box
            ToolFitTypeLb = uilabel(ParamTabToolfitParamGl,'Text','Tool fitting type');
            ToolFitTypeLb.Layout.Row = 1;
            ToolFitTypeLb.Layout.Column = 1;
            app.ToolFitTypeDd = uidropdown(ParamTabToolfitParamGl,'Items',{'onlyArc','arcRansac','lineArc'}, ...
                'Value',app.toolFitTypeDefault,'BackgroundColor',[1,1,1]);
            app.ToolFitTypeDd.Layout.Row = 1;
            app.ToolFitTypeDd.Layout.Column = 2;
            app.ToolFitTypeDd.ValueChangedFcn = createCallbackFcn(app,@ToolFitTypeDdValueChanged,true);
            app.ToolFitTypeDd.DropDownOpeningFcn = createCallbackFcn(app,@ToolFitTypeDdOpening,true);
            
            % Create the arc fitting method label and dropdown box
            ArcFitMethodLb = uilabel(ParamTabToolfitParamGl,'Text','Arc fitting method');
            ArcFitMethodLb.Layout.Row = 2;
            ArcFitMethodLb.Layout.Column = 1;
            app.ArcFitMethodDd = uidropdown(ParamTabToolfitParamGl, ...
                'Items',{'gradient-decent','normal-equation','levenberg-marquardt'}, ...
                'Value',app.arcFitMethodDefault);
            app.ArcFitMethodDd.Layout.Row = 2;
            app.ArcFitMethodDd.Layout.Column = 2;
            app.ArcFitMethodDd.BackgroundColor = [1 1 1];

            % Create the arc ransac fitting max distance label and dropdown box
            ArcRansacMaxDistLb = uilabel(ParamTabToolfitParamGl,'Text',{'Arc ransac fitting','max distance'});
            ArcRansacMaxDistLb.Layout.Row = 3;
            ArcRansacMaxDistLb.Layout.Column = 1;
            app.ArcRansacMaxDistEf = uieditfield(ParamTabToolfitParamGl,'numeric','Limits',[0 inf], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.arcRansacMaxDistDefault);
            app.ArcRansacMaxDistEf.Layout.Row = 3;
            app.ArcRansacMaxDistEf.Layout.Column = 2;
            app.ArcRansacMaxDistEf.Enable = 'off';
            app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
            % app.ArcRansacMaxDistEf.ValueChangedFcn = createCallbackFcn(app,@FontSizeEfValueChanged,true);

            % Create the line fitting method label and dropdown box
            LineFitMethodLb = uilabel(ParamTabToolfitParamGl,'Text','Line fitting method');
            LineFitMethodLb.Layout.Row = 4;
            LineFitMethodLb.Layout.Column = 1;
            app.LineFitMethodDd = uidropdown(ParamTabToolfitParamGl, ...
                'Items',{'polyfit','ransac'}, ...
                'Value',app.lineFitMethodDefault);
            app.LineFitMethodDd.Layout.Row = 4;
            app.LineFitMethodDd.Layout.Column = 2;
            app.LineFitMethodDd.Enable = "off";
            app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];

            % Create the line fitting max distance label and editfield box
            LineFitMaxDistEfLb = uilabel(ParamTabToolfitParamGl,'Text',{'Line fitting', 'max distance'});
            LineFitMaxDistEfLb.Layout.Row = 5;
            LineFitMaxDistEfLb.Layout.Column = 1;
            app.LineFitMaxDistEf = uieditfield(ParamTabToolfitParamGl,'numeric','Limits',[0 inf], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.lineFitMaxDistDefault);
            app.LineFitMaxDistEf.Layout.Row = 5;
            app.LineFitMaxDistEf.Layout.Column = 2;
            app.LineFitMaxDistEf.Enable = "off";
            app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
            
            Radius0Lb = uilabel(ParamTabToolfitParamGl,'Text',{'Theo Radius'});
            Radius0Lb.Layout.Row = 6;
            Radius0Lb.Layout.Column = 1;
            app.Radius0Ef = uieditfield(ParamTabToolfitParamGl,'numeric','Limits',[0 inf], ...
                'HorizontalAlignment','center','Value',app.radius0Default);
            app.Radius0Ef.Layout.Row = 6;
            app.Radius0Ef.Layout.Column = 2;

            %%%%%%%%% Create the reset button for tool fitting process
            app.ToolfitResetBtn = uibutton(ParamTabToolfitParamGl,'push','Text','Reset','Visible','on');
            app.ToolfitResetBtn.Layout.Row = 7;
            app.ToolfitResetBtn.Layout.Column = 1;
            app.ToolfitResetBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolFitResetBtnPushed,true);

            %%%%%%%%% Create the update button for tool fitting process
            app.ToolfitUpdateBtn = uibutton(ParamTabToolfitParamGl,'push','Text','Update','Visible','on');
            app.ToolfitUpdateBtn.Layout.Row = 7;
            app.ToolfitUpdateBtn.Layout.Column = 2;
            app.ToolfitUpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolFitUpdateBtnPushed,true);

            ParamTabToolfitExeGl = uigridlayout(ParamTabToolfitGl,[1,3]);
            ParamTabToolfitExeGl.Layout.Row = 2;
            ParamTabToolfitExeGl.RowHeight = {'fit'};
            ParamTabToolfitExeGl.ColumnWidth = {'1x','1x','1x'};
            ParamTabToolfitExeGl.Scrollable = 'on';

            %%%%%%%%% button to execute s1_toolExtract_2Dline.m
            app.S1ToolExtract2DLineBtn = uibutton(ParamTabToolfitExeGl,'push','WordWrap','on', ...
                'Text',{'Extract&Fit','(2D-Line)'},'FontWeight','bold');
            app.S1ToolExtract2DLineBtn.Layout.Column = 1;
            app.S1ToolExtract2DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract2DLineBtnPushed,true);

            %%%%%%%%% button to execute s1_toolExtract_3Dline.m
            app.S1ToolExtract3DLineBtn = uibutton(ParamTabToolfitExeGl,'push','WordWrap','on', ...
                'Text',{'Extract&Fit','(3D-Line)'});
            app.S1ToolExtract3DLineBtn.Layout.Column = 2;
            app.S1ToolExtract3DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract3DLineBtnPushed,true);

            %%%%%%%%% button to execute s1_toolExtract_surf.m
            app.S1ToolExtractSurfBtn = uibutton(ParamTabToolfitExeGl,'push','WordWrap','on', ...
                'Text',{'Extract&Fit','(3D-Surf)'});
            app.S1ToolExtractSurfBtn.Layout.Column = 3;
            app.S1ToolExtractSurfBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtractSurfBtnPushed,true);

            %%%%%%%%% button to plot the tool data
            app.ToolfitPlotBtn = uibutton(ParamTabToolfitGl,'push','Text','Plot (Fit)');
            app.ToolfitPlotBtn.Layout.Row = 3;
            app.ToolfitPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolfitPlotBtnPushed,true);

            % ---------Tool Edge Modification parameter selection---------
            app.ParamTabToolmod = uitab(ParamTabgroup,'Title','E.g. Fit','Scrollable','on');
            ParamTabToolmodGl = uigridlayout(app.ParamTabToolmod,[2,2]);
            ParamTabToolmodGl.RowHeight = {'fit','fit'};
            ParamTabToolmodGl.ColumnWidth = {'1x','1x'};
            ParamTabToolmodGl.Scrollable = 'on';

            % text
            toolModLb = uilabel(ParamTabToolmodGl,'Text','No parameters!');
            toolModLb.Layout.Row = 1;
            toolModLb.Layout.Column = [1,2];

            %%%%%%%%% button to execute s1_tool2D.m
            app.S1Tool2DBtn = uibutton(ParamTabToolmodGl,'push','WordWrap','on', ...
                'Text',{'Example Fit','(2D Data)'});
            app.S1Tool2DBtn.Layout.Row = 2;
            app.S1Tool2DBtn.Layout.Column = 1;
            app.S1Tool2DBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1Tool2DBtnPushed,true);

            %%%%%%%%% button to execute s1_tool3D.m
            app.S1Tool3DBtn = uibutton(ParamTabToolmodGl,'push','WordWrap','on', ...
                'Text',{'Example Fit','(3D Data)'});
            app.S1Tool3DBtn.Layout.Row = 2;
            app.S1Tool3DBtn.Layout.Column = 2;
            app.S1Tool3DBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1Tool3DBtnPushed,true);

            % ---------Tool Interpolation parameter selection---------
            app.ParamTabToolinterp = uitab(ParamTabgroup,'Title','Tool Interpolation','Scrollable','on');
            ParamTabToolinterpGl = uigridlayout(app.ParamTabToolinterp,[4,2]);
            ParamTabToolinterpGl.RowHeight = {'fit','fit','fit','fit'};
            ParamTabToolinterpGl.ColumnWidth = {'1x','1x'};
            ParamTabToolinterpGl.Scrollable = 'on';

            % Create the B-spline parametric method label and dropdown box
            ParamMethodLb = uilabel(ParamTabToolinterpGl,'Text',{'B-spline', 'param-method'});
            ParamMethodLb.Layout.Row = 1;
            ParamMethodLb.Layout.Column = 1;
            app.ParamMethodDd = uidropdown(ParamTabToolinterpGl, ...
                'Items',{'uniform','chord','centripetal'}, ...
                'Value',app.paramMethodDefault,'BackgroundColor',[1,1,1]);
            app.ParamMethodDd.Layout.Row = 1;
            app.ParamMethodDd.Layout.Column = 2;

            %%%%%%%%% Create the reset button
            app.ToolinterpResetBtn = uibutton(ParamTabToolinterpGl,'push','Text','Reset','Visible','on');
            app.ToolinterpResetBtn.Layout.Row = 2;
            app.ToolinterpResetBtn.Layout.Column = 1;
            app.ToolinterpResetBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolinterpResetBtnPushed,true);

            %%%%%%%%% Create the update button
            app.ToolinterpUpdateBtn = uibutton(ParamTabToolinterpGl,'push','Text','Update','Visible','on');
            app.ToolinterpUpdateBtn.Layout.Row = 2;
            app.ToolinterpUpdateBtn.Layout.Column = 2;
            app.ToolinterpUpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolinterpUpdateBtnPushed,true);

            % button to execute s1_toolModel.m
            app.S1ToolModelBtn = uibutton(ParamTabToolinterpGl,'push','WordWrap','on', ...
                'Text','Interpolation','Enable','off');
            app.S1ToolModelBtn.Layout.Row = 3;
            app.S1ToolModelBtn.Layout.Column = [1,2];
            app.S1ToolModelBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolModelBtnPushed,true);

            % button to plot the interpolation result
            app.ToolinterpPlotBtn = uibutton(ParamTabToolinterpGl,'Text','Plot Interpolation');
            app.ToolinterpPlotBtn.Layout.Row = 4;
            app.ToolinterpPlotBtn.Layout.Column = [1,2];
            app.ToolinterpPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolinterpPlotBtnPushed,true);


            % ------------------------Ploting axes------------------------
            app.ToolDataAxes = uiaxes(ToolTbGl);
            title(app.ToolDataAxes,'tool original data','FontSize',16);
            % xlabel(app.UIAxes, 'X');
            % ylabel(app.UIAxes, 'Y');
            % zlabel(app.UIAxes, 'Z');
            app.ToolDataAxes.Layout.Row = 1;
            app.ToolDataAxes.Layout.Column = [2,4];
            app.ToolDataAxes.ContextMenu = uicontextmenu(app.UIFigure);
            ToolDataAxesMenu1 = uimenu(app.ToolDataAxes.ContextMenu,'Text','Save', ...
                'MenuSelectedFcn',createCallbackFcn(app,@ToolDataAxesSaveMenuSelected,true));

            app.ToolDataAxesClearBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Cancel');
            app.ToolDataAxesClearBtn.Layout.Row = 2;
            app.ToolDataAxesClearBtn.Layout.Column = 4;
            app.ToolDataAxesClearBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolDataAxesClearBtnPushed,true);

            % ------------------------------------------------------------------------
            % ---------------------------Surface Processing---------------------------
            % ------------------------------------------------------------------------

            app.SurfaceTb = uitab(app.FigureTbGp,'Title','Surface Load','Scrollable','on');
            app.SurfaceTb.ButtonDownFcn = createCallbackFcn(app,@SurfaceTbButtonDown,true);

            SurfaceTbGl = uigridlayout(app.SurfaceTb,[6,4]);
            SurfaceTbGl.RowHeight = {'fit','fit','fit','1x','fit','fit'};
            SurfaceTbGl.ColumnWidth = {'2x','2x','3x','3x'};

            app.AddSurfaceBtn = uibutton(SurfaceTbGl,'push','WordWrap','on', ...
                'Text','Add Geometry');
            app.AddSurfaceBtn.Layout.Row = 1;
            app.AddSurfaceBtn.Layout.Column = [1,2];
            app.AddSurfaceBtn.Icon = '@upm_toolpath_gui/Icon/AddSurf1.svg';
            app.AddSurfaceBtn.ButtonPushedFcn = createCallbackFcn(app,@AddSurfaceBtnPushed,true);

            SurfacePlotSparLb = uilabel(SurfaceTbGl,'Text','Plot Discretization:');
            SurfacePlotSparLb.Layout.Row = 2;
            SurfacePlotSparLb.Layout.Column = 1;
            app.SurfacePlotSparSpin = uispinner(SurfaceTbGl,'Value',app.surfPlotSparDefault, ...
                'Limits',[0 Inf]);
            app.SurfacePlotSparSpin.Layout.Row = 2;
            app.SurfacePlotSparSpin.Layout.Column = 2;

            SurfaceUnitLb = uilabel(SurfaceTbGl,'Text','Surface Unit');
            SurfaceUnitLb.Layout.Row = 3;
            SurfaceUnitLb.Layout.Column = 1;
            app.SurfaceUnitDd = uidropdown(SurfaceTbGl);
            app.SurfaceUnitDd.ValueChangedFcn = createCallbackFcn(app,@SurfaceUnitDdValueChanged,true);
            app.SurfaceUnitDd.Layout.Row = 3;
            app.SurfaceUnitDd.Layout.Column = 2;

            app.SurfaceDetailTa = uitextarea(SurfaceTbGl,'WordWrap','on', ...
                'Editable','off');
            app.SurfaceDetailTa.Layout.Row = [4,6];
            app.SurfaceDetailTa.Layout.Column = [1,2];
            app.SurfaceDetailTa.FontName = 'Times New Roman';
            app.SurfaceDetailTa.FontSize = 16;
            app.SurfaceDetailTa.FontWeight = 'normal';
            
            app.SurfaceDataAxes = uiaxes(SurfaceTbGl);
            title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            app.SurfaceDataAxes.Layout.Row = [1,4];
            app.SurfaceDataAxes.Layout.Column = [3,4];
            app.SurfaceDataAxes.ContextMenu = uicontextmenu(app.UIFigure);
            SurfaceDataAxesMenu1 = uimenu(app.SurfaceDataAxes.ContextMenu,'Text','Save', ...
                'MenuSelectedFcn',createCallbackFcn(app,@SurfaceDataAxesSaveMenuSelected,true));

            app.Surface2DPlotBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','2D Image');
            app.Surface2DPlotBtn.Layout.Row = 5;
            app.Surface2DPlotBtn.Layout.Column = 3;
            app.Surface2DPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@Surface2DPlotBtnPushed,true);

            app.Surface3DPlotBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','3D Image');
            app.Surface3DPlotBtn.Layout.Row = 5;
            app.Surface3DPlotBtn.Layout.Column = 4;
            app.Surface3DPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@Surface3DPlotBtnPushed,true);

            % ------------------------Ending process------------------------
            app.SurfaceSavedBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','Enter');
            app.SurfaceSavedBtn.Layout.Row = 6;
            app.SurfaceSavedBtn.Layout.Column = 3;
            app.SurfaceSavedBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfaceSavedBtnPushed,true);
            
            app.SurfaceCancelBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','Cancel');
            app.SurfaceCancelBtn.Layout.Row = 6;
            app.SurfaceCancelBtn.Layout.Column = 4;
            app.SurfaceCancelBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfaceCancelBtnPushed,true);

            % ------------------------------------------------------------------------
            % ------------------------- Optimization Process -------------------------
            % ------------------------------------------------------------------------

            app.OptimTb = uitab(app.FigureTbGp,'Title','Optimization');
            app.OptimTb.ButtonDownFcn = createCallbackFcn(app,@OptimTbButtonDown,true);
            app.OptimTb.Scrollable = 'on';

            OptimTbGl = uigridlayout(app.OptimTb,[2,2]);
            OptimTbGl.RowHeight = {'fit','fit'};
            OptimTbGl.ColumnWidth = {'1x','1x'};
            OptimTbGl.Scrollable = 'on';

            % ------------------------- Optim- Condition Check Area -------------------------
            CheckGl = uigridlayout(OptimTbGl,[1,2],'Padding',[0,0,0,0]);
            CheckGl.Layout.Row = 1;
            CheckGl.Layout.Column = [1,2];
            CheckGl.RowHeight = {'fit'};
            CheckGl.ColumnWidth = {'1x','1x'};

            CheckToolPn = uipanel(CheckGl);
            CheckToolPn.Layout.Row = 1;
            CheckToolPn.Layout.Column = 1;
            CheckToolPnGl = uigridlayout(CheckToolPn,[1,3]);
            CheckToolPnGl.RowHeight = {'fit'};
            CheckToolPnGl.ColumnWidth = {'fit','1x','fit'};
            CheckToolLb = uilabel(CheckToolPnGl,'Text','Check Tool');
            CheckToolLb.Layout.Row = 1;
            CheckToolLb.Layout.Column = 1;
            app.CheckToolLamp = uilamp(CheckToolPnGl,'Color','r');
            app.CheckToolLamp.Layout.Row = 1;
            app.CheckToolLamp.Layout.Column = 2;
%             app.CheckToolPath = uieditfield()

            CheckSurfPn = uipanel(CheckGl);
            CheckSurfPn.Layout.Row = 1;
            CheckSurfPn.Layout.Column = 2;
            CheckSurfPnGl = uigridlayout(CheckSurfPn,[1,2]);
            CheckSurfPnGl.RowHeight = {'fit'};
            CheckSurfPnGl.ColumnWidth = {'fit','1x'};
            CheckSurfLb = uilabel(CheckSurfPnGl,'Text','Check Surface');
            CheckSurfLb.Layout.Row = 1;
            CheckSurfLb.Layout.Column = 1;
            app.CheckSurfLamp = uilamp(CheckSurfPnGl,'Color','r');
            app.CheckSurfLamp.Layout.Row = 1;
            app.CheckSurfLamp.Layout.Column = 2;

            % ------------------------- Optim- Parameters Selection -------------------------
            OptimParamPn = uipanel(OptimTbGl,'Title','Edit Parameters', ...
                'TitlePosition','centertop');
            OptimParamPn.Layout.Row = 2;
            OptimParamPn.Layout.Column = 1;
            OptimParamPn.Scrollable = 'on';
            OptimParamPnGl = uigridlayout(OptimParamPn,[2,2]);
            OptimParamPnGl.RowHeight = {'1x','fit'};
            OptimParamPnGl.ColumnWidth = {'1x','1x'};
            OptimParamPnGl.Scrollable = 'on';
            
            OptimParamTabgroup = uitabgroup(OptimParamPnGl,'SelectedTab',app.OptimParamPathTab);
            OptimParamTabgroup.Layout.Row = 1;
            OptimParamTabgroup.Layout.Column = [1,2];

            app.OptimParamPathTab = uitab(OptimParamTabgroup,'Title','Tool Path');
            app.OptimParamPathTab.Scrollable = 'on';
            OptimParamPathTabGl = uigridlayout(app.OptimParamPathTab,[9,3]);
            OptimParamPathTabGl.Scrollable = 'on';
            OptimParamPathTabGl.RowHeight = {'fit','fit','fit','fit','fit','fit','fit','fit','fit'};
            OptimParamPathTabGl.ColumnWidth = {'fit','1x','fit'};

            CutDirectionLb = uilabel(OptimParamPathTabGl,'Text','Cut Direction');
            CutDirectionLb.Layout.Row = 1;
            CutDirectionLb.Layout.Column = 1;
            app.CutDirectionDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'Center to Edge','Edge to Center'}, ...
                'Value',app.cutDirectionDefault);
            app.CutDirectionDd.Layout.Row = 1;
            app.CutDirectionDd.Layout.Column = 2;

            SpindleDirectionLb = uilabel(OptimParamPathTabGl,'Text','Spindle Direction');
            SpindleDirectionLb.Layout.Row = 2;
            SpindleDirectionLb.Layout.Column = 1;
            app.StartDirectionDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'X Plus','X Minus'}, ...
                'Value',app.startDirectionDefault);
            app.StartDirectionDd.Layout.Row = 2;
            app.StartDirectionDd.Layout.Column = 2;

            RadialIncrementLb = uilabel(OptimParamPathTabGl,'Text','Radial Increment Type');
            RadialIncrementLb.Layout.Row = 3;
            RadialIncrementLb.Layout.Column = 1;
            app.RadialIncrementDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'On-Axis','Surface'}, ...
                'Value',app.radialIncrementDefault);
            app.RadialIncrementDd.Layout.Row = 3;
            app.RadialIncrementDd.Layout.Column = 2;
            app.RadialIncrementDd.Enable = 'off';

            AimResLb = uilabel(OptimParamPathTabGl,'Text','Aimed Residual Height');
            AimResLb.Layout.Row = 4;
            AimResLb.Layout.Column = 1;
            app.AimResEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.aimResDefault);
            app.AimResEf.Layout.Row = 4;
            app.AimResEf.Layout.Column = 2;
            AimResUnitLb = uilabel(OptimParamPathTabGl,'Interpreter','latex');
            AimResUnitLb.Layout.Row = 4;
            AimResUnitLb.Layout.Column = 3;
            AimResUnitLb.Text = app.unit;

            AngularDiscreteLb = uilabel(OptimParamPathTabGl,'Text','Angular Increment Type');
            AngularDiscreteLb.Layout.Row = 5;
            AngularDiscreteLb.Layout.Column = 1;
            app.AngularDiscreteDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'Constant Arc','Constant Angle'}, ...
                'Value',app.angularDiscreteDefault);
            app.AngularDiscreteDd.Layout.Row = 5;
            app.AngularDiscreteDd.Layout.Column = 2;
            app.AngularDiscreteDd.ValueChangedFcn = createCallbackFcn(app,@AngularDiscreteDdValueChanged,true);

            app.ArcLengthLb = uilabel(OptimParamPathTabGl,'Text','Arc Length');
            app.ArcLengthLb.Layout.Row = 6;
            app.ArcLengthLb.Layout.Column = 1;
            app.ArcLengthEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.arcLengthDefault);
            app.ArcLengthEf.Layout.Row = 6;
            app.ArcLengthEf.Layout.Column = 2;
            app.ArcLengthUnitLb = uilabel(OptimParamPathTabGl,'Interpreter','latex');
            app.ArcLengthUnitLb.Layout.Row = 6;
            app.ArcLengthUnitLb.Layout.Column = 3;
            app.ArcLengthUnitLb.Text = app.unit;

            app.MaxAngPtDistLb = uilabel(OptimParamPathTabGl,'Text','Max Angular Point Distance');
            app.MaxAngPtDistLb.Layout.Row = 7;
            app.MaxAngPtDistLb.Layout.Column = 1;
            app.MaxAngPtDistEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.maxAngPtDistDefault);
            app.MaxAngPtDistEf.Layout.Row = 7;
            app.MaxAngPtDistEf.Layout.Column = 2;
            app.MaxAngPtDistUnitLb = uilabel(OptimParamPathTabGl,'Interpreter','latex');
            app.MaxAngPtDistUnitLb.Layout.Row = 7;
            app.MaxAngPtDistUnitLb.Layout.Column = 3;
            app.MaxAngPtDistUnitLb.Text = 'rad';

            app.AngularLengthLb = uilabel(OptimParamPathTabGl,'Text','Angular Length');
            app.AngularLengthLb.Layout.Row = 6;
            app.AngularLengthLb.Layout.Column = 1;
            app.AngularLengthLb.Visible = 'off';
            app.AngularLengthEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.angularLengthDefault);
            app.AngularLengthEf.Layout.Row = 6;
            app.AngularLengthEf.Layout.Column = 2;
            app.AngularLengthEf.Visible = 'off';
            app.AngularLengthUnitLb = uilabel(OptimParamPathTabGl,'Text','rad');
            app.AngularLengthUnitLb.Layout.Row = 6;
            app.AngularLengthUnitLb.Layout.Column = 3;
            app.AngularLengthUnitLb.Visible = 'off';

            MaxIterLb = uilabel(OptimParamPathTabGl,'Text','Max No. of Iteration');
            MaxIterLb.Layout.Row = 8;
            MaxIterLb.Layout.Column = 1;
            app.MaxIterSpin = uispinner(OptimParamPathTabGl,'Value',app.maxIterDefault, ...
                'Limits',[1,inf],'ValueDisplayFormat','%d','RoundFractionalValues','on');
            app.MaxIterSpin.Layout.Row = 8;
            app.MaxIterSpin.Layout.Column = 2;
            app.MaxIterSpin.Enable = 'off';

            % initial annulus width when calculating concentric tool path
            RStepLb = uilabel(OptimParamPathTabGl,'Text','Initial Annulus Width (rStep)');
            RStepLb.Layout.Row = 9;
            RStepLb.Layout.Column = 1;
            app.RStepEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center');
            app.RStepEf.Layout.Row = 9;
            app.RStepEf.Layout.Column = 2;
            RStepUnitLb = uilabel(OptimParamPathTabGl,'Interpreter','latex');
            RStepUnitLb.Layout.Row = 9;
            RStepUnitLb.Layout.Column = 3;
            RStepUnitLb.Text = app.unit;

            ZAllowanceLb = uilabel(OptimParamPathTabGl,'Text','Surface Domain Allowance');
            ZAllowanceLb.Layout.Row = 10;
            ZAllowanceLb.Layout.Column = 1;
            app.ZAllowanceEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.zAllowanceDefault);
            app.ZAllowanceEf.Layout.Row = 10;
            app.ZAllowanceEf.Layout.Column = 2;

            % feed tab
            app.OptimParamFeedTab = uitab(OptimParamTabgroup,'Title','Feed');
            app.OptimParamFeedTab.Scrollable = 'on';
            OptimParamFeedTabGl = uigridlayout(app.OptimParamFeedTab,[1,3]);
            OptimParamFeedTabGl.Scrollable = 'on';
            OptimParamFeedTabGl.RowHeight = {'fit'};
            OptimParamFeedTabGl.ColumnWidth = {'fit','1x','fit'};

            SpiralMethodLb = uilabel(OptimParamFeedTabGl,'Text','Spiral Method');
            SpiralMethodLb.Layout.Row = 1;
            SpiralMethodLb.Layout.Column = 1;
            app.SpiralMethodDd = uidropdown(OptimParamFeedTabGl, ...
                'Items',{'Radius-Number','Radius-Angle'}, ...
                'Value',app.spiralMethodDefault);
            app.SpiralMethodDd.Layout.Row = 1;
            app.SpiralMethodDd.Layout.Column = 2;

            % ---Create the reset button---
            app.OptimResetBtn = uibutton(OptimParamPnGl,'push','Text','Reset','Visible','on');
            app.OptimResetBtn.Layout.Row = 2;
            app.OptimResetBtn.Layout.Column = 1;
            app.OptimResetBtn.ButtonPushedFcn = createCallbackFcn(app,@OptimResetBtnPushed,true);
            
            % ---Create the update button---
            app.OptimUpdateBtn = uibutton(OptimParamPnGl,'push','Text','Update','Visible','on');
            app.OptimUpdateBtn.Layout.Row = 2;
            app.OptimUpdateBtn.Layout.Column = 2;
            app.OptimUpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@OptimUpdateBtnPushed,true);

            % --------------
            OptimProcessGl = uigridlayout(OptimTbGl,[2,2]);
            OptimProcessGl.Scrollable = 'on';
            OptimProcessGl.RowHeight = {'1x','1x'};
            OptimProcessGl.ColumnWidth = {'1x','1x'};
            OptimProcessGl.Layout.Row = 2;
            OptimProcessGl.Layout.Column = 2;

            app.Optim2DSingleBtn = uibutton(OptimProcessGl,'push','WordWrap','on', ...
                'Text',{'tool path optimization';'2D aspheric surface';'single-interaction & solve'}, ...
                'Icon','@upm_toolpath_gui/Icon/s4_optim_2D_aspheric_actual_tip_single.svg', ...
                'IconAlignment','top');
            app.Optim2DSingleBtn.Layout.Row = 1;
            app.Optim2DSingleBtn.Layout.Column = 1;
            app.Optim2DSingleBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@Optim2DSingleBtnPushed,true);

            app.Optim2DSingleIterBtn = uibutton(OptimProcessGl,'push','WordWrap','on', ...
                'Text',{'tool path optimization';'2D aspheric surface';'single-interaction & iter'}, ...
                'Icon','@upm_toolpath_gui/Icon/s4_optim_2D_aspheric_actual_tip_single_iter.svg', ...
                'IconAlignment','top');
            app.Optim2DSingleIterBtn.Layout.Row = 1;
            app.Optim2DSingleIterBtn.Layout.Column = 2;
            app.Optim2DSingleIterBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@Optim2DSingleIterBtnPushed,true);

            app.Optim2DMultiBtn = uibutton(OptimProcessGl,'push','WordWrap','on', ...
                'Text',{'tool path optimization';'2D aspheric surface';'multi-interaction & solve'}, ...
                'Icon','@upm_toolpath_gui/Icon/s4_optim_2D_aspheric_actual_tip_multi.svg', ...
                'IconAlignment','top');
            app.Optim2DMultiBtn.Layout.Row = 2;
            app.Optim2DMultiBtn.Layout.Column = 1;
            app.Optim2DMultiBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@Optim2DMultiBtnPushed,true);

            app.Optim2DMultiIterBtn = uibutton(OptimProcessGl,'push','WordWrap','on', ...
                'Text',{'tool path optimization';'2D aspheric surface';'multi-interaction & iter'}, ...
                'Icon','@upm_toolpath_gui/Icon/s4_optim_2D_aspheric_actual_tip_multi_iter.svg', ...
                'IconAlignment','top');
            app.Optim2DMultiIterBtn.Layout.Row = 2;
            app.Optim2DMultiIterBtn.Layout.Column = 2;
            app.Optim2DMultiIterBtn.FontWeight = 'bold';
            app.Optim2DMultiIterBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@Optim2DMultiIterBtnPushed,true);

            % ------------------------------------------------------------------------
            % ---------------------------Program Processing---------------------------
            % ------------------------------------------------------------------------

            app.ProgramTb = uitab(app.FigureTbGp,'Title','Program','Scrollable','on');
            app.ProgramTb.ButtonDownFcn = createCallbackFcn(app,@ProgramTbButtonDown,true);
            app.ProgramTb.HandleVisibility = 'off';
            ProgramGl = uigridlayout(app.ProgramTb,[4,2]);
            ProgramGl.RowHeight = {'1x','1x','1x','fit'};
            ProgramGl.ColumnWidth = {'1x','1x'};

            % ---------------------------simulation process---------------------------
            SimulatePn = uipanel(ProgramGl,'Visible','on', ...
                'Title','','TitlePosition','centertop');
            SimulatePn.Layout.Row = 1;
            SimulatePn.Layout.Column = 2;
            SimulatePn.Scrollable = 'on';

            SimulateGl = uigridlayout(SimulatePn,[1,2]);
            SimulateGl.RowHeight = {'1x'};
            SimulateGl.ColumnWidth = {'1x','1x'};
            
            % button to execute s2_design_simul_aspheric_concentric.m
            app.S2DesignSimulAsphericConcentricBtn = uibutton(SimulateGl,'push', ...
                'WordWrap','on','Text','s2_design_simul_aspheric_concentric');
            app.S2DesignSimulAsphericConcentricBtn.Layout.Column = 1;
            app.S2DesignSimulAsphericConcentricBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S2DesignSimulAsphericConcentricBtnPushed,true);

            % button to execute s2_design_simul_freeform.m
            app.S2DesignSimulFreeformBtn = uibutton(SimulateGl,'push', ...
                'WordWrap','on','Text','s2_design_simul_freeform');
            app.S2DesignSimulFreeformBtn.Layout.Column = 2;
            app.S2DesignSimulFreeformBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S2DesignSimulFreeformBtnPushed,true);

            % concentric optimization process
            ConcentricOptimBtngp = uibuttongroup(ProgramGl);
            ConcentricOptimBtngp.Layout.Row = 2;
            ConcentricOptimBtngp.Layout.Column = 2;

            % ---------------------------Optimization process---------------------------
            OptimPn = uipanel(ProgramGl,'Visible','on', ...
                'Title','','TitlePosition','centertop');
            OptimPn.Layout.Row = 2;
            OptimPn.Layout.Column = 2;
            OptimPn.Scrollable = 'on';

            OptimGl = uigridlayout(OptimPn,[1,5]);
            OptimGl.RowHeight = {'1x'};
            OptimGl.ColumnWidth = {'fit','1x','fit','1x','fit'};

            OptimArrow1 = uiimage(OptimGl,'ImageSource','@upm_toolpath_gui/Icon/RightArrow6.svg', ...
                'ScaleMethod','stretch');
            OptimArrow1.Layout.Row = 1;
            OptimArrow1.Layout.Column = 2;

            OptimArrow2 = uiimage(OptimGl,'ImageSource','@upm_toolpath_gui/Icon/RightArrow6.svg', ...
                'ScaleMethod','stretch');
            OptimArrow2.Layout.Row = 1;
            OptimArrow2.Layout.Column = 4;

            % ------------------------------------------------------------------------
            % -------------------------Info displaying window-------------------------
            % ------------------------------------------------------------------------
            app.InfoTa = uitextarea(app.FigureGl,'WordWrap','on', ...
                'Editable','off','BackgroundColor',[0.96 0.96 0.96]);
            app.InfoTa.Layout.Row = 3;
            app.InfoTa.Layout.Column = 1;
            app.InfoTa.FontSize = 14;
            app.InfoTa.FontWeight = 'normal';
            app.InfoTa.Value = {'1 Welcome.'};
            app.InfoTa.ValueChangedFcn = createCallbackFcn(app,@InfoTaValueChanged,true);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            app.UIFigure.WindowStyle = "normal";
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = upm_toolpath_gui
            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app,app.UIFigure);

            % Execute the startup function
            runStartupFcn(app,@startupFcn);

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
            % rmpath(genpath('funcs');
        end
    end
end
