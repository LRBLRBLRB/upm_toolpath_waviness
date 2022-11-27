classdef upm_toolpath < matlab.apps.AppBase
    % to create a ui figure for tool data & parameters input

    % Default properties
    properties (Constant)
        workspaceDirDefault     = 'D:\Code\2021-11_ToolWaviness\upm_toolpath_waviness\workspace';
        unitDefault             = 'mm'
        fontNameDefault         = 'Times New Roman'
        fontSizeDefault         = 12
        toolFitTypeDefault      = 'lineArc'
        arcFitMethodDefault     = 'levenberg-marquardt'
        arcRansacMaxDistDefault = 1e-2
        lineFitMethodDefault    = 'polyfit'
        lineFitMaxDistDefault   = 1e-3
        paramMethodDefault      = 'centripetal'
        surfFuncsDefault        = 'C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2)';
        cutDirectionDefault     = 'Center to Edge'
        spindleDirectionDefault = 'Counterclockwise'
        angularDiscreteDefault  = 'Constant Arc'
        aimResDefault           = 50
        rStepDefault            = 200
        arcLengthDefault        = 30
        maxAngPtDistDefault     = 6*pi/180
        angularLengthDefault    = 6*pi/180
    end

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        FigureMenu              matlab.ui.container.Menu
        FigureToolbar           matlab.ui.container.Toolbar
        CloseAllFigurePushtool  matlab.ui.container.toolbar.PushTool
        BoldInfoToggletool      matlab.ui.container.toolbar.ToggleTool
        TopToggletool           matlab.ui.container.toolbar.ToggleTool
        FigureTbGp              matlab.ui.container.TabGroup
        ToolTb                  matlab.ui.container.Tab
        WorkspaceDirEf          matlab.ui.control.EditField
        WorkspaceDirBtn         matlab.ui.control.Button
        ToolFileEf              matlab.ui.control.EditField
        ToolFileBtn             matlab.ui.control.Button
        ParamTabCommon          matlab.ui.container.Tab
        UnitDd                  matlab.ui.control.DropDown
        FontNameDd              matlab.ui.control.DropDown
        FontSizeEf              matlab.ui.control.NumericEditField
        ParamTabToolfit         matlab.ui.container.Tab
        ToolFitTypeDd           matlab.ui.control.DropDown
        LineFitMethodDd         matlab.ui.control.DropDown
        LineFitMaxDistEf        matlab.ui.control.NumericEditField
        ArcFitMethodDd          matlab.ui.control.DropDown
        ArcRansacMaxDistEf      matlab.ui.control.NumericEditField
        ParamTabToolinterp      matlab.ui.container.Tab
        ParamMethodDd           matlab.ui.control.DropDown
        ToolResetBtn            matlab.ui.control.Button
        ToolUpdateBtn           matlab.ui.control.Button
        ToolDataAxes            matlab.ui.control.UIAxes
        ToolPlotBtn             matlab.ui.control.Button
        Tool2DLineBtnBtn        matlab.ui.control.Button
        Tool3DLineBtnBtn        matlab.ui.control.Button
        ToolCancelBtn           matlab.ui.control.Button
        SurfaceTb               matlab.ui.container.Tab
        AddSurfaceBtn           matlab.ui.control.Button
        SurfaceDetailTa         matlab.ui.control.TextArea
        SurfaceDataAxes         matlab.ui.control.UIAxes
        SurfacePlotBtn          matlab.ui.control.Button
        SurfaceSavedBtn         matlab.ui.control.Button
        SurfaceCancelBtn        matlab.ui.control.Button
        ProgramTb               matlab.ui.container.Tab
        OptimTb                 matlab.ui.container.Tab
        CheckToolLamp           matlab.ui.control.Lamp
        CheckSurfLamp           matlab.ui.control.Lamp
        OptimParamPathTab       matlab.ui.container.Tab
        CutDirectionDd          matlab.ui.control.DropDown
        SpindleDirectionDd      matlab.ui.control.DropDown
        AngularDiscreteDd       matlab.ui.control.DropDown
        AimResEf                matlab.ui.control.NumericEditField
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
        OptimParamFeedTab       matlab.ui.container.Tab
        OptimUpdateBtn          matlab.ui.control.Button
        OptimResetBtn           matlab.ui.control.Button
        InfoTa                  matlab.ui.control.TextArea
        MsgState                logical
        Msg                     char
        MsgNum                  int16
        S1Tool2DBtn                         matlab.ui.control.Button
        S1Tool3DBtn                         matlab.ui.control.Button
        S1ToolExtract2DLineBtn              matlab.ui.control.Button
        S1ToolExtract3DLineBtn              matlab.ui.control.Button
        S1ToolExtractSurfBtn                matlab.ui.control.Button
        S2DesignSimulAsphericConcentricBtn  matlab.ui.control.Button
        S2DesignSimulFreeformBtn            matlab.ui.control.Button
        S4OptimAsphericConcentricBtn        matlab.ui.control.Button
    end

    % properties for the opening app
    properties (Access = private)
        AddSurfaceApp
    end
    
    % properties that should be used in the .m program
    properties (Access = public)
        workspaceDir                string
        toolPathName                string
        unit                        char
        fontName                    string
        fontSize            (1,1)   double
        toolFitType                 string
        arcFitMethod                string
        arcRansacMaxDist    (1,1)   double
        lineFitMethod       (1,1)   string
        lineFitMaxDist      (1,1)   double
        paramMethod                 string
        toolOri            
        surfType                    char
        surfDomain          (2,2)   double
        surfPathName                char
        surfPt              (3,:)   double
        surfFuncs                   function_handle
        surfFx                      function_handle
        surfFy                      function_handle
        cutDirection                string
        spindleDirection            string
        angularDiscrete             string
        aimRes              (1,1)   double
        rStep               (1,1)   double
        arcLength           (1,1)   double
        maxAngPtDist        (1,1)   double
        angularLength       (1,1)   double
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
                    app.SurfaceDetailTa.Value{5} = ['    X: ',num2str(app.surfDomain(1,:)),' (um)'];
                    app.SurfaceDetailTa.Value{6} = ['    Y: ',num2str(app.surfDomain(2,:)),' (um)'];
                case 'Function-Based'
                    app.SurfaceDetailTa.Value{2} = '- Surface function: ';
                    app.SurfaceDetailTa.Value{3} = char(surfString);
                    app.SurfaceDetailTa.Value{4} = '- Surface Domain: ';
                    app.SurfaceDetailTa.Value{5} = ['    X: ',num2str(app.surfDomain(1,:)),' (um)'];
                    app.SurfaceDetailTa.Value{6} = ['    Y: ',num2str(app.surfDomain(2,:)),' (um)'];
                otherwise
                    app.SurfaceDetailTa.Value{2} = '- Surface Function: ';
                    % app.SurfaceDetailTa.Value{3} = ...
                    %     '$Z = \sqrt{C^{2}(D-\frac{x^{2}}{A^{2}}-\frac{y^{2}}{B^{2}})}$';
                    app.SurfaceDetailTa.Value{3} = char(surfString);
                    app.SurfaceDetailTa.Value{4} = '- Surface Domain: ';
                    app.SurfaceDetailTa.Value{5} = ['    X: ',num2str(app.surfDomain(1,:)),' (um)'];
                    app.SurfaceDetailTa.Value{6} = ['    Y: ',num2str(app.surfDomain(2,:)),' (um)'];
                    app.surfFuncs = matlabFunction(surfString);
            end

            % Re-enable the AddSurface buttion
            app.AddSurfaceBtn.Enable = 'on';
        end
    end
    
    % Functions that is used in the callbacks
    methods (Access = private)
        resetToolfitParams(app);

        resetOptimParams(app);
    end

    % Callbacks that handle component events
    methods (Access = private)
        % Code that executes after component creation
        function startupFcn(app)
            app.MsgNum = 1;
            app.Msg = 'APP is successfully initialized!';
            InfoTaValueChanged(app,true);
%             app.WorkspaceDirEf.Value = 'workspace';
%             [app.MsgState,app.Msg] = mkdir('workspace');
%             if ~app.MsgState
%                 InfoTaValueChanged(app,true)
%             end
            resetToolfitParams(app);
            resetOptimParams(app);
            % updateSurface(app,app.surfType,app.surfString);
            % addpath(genpath('funcs'));
            app.Msg = 'All the parameters have been reset.';
            InfoTaValueChanged(app,true);
            app.Msg = 'Please choose a directory for the workspace.';
            InfoTaValueChanged(app,true);
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
        end

        % CLicked toolbar toggletool: set the window style
        function TopToggletoolClicked(app,event)
            switch app.UIFigure.WindowStyle
                case 'normal'
                    app.UIFigure.WindowStyle = 'alwaysontop';
                case 'alwaysontop'
                    app.UIFigure.WindowStyle = 'normal';
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
        function ToolFileBtnPushed(app,event)
            if isempty(app.workspaceDir)
                uialert(app.UIFigure, ...
                    'Workspace directory doesn''t exist, please choose one first!', ...
                    'Alert Message','CloseFcn',createCallbackFcn(app,@WorkspaceDirBtnPushed,true));
                return
            end
            [fileName,dirName] = uigetfile({ ...
                '*.csv','Comma-Separated Values-files(*.csv)'; ...
                '*.mat','MAT-files(*.mat)'; ...
                '*.txt','text-files(*.txt)'; ...
                '*.*','all files(*.*)'...
                }, ...
                'Select One Tool Tip Measurement Data', ...
                app.workspaceDir, ...
                'MultiSelect','off');
            app.toolPathName = fullfile(dirName,fileName);
            app.ToolFileEf.Value = app.toolPathName;
            app.Msg = 'Please select the parameters and click ''Update''.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: DdFontName i.e., font type selection
        function FontNameDdValueChanged(app,event)
            app.fontName = app.FontNameDd.Value;
            set(app.ToolDataAxes,'FontName',app.fontName);
        end
        
        % Value changed function: font size selection
        function FontSizeEfValueChanged(app,event)
            app.fontSize = app.FontSizeEf.Value;
            set(app.ToolDataAxes,'FontSize',app.fontSize);
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

        % Value changed function: update the selection selection
        function ToolUpdateBtnPushed(app,event)
            app.unit = app.UnitDd.Value;
            app.fontName = app.FontNameDd.Value;
            app.fontSize = app.FontSizeEf.Value;

            app.toolFitType = app.ToolFitTypeDd.Value;
            app.arcFitMethod = app.ArcFitMethodDd.Value;
            app.arcRansacMaxDist = app.ArcRansacMaxDistEf.Value;
            app.lineFitMethod = app.LineFitMethodDd.Value;
            app.lineFitMaxDist = app.LineFitMaxDistEf.Value;

            app.paramMethod = app.ParamMethodDd.Value;

            app.Msg = 'All the parameters are set. Click ''Plot'' to plot the data.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: reset the parameters
        function ToolResetBtnPushed(app,event)
            % Reset all the values to the default
            resetToolfitParams(app);

            % Report the infomation
            app.Msg = ['All the parameters in Tool File Processing are reset.', ... 
                'Modify them and click ''Update'' to set.'];
            InfoTaValueChanged(app,true);
        end

        % Value changed function: plot the tool file data
        function ToolPlotBtnPushed(app,event)
            % ensure the workspace directory and tool file path have been difined
            if isempty(app.workspaceDir)
                uialert(app.UIFigure,{'Invalid workspace directory:', ...
                    'workspace directory should not be empty'},'Alert Message');
                return;
            elseif isempty(app.toolPathName)
                uialert(app.UIFigure,{'Invalid tool file name:', ...
                    'tool file name should not be empty'},'Alert Message');
                return;
            end
            % get rid of the header of the csv file
            numHeader = 0;
            tooltipFile = fopen(app.toolPathName);
            while ~feof(tooltipFile)
                tmpLine = fgets(tooltipFile);
                % if the line begins with %d%d or -%d, then break
                if ~isnan(str2double(tmpLine(1:2)))
                    break;
                end
                numHeader = numHeader + 1;
            end
            fclose(tooltipFile);
            app.toolOri = importdata(app.toolPathName,',',numHeader);
            app.toolOri = app.toolOri.data;
            app.toolOri(:,3) = [];
            app.toolOri = app.toolOri';
            % plot the importing results
            plot(app.ToolDataAxes,app.toolOri(1,:),app.toolOri(2,:),'.','MarkerSize',2);
            hold(app.ToolDataAxes,'on');
            grid(app.ToolDataAxes,'on');
            xlabel(app.ToolDataAxes,['x (',app.unit,')']);
            ylabel(app.ToolDataAxes,['y (',app.unit,')']);
            app.Msg = 'Tool data is successfully loaded.';
            InfoTaValueChanged(app,true);
        end

        % Value changed function: transfer the parameters to the main program
        function Tool2DLineBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure,'Continue to fit the tool tip?',...
                'Confirmation');
            switch selection
                case 'OK'
                    s1_toolExtract_2Dline;
                    app.Msg = 'Tool tip is fitted successfully.';
                    InfoTaValueChanged(app,true);
                    app.CheckToolLamp.Color = 'g';
                case 'Cancel'
                    return
            end
        end

        % Value changed function: cancel the process with nothing to be saved
        function ToolCancelBtnPushed(app,event)
            resetToolfitParams(app);
            clf(app.ToolDataAxes,'reset');
            title(app.ToolDataAxes,'tool original data');
            app.CheckToolLamp.Color = 'g';
        end

        % Code that update the infomation of the EfInfo window
        function InfoTaValueChanged(app,event)
            app.MsgNum = app.MsgNum + 1;
            app.InfoTa.Value{app.MsgNum} = char([num2str(app.MsgNum),' ',app.Msg]);
            scroll(app.InfoTa,"bottom");
        end

        % Button down: shift to the surface tab, and refresh the message
        function SurfaceTbButtonDown(app,event)
            app.Msg = ['Switch to the surface tab. ', ...
                'Press the corresponding button to finish the programming process.'];
            InfoTaValueChanged(app,true);
        end

        % Button down: open the surface adding page
        function AddSurfaceBtnPushed(app,event)
            % Disable Plot Options button while dialog is open
            app.AddSurfaceBtn.Enable = 'off';

            % Open the add surface dialog
            app.surfType = '3D Geometry';
            app.AddSurfaceApp = add_surface(app,app.surfType,'No surface is activated now.');
        end

        % Value changed function: plot the tool file data
        function SurfacePlotBtnPushed(app,event)
            % ensure the workspace directory and tool file path have been difined
            % if isempty(app.workspaceDir)
            %     uialert(app.UIFigure,{'Invalid workspace directory:', ...
            %         'workspace directory should not be empty'},'Alert Message');
            %     return;
            % elseif isempty(app.surfPathName)
            %     uialert(app.UIFigure,{'Invalid surface file name:', ...
            %         'surface file name should not be empty'},'Alert Message');
            %     return;
            % end
            % get rid of the header of the csv file
            % numHeader = 0;
            % surfaceFile = fopen(app.surfPathName);
            % while ~feof(surfaceFile)
            %     tmpLine = fgets(surfaceFile);
            %     % if the line begins with %d%d or -%d, then break
            %     if ~isnan(str2double(tmpLine(1:2)))
            %         break;
            %     end
            %     numHeader = numHeader + 1;
            % end
            % fclose(tooltipFile);
            % app.surfOri = importdata(app.surfPathName,',',numHeader);

            % process the surface function or point cloud
            switch app.surfType
                case {'Point Cloud','Mesh'}
                    uialert(app.UIFigure,'There has not been methods fot import surface. ', ...
                        'Warning','Icon','warning');

                case {'Aspheric'} % 2D Geometry
                    % partial differential function
                    syms x;
                    app.surfFx = matlabFunction(diff(app.surfFuncs,x));

                    % sampling
                    spar = 501;
                    conR = linspace(0,R/4,spar); % concentric radius vector
                    conTheta = linspace(0,2*pi,spar);
                    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
                    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
                    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
                    surfMesh(:,:,3) = app.surfFuncs(surfMesh(:,:,1),surfMesh(:,:,2));

                    % plot the importing results
                    surf(app.SurfaceDataAxes, ...
                        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
                        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                    hold(app.SurfaceDataAxes,'on');
                    grid(app.SurfaceDataAxes,'on');
                    set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                        'FontName',app.fontName,'ZDir','reverse');
                    xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                    ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                    ylabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                    app.Msg = 'Surface data is successfully loaded.';
                    InfoTaValueChanged(app,true);

                case {'Ellipsoid','Function-Based'} % 3D Geometry
                    % partial differential function
                    syms x y;
                    app.surfFx = matlabFunction(diff(app.surfFuncs,x));
                    app.surfFy = matlabFunction(diff(app.surfFuncs,y));

                    % sampling
                    spar = 501;
                    xx = linspace(app.surfDomain(1,1),app.surfDomain(1,2),spar);
                    yy = linspace(app.surfDomain(2,1),app.surfDomain(2,2),spar);
                    [surfMesh(:,:,1),surfMesh(:,:,2)] = meshgrid(xx,yy);
                    surfMesh(:,:,3) = app.surfFuncs(surfMesh(:,:,1),surfMesh(:,:,2));

                    % plot the importing results
                    surf(app.SurfaceDataAxes, ...
                        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
                        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                    hold(app.SurfaceDataAxes,'on');
                    grid(app.SurfaceDataAxes,'on');
                    set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                        'FontName',app.fontName,'ZDir','reverse');
                    xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                    ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                    ylabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                    app.Msg = 'Surface data is successfully loaded.';
                    InfoTaValueChanged(app,true);
            end
        end

        % Value changed function: save the surface data
        function SurfaceSavedBtnPushed(app,event)
            selection = uiconfirm(app.UIFigure, ...
                'Save the surface data?','Comfirmation');
            switch selection
                case 'OK'
                    app.Msg = 'Surface is corrected.';
                    InfoTaValueChanged(app,true);
                    app.CheckSurfLamp.Color = 'g';
                case 'Cancel'
                    return
            end
        end

        % Value changed function: cancel the process with nothing to be saved
        function SurfaceCancelBtnPushed(app,event)
            resetSurfaceParams(app);
            clf(app.SurfaceDataAxes,'reset');
            title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            app.CheckSurfLamp.Color = 'g';
        end

        % Button down: shift to the program tab, and refresh the message
        function ProgramTbButtonDown(app,event)
            app.Msg = ['Switch to the program tab. ', ...
                'Press the corresponding button to finish the programming process.'];
            InfoTaValueChanged(app,true);
        end

        % Button down: shift to the optimization tab, and refresh the message
        function OptimTbButtonDown(app,event)
            % if
            %     app.CheckToolLamp.Color = 'g';
            % else
            %     app.CheckToolLamp.Color = 'r';
            % end
            % if
            %     app.CheckSurfLamp.Color = 'g';
            % else
            %     app.CheckSurfLamp.Color = 'r';
            % end
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
            app.spindleDirection = app.SpindleDirectionDd.Value;
            app.angularDiscrete = app.AngularDiscreteDd.Value;
            app.aimRes = app.AimResEf.Value;
            app.rStep = app.RStepEf.Value;
            app.arcLength = app.ArcLengthEf.Value;
            app.maxAngPtDist = app.MaxAngPtDistEf.Value;
            app.angularLength = app.AngularLengthEf.Value;

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
        
        function S4OptimAsphericConcentricBtnPushed(app,event)
            s4_optim_aspheric_concentric;
        end

        % --------------------------Function Execution--------------------------
        function S1Tool2DBtnPushed(app,event)
            s1_tool2D
        end

        function S1Tool3DBtnPushed(app,event)
            s1_tool3D
        end

        function S1ToolExtract2DLineBtnPushed(app,event)
            s1_toolExtract_2Dline
        end

        function S1ToolExtract3DLineBtnPushed(app,event)
            s1_toolExtract_3Dline
        end

        function S1ToolExtractSurfBtnPushed(app,event)
            s1_toolExtract_surf
        end

        function S2DesignSimulAsphericConcentricBtnPushed(app,event)
            s2_design_simul_aspheric_concentric
        end

        function S2DesignSimulFreeformBtnPushed(app,event)
            s2_design_simul_freeform
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
                'WindowStyle','alwaysontop','WindowState','normal','Visible','off');
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app,@UIFigureCloseReq,true);
            app.UIFigure.Resize = "on";
            app.UIFigure.Position = [500,500,800,600];
            app.UIFigure.Scrollable = "on";

            % Manage menu
            app.FigureMenu = uimenu(app.UIFigure,'Text','Options');

            % Manage toolbar
            app.FigureToolbar = uitoolbar(app.UIFigure);

            app.CloseAllFigurePushtool = uipushtool(app.FigureToolbar, ...
                'Icon','resource/image/CloseAllFigure.svg');
            app.CloseAllFigurePushtool.ClickedCallback = createCallbackFcn( ...
                app,@CloseAllFigurePushtoolClicked,true);
            app.CloseAllFigurePushtool.Tooltip = 'Close all the figures';

            app.BoldInfoToggletool = uitoggletool(app.FigureToolbar, ...
                'Icon','resource/image/Bold.svg','Separator','on');
            app.BoldInfoToggletool.ClickedCallback = createCallbackFcn( ...
                app,@BoldInfoToggletoolClicked,true);
            app.BoldInfoToggletool.Tooltip = 'Bold/Normalize the infomation.';

            app.TopToggletool = uitoggletool(app.FigureToolbar, ...
                'Icon','resource/image/top.svg');
            app.TopToggletool.ClickedCallback = createCallbackFcn( ...
                app,@TopToggletoolClicked,true);
            app.CloseAllFigurePushtool.Tooltip = 'Put the APP on top';

            % Manage tab groups
            FigureGl = uigridlayout(app.UIFigure,[3,1],'Padding',[5,5,5,5]);
            FigureGl.RowHeight = {'fit','1x','fit'};
            FigureGl.ColumnWidth = {'1x'};

            % ------------------------Workspace directory------------------------
            WorkspaceDirGl = uigridlayout(FigureGl,[1,3]);
            WorkspaceDirGl.Layout.Row = 1;
            WorkspaceDirGl.RowHeight = {'fit'};
            WorkspaceDirGl.ColumnWidth = {'fit','1x','fit'};
            WorkspaceDirGl.Padding = [0,0,0,0];
            
            WorkspaceDirLb = uilabel(WorkspaceDirGl,'Text','Workspace directory:');
            WorkspaceDirLb.Layout.Row = 1;
            WorkspaceDirLb.Layout.Column = 1;
            
            app.WorkspaceDirEf = uieditfield(WorkspaceDirGl,'text');
            app.WorkspaceDirEf.Layout.Row = 1;
            app.WorkspaceDirEf.Layout.Column = 2;
            app.WorkspaceDirEf.ValueChangedFcn = createCallbackFcn(app,@WorkspaceDirEfValueChanged,true);
            
            app.WorkspaceDirBtn = uibutton(WorkspaceDirGl,'push','Text','Choose');
            app.WorkspaceDirBtn.Layout.Row = 1;
            app.WorkspaceDirBtn.Layout.Column = 3;
            app.WorkspaceDirBtn.ButtonPushedFcn = createCallbackFcn(app,@WorkspaceDirBtnPushed,true);

            % figure tab
            app.FigureTbGp = uitabgroup(FigureGl,'SelectedTab',app.ToolTb);
            app.FigureTbGp.Layout.Row = 2;
            app.FigureTbGp.Layout.Column = 1;
            
            % ------------------------------------------------------------------------
            % --------------------------Tool File Processing--------------------------
            % ------------------------------------------------------------------------

            app.ToolTb = uitab(app.FigureTbGp,'Title','Tool File Processing');
            app.ToolTb.ButtonDownFcn = createCallbackFcn(app,@ToolTbButtonDown,true);

            % Manage tool processing layout
            ToolTbGl = uigridlayout(app.ToolTb,[4,4]);
            ToolTbGl.RowHeight = {'fit','fit','fit','1x','fit','fit'};
            ToolTbGl.ColumnWidth = {'fit','1x','1x','1x'};
            % the elements of the cell array can be 'fit', fixed pixels, or '1x' '2x'
            %   'fit': to adjust the size to show the whole text, or basedon the default size
            %   fixed pixels: the size is fixed at the number of pixels
            %   variables ('2x'): to fill the remaining space, and the number is a
            %   weight for dividing up the remaining space among all variables
                 
            % ------------------------Tool data importing layout------------------------
            ToolFileGl = uigridlayout(ToolTbGl,[1,3]);
            ToolFileGl.Layout.Row = 1;
            ToolFileGl.Layout.Column = [1,4];
            ToolFileGl.ColumnWidth = {'fit','1x','fit'};
            ToolFileGl.Padding = [0,0,0,0];
            
            ToolFileLb = uilabel(ToolFileGl,'Text','Tool data path:');
            ToolFileLb.Layout.Row = 1;
            ToolFileLb.Layout.Column = 1;
            
            app.ToolFileEf = uieditfield(ToolFileGl,'text');
            app.ToolFileEf.Layout.Row = 1;
            app.ToolFileEf.Layout.Column = 2;
            app.ToolFileEf.ValueChangedFcn = createCallbackFcn(app,@ToolFileEfValueChanged,true);
            
            app.ToolFileBtn = uibutton(ToolFileGl,'push','Text','Choose');
            app.ToolFileBtn.Layout.Row = 1;
            app.ToolFileBtn.Layout.Column = 3;
            app.ToolFileBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolFileBtnPushed,true);

            % ------------------------Parameter editing tabgroup------------------------
            % Create the parameter selection tab group
            ParamPn = uipanel(ToolTbGl,'Title','Tool process', ...
                'TitlePosition','centertop','Scrollable','on');
            ParamPn.Layout.Row = [2,3];
            ParamPn.Layout.Column = 1;
            ParamPn.Scrollable = 'on';
            ParamPnGl = uigridlayout(ParamPn,[2,2]);
            ParamPnGl.ColumnWidth = {'1x','1x'};
            ParamPnGl.RowHeight = {'1x','fit'};

            ParamTabgroup = uitabgroup(ParamPnGl,'SelectedTab',app.ParamTabCommon);
            ParamTabgroup.Layout.Row = 1;
            ParamTabgroup.Layout.Column = [1,2];

            % ---Common parameter selection---
            app.ParamTabCommon = uitab(ParamTabgroup,'Title','Common','Scrollable','on');
            ParamTabCommonGl = uigridlayout(app.ParamTabCommon,[3,2]);
            ParamTabCommonGl.RowHeight = {'fit','fit','fit'};
            ParamTabCommonGl.ColumnWidth = {'fit','1x'};
            ParamTabCommonGl.Scrollable = 'on';
            
            % Create the unit selection label and dropdown box
            UnitLB = uilabel(ParamTabCommonGl,'Text','Data unit');
            UnitLB.Layout.Row = 1;
            UnitLB.Layout.Column = 1;
            app.UnitDd = uidropdown(ParamTabCommonGl,'Items',{'m','mm','/mum','nm'},'Editable','on', ...
                'Value',app.unitDefault,'BackgroundColor',[1,1,1]);
            app.UnitDd.Layout.Row = 1;
            app.UnitDd.Layout.Column = 2;
            
            % Create the font name selection label and dropdown box
            FontNameLb = uilabel(ParamTabCommonGl,'Text','Ploting font type');
            FontNameLb.Layout.Row = 2;
            FontNameLb.Layout.Column = 1;
            app.FontNameDd = uidropdown(ParamTabCommonGl,'Items',listfonts,'Value',app.fontNameDefault, ...
                'BackgroundColor',[1,1,1],'Visible','on');
            app.FontNameDd.Layout.Row = 2;
            app.FontNameDd.Layout.Column = 2;
            app.FontNameDd.ValueChangedFcn = createCallbackFcn(app,@FontNameDdValueChanged,true);
            
            % Create the font size selection label and editfield box
            FontSizeLb = uilabel(ParamTabCommonGl,'Text','Ploting font size');
            FontSizeLb.Layout.Row = 3;
            FontSizeLb.Layout.Column = 1;
            app.FontSizeEf = uieditfield(ParamTabCommonGl,'numeric','Limits',[6 72], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.fontSizeDefault);
            app.FontSizeEf.Layout.Row = 3;
            app.FontSizeEf.Layout.Column = 2;
            app.FontSizeEf.ValueChangedFcn = createCallbackFcn(app,@FontSizeEfValueChanged,true);

            % ---Tool fitting parameter selection---
            app.ParamTabToolfit = uitab(ParamTabgroup,'Title','Tool fit','Scrollable','on');
            ParamTabToolfitGl = uigridlayout(app.ParamTabToolfit,[5,2]);
            ParamTabToolfitGl.RowHeight = {'fit','fit','fit','fit','fit'};
            ParamTabToolfitGl.ColumnWidth = {'fit','1x'};
            ParamTabToolfitGl.Scrollable = 'on';

            % Create the tool fitting type label and dropdown box
            ToolFitTypeLb = uilabel(ParamTabToolfitGl,'Text','Tool fitting type');
            ToolFitTypeLb.Layout.Row = 1;
            ToolFitTypeLb.Layout.Column = 1;
            app.ToolFitTypeDd = uidropdown(ParamTabToolfitGl,'Items',{'onlyArc','arcRansac','lineArc'}, ...
                'Value',app.toolFitTypeDefault,'BackgroundColor',[1,1,1]);
            app.ToolFitTypeDd.Layout.Row = 1;
            app.ToolFitTypeDd.Layout.Column = 2;
            app.ToolFitTypeDd.ValueChangedFcn = createCallbackFcn(app,@ToolFitTypeDdValueChanged,true);
            app.ToolFitTypeDd.DropDownOpeningFcn = createCallbackFcn(app,@ToolFitTypeDdOpening,true);
            
            % Create the arc fitting method label and dropdown box
            ArcFitMethodLb = uilabel(ParamTabToolfitGl,'Text','Arc fitting method');
            ArcFitMethodLb.Layout.Row = 2;
            ArcFitMethodLb.Layout.Column = 1;
            app.ArcFitMethodDd = uidropdown(ParamTabToolfitGl, ...
                'Items',{'gradient-decent','normal-equation','levenberg-marquardt'}, ...
                'Value',app.arcFitMethodDefault);
            app.ArcFitMethodDd.Layout.Row = 2;
            app.ArcFitMethodDd.Layout.Column = 2;
            app.ArcFitMethodDd.BackgroundColor = [1 1 1];

            % Create the arc ransac fitting max distance label and dropdown box
            ArcRansacMaxDistLb = uilabel(ParamTabToolfitGl,'Text',{'Arc ransac fitting','max distance'});
            ArcRansacMaxDistLb.Layout.Row = 3;
            ArcRansacMaxDistLb.Layout.Column = 1;
            app.ArcRansacMaxDistEf = uieditfield(ParamTabToolfitGl,'numeric','Limits',[0 inf], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.arcRansacMaxDistDefault);
            app.ArcRansacMaxDistEf.Layout.Row = 3;
            app.ArcRansacMaxDistEf.Layout.Column = 2;
            app.ArcRansacMaxDistEf.Enable = 'off';
            app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
            % app.ArcRansacMaxDistEf.ValueChangedFcn = createCallbackFcn(app,@FontSizeEfValueChanged,true);

            % Create the line fitting method label and dropdown box
            LineFitMethodLb = uilabel(ParamTabToolfitGl,'Text','Line fitting method');
            LineFitMethodLb.Layout.Row = 4;
            LineFitMethodLb.Layout.Column = 1;
            app.LineFitMethodDd = uidropdown(ParamTabToolfitGl, ...
                'Items',{'polyfit','ransac'}, ...
                'Value',app.lineFitMethodDefault);
            app.LineFitMethodDd.Layout.Row = 4;
            app.LineFitMethodDd.Layout.Column = 2;
            app.LineFitMethodDd.Enable = "off";
            app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];

            % Create the line fitting max distance label and editfield box
            LineFitMaxDistEfLb = uilabel(ParamTabToolfitGl,'Text',{'Line fitting', 'max distance'});
            LineFitMaxDistEfLb.Layout.Row = 5;
            LineFitMaxDistEfLb.Layout.Column = 1;
            app.LineFitMaxDistEf = uieditfield(ParamTabToolfitGl,'numeric','Limits',[0 inf], ...
                'HorizontalAlignment','center','ValueDisplayFormat','%d','Value',app.lineFitMaxDistDefault);
            app.LineFitMaxDistEf.Layout.Row = 5;
            app.LineFitMaxDistEf.Layout.Column = 2;
            app.LineFitMaxDistEf.Enable = "off";
            app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];

            % ---Tool Interpolation parameter selection---
            app.ParamTabToolinterp = uitab(ParamTabgroup,'Title','Tool Interpolation','Scrollable','on');
            ParamTabToolinterpGl = uigridlayout(app.ParamTabToolinterp,[3,2]);
            ParamTabToolinterpGl.RowHeight = {'fit','fit','fit'};
            ParamTabToolinterpGl.ColumnWidth = {'fit','1x'};
            ParamTabToolinterpGl.Scrollable = 'on';

            % Create the B-spline parametric method label and dropdown box
            ParamMethodLb = uilabel(ParamTabToolinterpGl,'Text',{'B-spline', 'param-method'});
            ParamMethodLb.Layout.Row = 2;
            ParamMethodLb.Layout.Column = 1;
            app.ParamMethodDd = uidropdown(ParamTabToolinterpGl, ...
                'Items',{'uniform','chord','centripetal'}, ...
                'Value',app.paramMethodDefault,'BackgroundColor',[1,1,1]);
            app.ParamMethodDd.Layout.Row = 2;
            app.ParamMethodDd.Layout.Column = 2;
            
            % ---Create the reset button---
            app.ToolResetBtn = uibutton(ParamPnGl,'push','Text','Reset','Visible','on');
            app.ToolResetBtn.Layout.Row = 2;
            app.ToolResetBtn.Layout.Column = 1;
            app.ToolResetBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolResetBtnPushed,true);
            
            % ---Create the update button---
            app.ToolUpdateBtn = uibutton(ParamPnGl,'push','Text','Update','Visible','on');
            app.ToolUpdateBtn.Layout.Row = 2;
            app.ToolUpdateBtn.Layout.Column = 2;
            app.ToolUpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolUpdateBtnPushed,true);

            % ------------------------Ploting axes------------------------
            app.ToolDataAxes = uiaxes(ToolTbGl);
            title(app.ToolDataAxes,'tool original data','FontSize',16);
            % xlabel(app.UIAxes, 'X');
            % ylabel(app.UIAxes, 'Y');
            % zlabel(app.UIAxes, 'Z');
            app.ToolDataAxes.Layout.Row = 2;
            app.ToolDataAxes.Layout.Column = [2,4];

            app.ToolPlotBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Extract & Plot');
            app.ToolPlotBtn.Layout.Row = 3;
            app.ToolPlotBtn.Layout.Column = 2;
            app.ToolPlotBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolPlotBtnPushed,true);

            % ------------------------Ending process------------------------
            app.Tool2DLineBtnBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Enter');
            app.Tool2DLineBtnBtn.Layout.Row = 3;
            app.Tool2DLineBtnBtn.Layout.Column = 3;
            app.Tool2DLineBtnBtn.ButtonPushedFcn = createCallbackFcn(app,@Tool2DLineBtnPushed,true);
            
            app.Tool3DLineBtnBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Enter');
            app.Tool3DLineBtnBtn.Layout.Row = 3;
            app.Tool3DLineBtnBtn.Layout.Column = 3;
            app.Tool3DLineBtnBtn.Visible = 'off';
            app.Tool3DLineBtnBtn.ButtonPushedFcn = createCallbackFcn(app,@Tool2DLineBtnPushed,true);

            app.ToolCancelBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Cancel');
            app.ToolCancelBtn.Layout.Row = 3;
            app.ToolCancelBtn.Layout.Column = 4;
            app.ToolCancelBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolCancelBtnPushed,true);

            % ------------------------------------------------------------------------
            % ---------------------------Surface Processing---------------------------
            % ------------------------------------------------------------------------

            app.SurfaceTb = uitab(app.FigureTbGp,'Title','Surface Load');
            app.SurfaceTb.ButtonDownFcn = createCallbackFcn(app,@SurfaceTbButtonDown,true);

            SurfaceTbGl = uigridlayout(app.SurfaceTb,[3,4]);
            SurfaceTbGl.RowHeight = {'fit','1x','fit'};
            SurfaceTbGl.ColumnWidth = {'2x','1x','1x','1x'};

            app.AddSurfaceBtn = uibutton(SurfaceTbGl,'push','WordWrap','on', ...
                'Text','Add Geometry');
            app.AddSurfaceBtn.Layout.Row = 1;
            app.AddSurfaceBtn.Layout.Column = 1;
            app.AddSurfaceBtn.Icon = 'resource/image/AddSurf1.svg';
            app.AddSurfaceBtn.ButtonPushedFcn = createCallbackFcn(app,@AddSurfaceBtnPushed,true);

            app.SurfaceDetailTa = uitextarea(SurfaceTbGl,'WordWrap','on', ...
                'Editable','off');
            app.SurfaceDetailTa.Layout.Row = [2,3];
            app.SurfaceDetailTa.Layout.Column = 1;
            app.SurfaceDetailTa.FontName = 'Times New Roman';
            app.SurfaceDetailTa.FontSize = 16;
            app.SurfaceDetailTa.FontWeight = 'normal';
            
            app.SurfaceDataAxes = uiaxes(SurfaceTbGl);
            title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            app.SurfaceDataAxes.Layout.Row = [1,2];
            app.SurfaceDataAxes.Layout.Column = [2,4];

            app.SurfacePlotBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','Generate & Plot');
            app.SurfacePlotBtn.Layout.Row = 3;
            app.SurfacePlotBtn.Layout.Column = 2;
            app.SurfacePlotBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfacePlotBtnPushed,true);

            % ------------------------Ending process------------------------
            app.SurfaceSavedBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','Enter');
            app.SurfaceSavedBtn.Layout.Row = 3;
            app.SurfaceSavedBtn.Layout.Column = 3;
            app.SurfaceSavedBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfaceSavedBtnPushed,true);
            
            app.SurfaceCancelBtn = uibutton(SurfaceTbGl,'push','WordWrap','on','Text','Cancel');
            app.SurfaceCancelBtn.Layout.Row = 3;
            app.SurfaceCancelBtn.Layout.Column = 4;
            app.SurfaceCancelBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfaceCancelBtnPushed,true);


            % ------------------------------------------------------------------------
            % ---------------------------Program Processing---------------------------
            % ------------------------------------------------------------------------

            app.ProgramTb = uitab(app.FigureTbGp,'Title','Program');
            app.ProgramTb.ButtonDownFcn = createCallbackFcn(app,@ProgramTbButtonDown,true);

            ProgramGl = uigridlayout(app.ProgramTb,[4,2]);
            ProgramGl.RowHeight = {'1x','1x','1x','fit'};
            ProgramGl.ColumnWidth = {'1x','1x'};

            % ---------------------------tool fitting process---------------------------
            ToolFitPn = uipanel(ProgramGl,'Visible','on', ...
                'Title','','TitlePosition','centertop');
            ToolFitPn.Layout.Row = [1,3];
            ToolFitPn.Layout.Column = 1;
            ToolFitPn.Scrollable = 'on';

            ToolFitGl = uigridlayout(ToolFitPn,[5,1]);
            ToolFitGl.RowHeight = {'1x','1x','1x','1x','1x'};
            ToolFitGl.ColumnWidth = {'1x'};

            % button to execute s1_tool2D.m
            app.S1Tool2DBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_tool2D');
            app.S1Tool2DBtn.Layout.Row = 1;
            app.S1Tool2DBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1Tool2DBtnPushed,true);

            % button to execute s1_tool3D.m
            app.S1Tool3DBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_tool3D');
            app.S1Tool3DBtn.Layout.Row = 2;
            app.S1Tool3DBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1Tool3DBtnPushed,true);

            % button to execute s1_toolExtract_2Dline.m
            app.S1ToolExtract2DLineBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_2DLine');
            app.S1ToolExtract2DLineBtn.Layout.Row = 3;
            app.S1ToolExtract2DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract2DLineBtnPushed,true);

            % button to execute s1_toolExtract_3Dline.m
            app.S1ToolExtract3DLineBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_3Dline');
            app.S1ToolExtract3DLineBtn.Layout.Row = 4;
            app.S1ToolExtract3DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract3DLineBtnPushed,true);

            % button to execute s1_toolExtract_surf.m
            app.S1ToolExtractSurfBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_surf');
            app.S1ToolExtractSurfBtn.Layout.Row = 5;
            app.S1ToolExtractSurfBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtractSurfBtnPushed,true);

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

            OptimArrow1 = uiimage(OptimGl,'ImageSource','resource/image/RightArrow6.svg', ...
                'ScaleMethod','stretch');
            OptimArrow1.Layout.Row = 1;
            OptimArrow1.Layout.Column = 2;

            OptimArrow2 = uiimage(OptimGl,'ImageSource','resource/image/RightArrow6.svg', ...
                'ScaleMethod','stretch');
            OptimArrow2.Layout.Row = 1;
            OptimArrow2.Layout.Column = 4;

            % ------------------------------------------------------------------------
            % --------------------------Optimization Process--------------------------
            % ------------------------------------------------------------------------

            app.OptimTb = uitab(app.FigureTbGp,'Title','Optimization');
            app.OptimTb.ButtonDownFcn = createCallbackFcn(app,@OptimTbButtonDown,true);
            % app.OptimTb.Scrollable = 'on';

            OptimTbGl = uigridlayout(app.OptimTb,[2,4]);
            OptimTbGl.RowHeight = {'fit','fit'};
            OptimTbGl.ColumnWidth = {'fit','1x','fit','1x'};

            CheckToolPn = uipanel(OptimTbGl);
            CheckToolPn.Layout.Row = 1;
            CheckToolPn.Layout.Column = [1,2];
            CheckToolPnGl = uigridlayout(CheckToolPn,[1,2]);
            CheckToolPnGl.RowHeight = {'fit'};
            CheckToolPnGl.ColumnWidth = {'fit','1x'};
            CheckToolLb = uilabel(CheckToolPnGl,'Text','Check Tool');
            CheckToolLb.Layout.Row = 1;
            CheckToolLb.Layout.Column = 1;
            app.CheckToolLamp = uilamp(CheckToolPnGl,'Color','r');
            app.CheckToolLamp.Layout.Row = 1;
            app.CheckToolLamp.Layout.Column = 2;

            CheckSurfPn = uipanel(OptimTbGl);
            CheckSurfPn.Layout.Row = 1;
            CheckSurfPn.Layout.Column = [3,4];
            CheckSurfPnGl = uigridlayout(CheckSurfPn,[1,2]);
            CheckSurfPnGl.RowHeight = {'fit'};
            CheckSurfPnGl.ColumnWidth = {'fit','1x'};
            CheckSurfLb = uilabel(CheckSurfPnGl,'Text','Check Surface');
            CheckSurfLb.Layout.Row = 1;
            CheckSurfLb.Layout.Column = 1;
            app.CheckSurfLamp = uilamp(CheckSurfPnGl,'Color','r');
            app.CheckSurfLamp.Layout.Row = 1;
            app.CheckSurfLamp.Layout.Column = 2;

            OptimParamPn = uipanel(OptimTbGl,'Title','Edit Parameters', ...
                'TitlePosition','centertop');
            OptimParamPn.Layout.Row = 2;
            OptimParamPn.Layout.Column = [1,3];
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
            OptimParamPathTabGl = uigridlayout(app.OptimParamPathTab,[7,3]);
            OptimParamPathTabGl.Scrollable = 'on';
            OptimParamPathTabGl.RowHeight = {'fit','fit','fit','fit','fit','fit','fit'};
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
            app.SpindleDirectionDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'Counterclockwise','Clockwise'}, ...
                'Value',app.spindleDirectionDefault);
            app.SpindleDirectionDd.Layout.Row = 2;
            app.SpindleDirectionDd.Layout.Column = 2;

            AngularDiscreteLb = uilabel(OptimParamPathTabGl,'Text','Spindle Direction');
            AngularDiscreteLb.Layout.Row = 3;
            AngularDiscreteLb.Layout.Column = 1;
            app.AngularDiscreteDd = uidropdown(OptimParamPathTabGl, ...
                'Items',{'Constant Arc','Constant Angle'}, ...
                'Value',app.angularDiscreteDefault);
            app.AngularDiscreteDd.Layout.Row = 3;
            app.AngularDiscreteDd.Layout.Column = 2;
            app.AngularDiscreteDd.ValueChangedFcn = createCallbackFcn(app,@AngularDiscreteDdValueChanged,true);

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
            AimResUnitLb.Text = '$\mu$m';

            RStepLb = uilabel(OptimParamPathTabGl,'Text','Initial Annulus Width');
            RStepLb.Layout.Row = 5;
            RStepLb.Layout.Column = 1;
            app.RStepEf = uieditfield(OptimParamPathTabGl,'numeric','Limits',[0,inf], ...
                'HorizontalAlignment','center','Value',app.rStepDefault);
            app.RStepEf.Layout.Row = 5;
            app.RStepEf.Layout.Column = 2;
            RStepUnitLb = uilabel(OptimParamPathTabGl,'Interpreter','latex');
            RStepUnitLb.Layout.Row = 5;
            RStepUnitLb.Layout.Column = 3;
            RStepUnitLb.Text = '$\mu$m';

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
            app.ArcLengthUnitLb.Text = '$\mu$m';

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

            % feed tab
            app.OptimParamFeedTab = uitab(OptimParamTabgroup,'Title','Feed');
            app.OptimParamFeedTab.Scrollable = 'on';

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

            app.S4OptimAsphericConcentricBtn = uibutton(OptimTbGl,'push','WordWrap','on', ...
                'Text','s4_optim_aspheric_concentric');
            app.S4OptimAsphericConcentricBtn.Layout.Row = 2;
            app.S4OptimAsphericConcentricBtn.Layout.Column = 4;
            app.S4OptimAsphericConcentricBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S4OptimAsphericConcentricBtnPushed,true);

            % ------------------------Info displaying window------------------------
            app.InfoTa = uitextarea(FigureGl,'WordWrap','on', ...
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
        function app = upm_toolpath
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
