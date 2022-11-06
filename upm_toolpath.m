classdef upm_toolpath < matlab.apps.AppBase
    % to create a ui figure for tool data & parameters input

    % Default properties
    properties (Constant)
        workspaceDirDefault = 'D:\Code\2021-11_ToolWaviness\upm_toolpath_waviness\workspace';
        toolFitTypeDefault = 'lineArc'
        arcFitMethodDefault = 'levenberg-marquardt'
        arcRansacMaxDistDefault = 1e-2
        lineFitMethodDefault = 'polyfit'
        lineFitMaxDistDefault = 1e-3
        paramMethodDefault = 'centripetal'
        unitDefault        = 'mm'
        fontNameDefault    = 'Times New Roman'
        fontSizeDefault    = 12
    end

    % Properties that correspond to app components
    properties (Access = public)
        PlotUIFigure            matlab.ui.Figure
        FigureMn                matlab.ui.container.Menu
        FigureToolbar           matlab.ui.container.Toolbar
        CloseAllFigurePushtool  matlab.ui.container.toolbar.PushTool
        BoldInfoToggletool      matlab.ui.container.toolbar.ToggleTool
        FigureTbGp              matlab.ui.container.TabGroup
        ToolTb                  matlab.ui.container.Tab
        SurfaceTb               matlab.ui.container.Tab
        ProgramTb               matlab.ui.container.Tab
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
        ResetBtn                matlab.ui.control.Button
        UpdateBtn               matlab.ui.control.Button
        ToolDataAxes            matlab.ui.control.UIAxes
        PlotBtn                 matlab.ui.control.Button
        EnterBtn                matlab.ui.control.Button
        CancelBtn               matlab.ui.control.Button
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
    end

    % properties that should be used in the .m program
    properties (Access = public)
        workspaceDir        string
        toolPathName        string
        unit                char
        fontName            string
        fontSize            double
        toolFitType         string
        arcFitMethod        string
        arcRansacMaxDist    double
        lineFitMethod       string
        lineFitMaxDist      double
        paramMethod         string
        toolOri
    end
    
    % Functions that is used in the callbacks
    methods (Access = private)
        function resetToolfitParams(app)
            app.UnitDd.Value = app.unitDefault;
            app.unit = app.UnitDd.Value;
            app.FontNameDd.Value = app.fontNameDefault;
            app.fontName = app.FontNameDd.Value;
            app.FontSizeEf.Value = app.fontSizeDefault;
            app.fontSize = app.FontSizeEf.Value;

            app.ToolFitTypeDd.Value = app.toolFitTypeDefault;
            app.toolFitType = app.ToolFitTypeDd.Value;
            app.ArcFitMethodDd.Value = app.arcFitMethodDefault;
            app.arcFitMethod = app.ArcFitMethodDd.Value;
            app.ArcRansacMaxDistEf.Value = app.arcRansacMaxDistDefault;
            app.arcRansacMaxDist = app.ArcRansacMaxDistEf.Value;
            app.ArcRansacMaxDistEf.Enable = "off";
            app.ArcRansacMaxDistEf.BackgroundColor = [0.96 0.96 0.96];
            app.LineFitMethodDd.Value = app.lineFitMethodDefault;
            app.lineFitMethod = app.LineFitMethodDd.Value;
            app.LineFitMethodDd.Enable = "off";
            app.LineFitMethodDd.BackgroundColor = [0.96 0.96 0.96];
            app.LineFitMaxDistEf.Value = app.lineFitMaxDistDefault;
            app.lineFitMaxDist = app.LineFitMaxDistEf.Value;
            app.LineFitMaxDistEf.Enable = "off";
            app.LineFitMaxDistEf.BackgroundColor = [0.96 0.96 0.96];

            app.ParamMethodDd.Value = app.paramMethodDefault;
            app.paramMethod = app.ParamMethodDd.Value;
        end

        function resetSurfaceParams(app)
        end

        function resetOptimParams(app)
        end
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
            resetSurfaceParams(app);
            resetOptimParams(app);
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

        % Button down: shift to the toolbar tab, and refresh the message
        function ToolTbButtonDown(app,event)
            app.Msg = 'Please choose a directory for the workspace.';
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
        function WorkspaceDirBtnValueChanged(app,event)
            app.workspaceDir = uigetdir(app.workspaceDirDefault, ...
                'Select the Workspace Directory');
            if isempty(app.workspaceDir)
                uialert(app.PlotUIFigure,{'Invalid workspace directory:', ...
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
                uialert(app.PlotUIFigure, ...
                    'Workspace directory doesn''t exist, please choose one first!', ...
                    'Alert Message','CloseFcn',createCallbackFcn(app,@WorkspaceDirBtnValueChanged,true));
                return
            end
            app.toolPathName = app.ToolFileEf.Value;
            app.Msg = 'Please select the parameters and click Update.';
            InfoTaValueChanged(app,true);
        end
        
        % Value changed function: to choose the tool file path
        function ToolFileBtnValueChanged(app,event)
            if isempty(app.workspaceDir)
                uialert(app.PlotUIFigure, ...
                    'Workspace directory doesn''t exist, please choose one first!', ...
                    'Alert Message','CloseFcn',createCallbackFcn(app,@WorkspaceDirBtnValueChanged,true));
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
        function UpdateBtnValueChanged(app,event)
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
        function ResetBtnValueChanged(app,event)
            % Reset all the values to the default
            resetToolfitParams(app);

            % Report the infomation
            app.Msg = ['All the parameters in Tool File Processing are reset.', ... 
                'Modify them and click ''Update'' to set.'];
            InfoTaValueChanged(app,true);
        end

        % Value changed function: plot the tool file data
        function PlotBtnValueChanged(app,event)
            % ensure the workspace directory and tool file path have been difined
            if isempty(app.workspaceDir)
                uialert(app.PlotUIFigure,{'Invalid workspace directory:', ...
                    'workspace directory should not be empty'},'Alert Message');
                return;
            elseif isempty(app.toolPathName)
                uialert(app.PlotUIFigure,{'Invalid tool file name:', ...
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
        function EnterBtnValueChanged(app,event)
            selection = uiconfirm(app.PlotUIFigure, ...
                'Save the original tool data and close the window?','Comfirmation');
            switch selection
                case 'OK'
                    delete(app.PlotUIFigure);
                case 'Cancel'
                    return
            end
        end

        % Value changed function: cancel the process with nothing to be saved
        function CancelBtnValueChanged(app,event)
            UIFigureCloseReq(app,true);
        end

        % Code that update the infomation of the EfInfo window
        function InfoTaValueChanged(app,event)
            app.MsgNum = app.MsgNum + 1;
            app.InfoTa.Value{app.MsgNum} = char([num2str(app.MsgNum),' ',app.Msg]);
            scroll(app.InfoTa,"bottom");
        end

        % --------------------------Function Execution--------------------------
        function S1Tool2DBtnValueChanged(app,event)
            s1_tool2D
        end

        function S1Tool3DBtnValueChanged(app,event)
            s1_tool3D
        end

        function S1ToolExtract2DLineBtnValueChanged(app,event)
            s1_toolExtract_2Dline
        end

        function S1ToolExtract3DLineBtnValueChanged(app,event)
            s1_toolExtract_3Dline
        end

        function S1ToolExtractSurfBtnValueChanged(app,event)
            s1_toolExtract_surf
        end

        % Button down: shift to the surface tab, and refresh the message
        function SurfaceTbButtonDown(app,event)
            app.Msg = '';
            InfoTaValueChanged(app,true);
        end

        % Button down: shift to the program tab, and refresh the message
        function ProgramTbButtonDown(app,event)
            app.Msg = 'Press the corresponding button to finish the programming process.';
            InfoTaValueChanged(app,true);
        end

        function S2DesignSimulAsphericConcentricBtnValueChanged(app,event)
            s2_design_simul_aspheric_concentric
        end

        function S2DesignSimulFreeformBtnValueChanged(app,event)
            s2_design_simul_freeform
        end

        % Cross pushed function: execute when pushing the cross of the UIFigure
        function UIFigureCloseReq(app,event)
            selection = uiconfirm(app.PlotUIFigure,'Close the figure window?',...
                'Confirmation');
            switch selection
                case 'OK'
                    delete(app.PlotUIFigure)
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
            app.PlotUIFigure = uifigure('Name','Tool Data & Parameters Input', ...
                'WindowStyle','alwaysontop','WindowState','normal','Visible','off');
            app.PlotUIFigure.CloseRequestFcn = createCallbackFcn(app,@UIFigureCloseReq,true);
            app.PlotUIFigure.Resize = "on";
            app.PlotUIFigure.Position = [1000,1000,1000,600];
            app.PlotUIFigure.Scrollable = "on";

            % Manage menu
            app.FigureMn = uimenu(app.PlotUIFigure,'Text','Options');

            % Manage toolbar
            app.FigureToolbar = uitoolbar(app.PlotUIFigure);

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

            % Manage tab groups
            FigureGl = uigridlayout(app.PlotUIFigure,[2,1],'Padding',[5,5,5,5]);
            FigureGl.RowHeight = {'1x','fit'};
            FigureGl.ColumnWidth = {'1x'};
            app.FigureTbGp = uitabgroup(FigureGl,'SelectedTab',app.ToolTb);
            app.FigureTbGp.Layout.Row = 1;
            app.FigureTbGp.Layout.Column = 1;
            
            % ------------------------------------------------------------------------
            % --------------------------Tool File Processing--------------------------
            % ------------------------------------------------------------------------

            app.ToolTb = uitab(app.FigureTbGp,'Title','Tool File Processing');
            app.ToolTb.ButtonDownFcn = createCallbackFcn(app,@ToolTbButtonDown,true);

            % Manage tool processing layout
            ToolTbGl = uigridlayout(app.ToolTb,[6,4]);
            ToolTbGl.RowHeight = {'fit','fit','fit','1x','fit','fit'};
            ToolTbGl.ColumnWidth = {'fit','1x','1x','1x'};
            % the elements of the cell array can be 'fit', fixed pixels, or '1x' '2x'
            %   'fit': to adjust the size to show the whole text, or basedon the default size
            %   fixed pixels: the size is fixed at the number of pixels
            %   variables ('2x'): to fill the remaining space, and the number is a
            %   weight for dividing up the remaining space among all variables
            
            % ------------------------Workspace directory------------------------
            WorkspaceDirGl = uigridlayout(ToolTbGl,[1,3]);
            WorkspaceDirGl.Layout.Row = 2;
            WorkspaceDirGl.Layout.Column = [1,4];
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
            app.WorkspaceDirBtn.ButtonPushedFcn = createCallbackFcn(app,@WorkspaceDirBtnValueChanged,true);
                 
            % ------------------------Tool data importing layout------------------------
            ToolFileGl = uigridlayout(ToolTbGl,[1,3]);
            ToolFileGl.Layout.Row = 3;
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
            app.ToolFileBtn.ButtonPushedFcn = createCallbackFcn(app,@ToolFileBtnValueChanged,true);

            % ------------------------Parameter editing tabgroup------------------------
            % Create the parameter selection tab group
            ParamPn = uipanel(ToolTbGl,'Title','Parameters', ...
                'TitlePosition','centertop','Scrollable','on');
            ParamPn.Layout.Row = [4,5];
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

            % Create the line fitting max distance label and dropdown box
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
            app.ResetBtn = uibutton(ParamPnGl,'push','Text','Reset','Visible','on');
            app.ResetBtn.Layout.Row = 2;
            app.ResetBtn.Layout.Column = 1;
            app.ResetBtn.ButtonPushedFcn = createCallbackFcn(app,@ResetBtnValueChanged,true);
            
            % ---Create the update button---
            app.UpdateBtn = uibutton(ParamPnGl,'push','Text','Update','Visible','on');
            app.UpdateBtn.Layout.Row = 2;
            app.UpdateBtn.Layout.Column = 2;
            app.UpdateBtn.ButtonPushedFcn = createCallbackFcn(app,@UpdateBtnValueChanged,true);

            % ------------------------Ploting axes------------------------
            app.ToolDataAxes = uiaxes(ToolTbGl);
            title(app.ToolDataAxes,'tool original data');
            % xlabel(app.UIAxes, 'X');
            % ylabel(app.UIAxes, 'Y');
            % zlabel(app.UIAxes, 'Z');
            app.ToolDataAxes.Layout.Row = 4;
            app.ToolDataAxes.Layout.Column = [2,4];

            app.PlotBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Extract & Plot');
            app.PlotBtn.Layout.Row = 5;
            app.PlotBtn.Layout.Column = 2;
            app.PlotBtn.ButtonPushedFcn = createCallbackFcn(app,@PlotBtnValueChanged,true);

            % ------------------------Ending process------------------------
            app.EnterBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Enter');
            app.EnterBtn.Layout.Row = 5;
            app.EnterBtn.Layout.Column = 3;
            app.EnterBtn.ButtonPushedFcn = createCallbackFcn(app,@EnterBtnValueChanged,true);
            
            app.CancelBtn = uibutton(ToolTbGl,'push','WordWrap','on','Text','Cancel');
            app.CancelBtn.Layout.Row = 5;
            app.CancelBtn.Layout.Column = 4;
            app.CancelBtn.ButtonPushedFcn = createCallbackFcn(app,@CancelBtnValueChanged,true);

            % ------------------------------------------------------------------------
            % ---------------------------Surface Processing---------------------------
            % ------------------------------------------------------------------------

            app.SurfaceTb = uitab(app.FigureTbGp,'Title','Surface Load');
            app.SurfaceTb.ButtonDownFcn = createCallbackFcn(app,@SurfaceTbButtonDown,true);

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
                app,@S1Tool2DBtnValueChanged,true);

            % button to execute s1_tool3D.m
            app.S1Tool3DBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_tool3D');
            app.S1Tool3DBtn.Layout.Row = 2;
            app.S1Tool3DBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1Tool3DBtnValueChanged,true);

            % button to execute s1_toolExtract_2Dline.m
            app.S1ToolExtract2DLineBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_2DLine');
            app.S1ToolExtract2DLineBtn.Layout.Row = 3;
            app.S1ToolExtract2DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract2DLineBtnValueChanged,true);

            % button to execute s1_toolExtract_3Dline.m
            app.S1ToolExtract3DLineBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_3Dline');
            app.S1ToolExtract3DLineBtn.Layout.Row = 4;
            app.S1ToolExtract3DLineBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtract3DLineBtnValueChanged,true);

            % button to execute s1_toolExtract_surf.m
            app.S1ToolExtractSurfBtn = uibutton(ToolFitGl,'push','WordWrap','on', ...
                'Text','s1_toolExtract_surf');
            app.S1ToolExtractSurfBtn.Layout.Row = 5;
            app.S1ToolExtractSurfBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S1ToolExtractSurfBtnValueChanged,true);

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
                app,@S2DesignSimulAsphericConcentricBtnValueChanged,true);

            % button to execute s2_design_simul_freeform.m
            app.S2DesignSimulFreeformBtn = uibutton(SimulateGl,'push', ...
                'WordWrap','on','Text','s2_design_simul_freeform');
            app.S2DesignSimulFreeformBtn.Layout.Column = 2;
            app.S2DesignSimulFreeformBtn.ButtonPushedFcn = createCallbackFcn( ...
                app,@S2DesignSimulFreeformBtnValueChanged,true);

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

            % ------------------------Info displaying window------------------------
            app.InfoTa = uitextarea(FigureGl,'WordWrap','on', ...
                'Editable','off','BackgroundColor',[0.96 0.96 0.96]);
            app.InfoTa.Layout.Row = 2;
            app.InfoTa.Layout.Column = 1;
            app.InfoTa.FontSize = 14;
            app.InfoTa.FontWeight = 'normal';
            app.InfoTa.Value = {'1 Welcome.'};
            app.InfoTa.ValueChangedFcn = createCallbackFcn(app,@InfoEfValueChanged,true);

            % Show the figure after all components are created
            app.PlotUIFigure.Visible = 'on';
            app.PlotUIFigure.WindowStyle = "normal";
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = upm_toolpath
            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app, app.PlotUIFigure);

            % Execute the startup function
            runStartupFcn(app, @startupFcn);

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.PlotUIFigure)
            % rmpath(genpath('funcs');
        end
    end
end
