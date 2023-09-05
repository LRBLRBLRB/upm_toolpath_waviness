classdef add_surface < matlab.apps.AppBase
    % the figure to add surfaces

    properties (Constant)
        ButtonRowInterval       = 15
        ButtonColumnInterval    = 20
        ButtonWidth             = 120
        ButtonHeight            = 20
        DomainDefault           = 2.5
    end
    
    properties (Access = public)
        SurfaceUI               matlab.ui.Figure
        SurfaceGl               matlab.ui.container.GridLayout
        SurfaceDd               matlab.ui.control.DropDown
        ImportBtnGp             matlab.ui.container.ButtonGroup
        PointCloudRadioBtn      matlab.ui.control.RadioButton
        PointCloudBtnGl         matlab.ui.container.GridLayout
        PointCloudLb            matlab.ui.control.Label
        PointCloudEf            matlab.ui.control.EditField
        PointCloudActivateBtn   matlab.ui.control.Button
        FuncsRadioBtn           matlab.ui.control.RadioButton
        FuncsBtnGl              matlab.ui.container.GridLayout
        FuncsLb                 matlab.ui.control.Label
        FuncsEf                 matlab.ui.control.EditField
        FuncsActivateBtn        matlab.ui.control.Button
        Geometry2DBtnGp         matlab.ui.container.ButtonGroup
        ParaboloidRadioBtn      matlab.ui.control.RadioButton
        ParaboloidBtnGl         matlab.ui.container.GridLayout
        ParaboloidFuncLb        matlab.ui.control.Label
        ParaboloidActivateBtn   matlab.ui.control.Button
        ParaboloidALb           matlab.ui.control.Label
        ParaboloidASpin         matlab.ui.control.Spinner
        ParaboloidAUnitLb       matlab.ui.control.Label
        ParaboloidBLb           matlab.ui.control.Label
        ParaboloidBSpin         matlab.ui.control.Spinner
        AsphericRadioBtn        matlab.ui.control.RadioButton
        AsphericBtnGl           matlab.ui.container.GridLayout
        AsphericFuncLb          matlab.ui.control.Label
        AsphericActivateBtn     matlab.ui.control.Button
        AsphericCurvatureLb     matlab.ui.control.Label
        AsphericCurvatureSpin   matlab.ui.control.Spinner
        AsphericRadiusLb        matlab.ui.control.Label
        AsphericRadiusSpin      matlab.ui.control.Spinner
        AsphericRadiusUnitLb    matlab.ui.control.Label
        AsphericConicLb         matlab.ui.control.Label
        AsphericConicSpin       matlab.ui.control.Spinner
        AsphericOffsetLb        matlab.ui.control.Label
        AsphericOffsetSpin      matlab.ui.control.Spinner
        AsphericOffsetUnitLb    matlab.ui.control.Label
        Geometry3DBtnGp         matlab.ui.container.ButtonGroup
        EllipsoidRadioBtn       matlab.ui.control.RadioButton
        EllipsoidBtnGl          matlab.ui.container.GridLayout
        EllipsoidFuncLb         matlab.ui.control.Label
        EllipsoidActivateBtn    matlab.ui.control.Button
        EllipsoidALb            matlab.ui.control.Label
        EllipsoidASpin          matlab.ui.control.Spinner
        EllipsoidBLb            matlab.ui.control.Label
        EllipsoidBSpin          matlab.ui.control.Spinner
        EllipsoidCLb            matlab.ui.control.Label
        EllipsoidCSpin          matlab.ui.control.Spinner
        EllipsoidDLb            matlab.ui.control.Label
        EllipsoidDSpin          matlab.ui.control.Spinner
        BiconicRadioBtn         matlab.ui.control.RadioButton
        BiconicBtnGl            matlab.ui.container.GridLayout
        SurfaceDomainPn         matlab.ui.container.Panel
        SurfaceDomainPnGl       matlab.ui.container.GridLayout
        SurfaceDomainDd         matlab.ui.control.DropDown
        SurfaceXLb              matlab.ui.control.Label
        SurfaceX1Spin           matlab.ui.control.Spinner
        SurfaceX2Spin           matlab.ui.control.Spinner
        SurfaceXUnitLb          matlab.ui.control.Label
        SurfaceYLb              matlab.ui.control.Label
        SurfaceY1Spin           matlab.ui.control.Spinner
        SurfaceY2Spin           matlab.ui.control.Spinner
        SurfaceYUnitLb          matlab.ui.control.Label
        SurfaceRLb              matlab.ui.control.Label
        SurfaceR1Spin           matlab.ui.control.Spinner
        SurfaceR2Spin           matlab.ui.control.Spinner
        SurfaceRUnitLb          matlab.ui.control.Label
        SurfaceFuncsEf          matlab.ui.control.EditField
        SurfaceEnterBtn         matlab.ui.control.Button

        surfFuncs               
        surfPath
        unit
    end

    % relation to the main app
    properties (Access = private)
        CallingApp % main app object
    end

    methods (Access = private)
        function paramGpIni(app)
            app.PointCloudBtnGl.Visible = 'off';
            app.FuncsBtnGl.Visible = 'off';
            app.ParaboloidBtnGl.Visible = 'off';
            app.AsphericBtnGl.Visible = 'off';
            app.EllipsoidBtnGl.Visible = 'off';
        end
    end

    methods (Access = private)
        function startupFcn(app,mainapp,surfName,unit)
            % Store main app in property for CloseRequestFcn to use
            app.CallingApp = mainapp;

            % Update UI with the input values
            app.SurfaceDd.Value = surfName;
            SurfaceDdValueChanged(app,true);
            % switch surfName
            %     case 'Import'
            %         app.SurfaceDd.Value = 'Import';
            %         app.ImportBtnGp.Visible = 'on';
            %         app.PointCloudRadioBtn.Value = true;
            %         app.PointCloudBtnGl.Visible = 'on';
            %         app.PointCloudRadioBtn.Position = [app.ButtonColumnInterval, ...
            %             app.ImportBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
            %             app.ButtonWidth,app.ButtonHeight];
            %     case '2D Geometry'
            %         app.SurfaceDd.Value = '2D Geometry';
            %         app.Geometry2DBtnGp.Visible = 'on';
            %         app.AsphericRadioBtn.Value = true;
            %         app.AsphericBtnGl.Visible = 'on';
            %         app.AsphericRadioBtn.Position = [app.ButtonColumnInterval, ...
            %             app.Geometry2DBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
            %             app.ButtonWidth,app.ButtonHeight];
            %     case '3D Geometry'
            %         app.SurfaceDd.Value = '3D Geometry';
            %         app.Geometry3DBtnGp.Visible = 'on';
            %         app.EllipsoidRadioBtn.Value = true;
            %         app.EllipsoidBtnGl.Visible = 'on';
            %         app.FuncsRadioBtn.Position = [app.ButtonColumnInterval, ...
            %             app.Geometry3DBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
            %             app.ButtonWidth,app.ButtonHeight];
            %         app.EllipsoidRadioBtn.Position = [app.ButtonColumnInterval, ...
            %             app.Geometry3DBtnGp.Position(4) - 2*app.ButtonHeight - 2*app.ButtonRowInterval, ...
            %             app.ButtonWidth,app.ButtonHeight];
            % end
            
            app.unit = unit;
            if strcmp(app.unit,'\mum')
                app.unit = '$\mu$m';
            end
            app.AsphericRadiusUnitLb.Text = app.unit;
            app.AsphericOffsetUnitLb.Text = app.unit;
            app.SurfaceXUnitLb.Text = app.unit;
            app.SurfaceYUnitLb.Text = app.unit;
            app.SurfaceRUnitLb.Text = app.unit;
        end

        % uifigure closed request callback
        function AddSurfaceUICloseReq(app,event)
            app.CallingApp.AddSurfaceBtn.Enable = 'on';
            delete(app);
        end

        % Value changed function: select the proper surface category
        function SurfaceDdValueChanged(app,event)
            switch app.SurfaceDd.Value
                case 'Import'
                    app.ImportBtnGp.Visible = 'on';
                    app.Geometry2DBtnGp.Visible = 'off';
                    app.Geometry3DBtnGp.Visible = 'off';
                    paramGpIni(app);
                    app.PointCloudRadioBtn.Value = true;
                    app.PointCloudBtnGl.Visible = 'on';
                    app.PointCloudRadioBtn.Position = [app.ButtonColumnInterval, ...
                        app.ImportBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
                        app.ButtonWidth,app.ButtonHeight];
                    app.SurfaceDomainDd.Enable = 'on';
                case '2D Geometry'
                    app.ImportBtnGp.Visible = 'off';
                    app.Geometry2DBtnGp.Visible = 'on';
                    app.Geometry3DBtnGp.Visible = 'off';
                    paramGpIni(app);
                    app.ParaboloidRadioBtn.Value = true;
                    app.ParaboloidBtnGl.Visible = 'on';
                    app.ParaboloidRadioBtn.Position = [app.ButtonColumnInterval, ...
                        app.Geometry2DBtnGp.Position(4) - 2*app.ButtonHeight - app.ButtonRowInterval, ...
                        app.ButtonWidth,2*app.ButtonHeight];
                    app.AsphericRadioBtn.Position = [app.ButtonColumnInterval, ...
                        app.Geometry2DBtnGp.Position(4) - 3*app.ButtonHeight - 2*app.ButtonRowInterval, ...
                        app.ButtonWidth,app.ButtonHeight];
                    app.SurfaceDomainDd.Value = 'Polar';
                    app.SurfaceDomainDd.Enable = 'off';
                case '3D Geometry'
                    app.ImportBtnGp.Visible = 'off';
                    app.Geometry2DBtnGp.Visible = 'off';
                    app.Geometry3DBtnGp.Visible = 'on';
                    paramGpIni(app);
                    app.EllipsoidRadioBtn.Value = true;
                    app.EllipsoidBtnGl.Visible = 'on';
                    app.FuncsRadioBtn.Position = [app.ButtonColumnInterval, ...
                        app.Geometry3DBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
                        app.ButtonWidth,app.ButtonHeight];
                    app.EllipsoidRadioBtn.Position = [app.ButtonColumnInterval, ...
                        app.Geometry3DBtnGp.Position(4) - 2*app.ButtonHeight - 2*app.ButtonRowInterval, ...
                        app.ButtonWidth,app.ButtonHeight];
                    app.SurfaceDomainDd.Enable = 'on';
            end
        end

        % Selection changed function: change the radio buttion of the import buttongroup
        function ImportBtnGpSelectionChanged(app,true)
            switch app.ImportBtnGp.SelectedObject.Text
                case 'Point Cloud'
                    paramGpIni(app);
                    app.PointCloudBtnGl.Visible = 'on';
            end
        end

        % Button pushed function: change the add_functions edit field
        function PointCloudActivateButtonPushed(app,event)
            app.SurfaceFuncsEf.Value = app.PointCloudEf.Value;
        end

        % Selection changed function: change the radio buttion of the 2-D geometry buttongroup
        function Geometry2DBtnGpSelectionChanged(app,event)
            switch app.Geometry2DBtnGp.SelectedObject.Text
                case 'Rotating Paraboloid'
                    paramGpIni(app);
                    app.ParaboloidBtnGl.Visible = 'on';
                case 'Aspheric'
                    paramGpIni(app);
                    app.AsphericBtnGl.Visible = 'on';
            end
        end

        % Value changed function: change the spinner of the aspheric surface
        function ParaboloidActivateBtnButtonPushed(app,event)
            syms x y;
            a = app.ParaboloidASpin.Value;
            b = app.ParaboloidBSpin.Value;
            app.surfFuncs = a*(x^2 + y^2) + b;
            app.SurfaceFuncsEf.Value = char(app.surfFuncs);
        end

        function AsphericCurvatureSpinValueChanged(app,event)
            app.AsphericRadiusSpin.Value = 1/app.AsphericCurvatureSpin.Value;
        end

        function AsphericRadiusSpinValueChanged(app,event)
            app.AsphericCurvatureSpin.Value = 1/app.AsphericRadiusSpin.Value;
        end

        % Value changed function: change the spinner of the aspheric surface
        function AsphericActivateBtnButtonPushed(app,event)
            syms x y;
            c = app.AsphericCurvatureSpin.Value;
            k = app.AsphericConicSpin.Value;
            x0 = app.AsphericOffsetSpin.Value;
            app.surfFuncs = (c*((x - x0)^2 + (y - x0)^2))/(1 + sqrt(1 - (1 + k)*c^2*((x - x0)^2 + (y - x0)^2)));
            app.SurfaceFuncsEf.Value = char(app.surfFuncs);
        end
    
        % Selection changed function: change the radio buttion of the 3-D geometry buttongroup
        function Geometry3DBtnGpSelectionChanged(app,event)
            switch app.Geometry3DBtnGp.SelectedObject.Text
                case 'Function-Based'
                    paramGpIni(app);
                    app.FuncsBtnGl.Visible = 'on';
                case 'Ellipsoid'
                    paramGpIni(app);
                    app.EllipsoidBtnGl.Visible = 'on';
            end
        end

        % Button pushed function: change the add_functions edit field
        function FuncsActivateButtonPushed(app,event)
            app.SurfaceFuncsEf.Value = app.FuncsEf.Value;
            app.surfFuncs
        end

        % Value changed function: change the spinner of the Ellipsoid surface
        function EllipsoidActivateBtnButtonPushed(app,event)
            syms x y;
            A = app.EllipsoidASpin.Value;
            B = app.EllipsoidBSpin.Value;
            C = app.EllipsoidCSpin.Value;
            D = app.EllipsoidDSpin.Value;
            app.surfFuncs = C*sqrt(D.^2 - x.^2/A^2 - y.^2/B^2);
            app.SurfaceFuncsEf.Value = char(app.surfFuncs);
        end

        % Value changed funciton: change the domain type
        function SurfaceDomainDdValueChanged(app,event)
            switch app.SurfaceDomainDd.Value
                case 'XY'
                    app.SurfaceXLb.Visible = 'on';
                    app.SurfaceX1Spin.Visible = 'on';
                    app.SurfaceX2Spin.Visible = 'on';
                    app.SurfaceXUnitLb.Visible = 'on';
                    app.SurfaceYLb.Visible = 'on';
                    app.SurfaceY1Spin.Visible = 'on';
                    app.SurfaceY2Spin.Visible = 'on';
                    app.SurfaceYUnitLb.Visible = 'on';
                    app.SurfaceRLb.Visible = 'off';
                    app.SurfaceR1Spin.Visible = 'off';
                    app.SurfaceR2Spin.Visible = 'off';
                    app.SurfaceRUnitLb.Visible = 'off';
                case 'Ploar'
                    app.SurfaceXLb.Visible = 'off';
                    app.SurfaceX1Spin.Visible = 'off';
                    app.SurfaceX2Spin.Visible = 'off';
                    app.SurfaceXUnitLb.Visible = 'off';
                    app.SurfaceYLb.Visible = 'off';
                    app.SurfaceY1Spin.Visible = 'off';
                    app.SurfaceY2Spin.Visible = 'off';
                    app.SurfaceYUnitLb.Visible = 'off';
                    app.SurfaceRLb.Visible = 'on';
                    app.SurfaceR1Spin.Visible = 'on';
                    app.SurfaceR2Spin.Visible = 'on';
                    app.SurfaceRUnitLb.Visible = 'on';
            end
        end

        % Value changed function: ensure R1 and R2 to be opposite numbers
        function SurfaceR1SpinValueChanged(app,true)
            app.SurfaceR2Spin.Value = -app.SurfaceR1Spin.Value;
        end
        function SurfaceR2SpinValueChanged(app,true)
            app.SurfaceR1Spin.Value = -app.SurfaceR2Spin.Value;
        end

        % Button pushed function: correct the surface and sent it back
        function SurfaceEnterButtonPushed(app,event)
            switch app.SurfaceDomainDd.Value
                case 'XY'
                    surfDomain = [app.SurfaceX1Spin.Value,app.SurfaceX2Spin.Value;
                        app.SurfaceY1Spin.Value,app.SurfaceY2Spin.Value];
                case 'Polar'
                    surfDomain = [app.SurfaceR1Spin.Value,app.SurfaceR2Spin.Value];
            end
            switch app.SurfaceDd.Value
                case 'Import'
                    % Call main app's public function
                    updateSurface(app.CallingApp, ...
                        app.ImportBtnGp.SelectedObject.Text, ...
                        app.surfPath.Value,surfDomain);
                case '2D Geometry'
                    % Call main app's public function
                    updateSurface(app.CallingApp, ...
                        app.Geometry2DBtnGp.SelectedObject.Text, ...
                        app.surfFuncs,surfDomain);
                case '3D Geometry'
                    % Call main app's public function
                    updateSurface(app.CallingApp, ...
                        app.Geometry3DBtnGp.SelectedObject.Text, ...
                        app.surfFuncs,surfDomain);
            end

            % Delete the dialog box
            delete(app);
        end

        function createComponents(app)
            % Create the add surface window
            app.SurfaceUI = uifigure('Name','Add Surface', ...
                'WindowStyle','modal','WindowState','normal','Visible','off');
            app.SurfaceUI.CloseRequestFcn = createCallbackFcn(app,@AddSurfaceUICloseReq,true);
            app.SurfaceUI.Resize = "on";
            app.SurfaceUI.Position = [600,600,600,400];
            app.SurfaceUI.Scrollable = "on";
        
            app.SurfaceGl = uigridlayout(app.SurfaceUI,[4,3]);
            app.SurfaceGl.RowHeight = {'fit','2x','fit','fit'};
            app.SurfaceGl.ColumnWidth = {app.ButtonWidth + 2*app.ButtonColumnInterval,'2x','fit'};
        
            app.SurfaceDd = uidropdown(app.SurfaceGl, ...
                'Items',{'Import','2D Geometry','3D Geometry'}, ...
                'Value','3D Geometry','BackgroundColor',[1,1,1]);
            app.SurfaceDd.Layout.Row = 1;
            app.SurfaceDd.Layout.Column = 1;
            app.SurfaceDd.ValueChangedFcn = createCallbackFcn(app,@SurfaceDdValueChanged,true);

            % ------------------------------------------------------------------
            % ------------------------- Import Surface -------------------------
            % ------------------------------------------------------------------
            app.ImportBtnGp = uibuttongroup(app.SurfaceGl,'BorderType','line', ...
                'Title','Import Selection','TitlePosition','centertop', ...
                'FontSize',16,'FontName','Times New Roman');
            app.ImportBtnGp.Layout.Row = [2,3];
            app.ImportBtnGp.Layout.Column = 1;
            app.ImportBtnGp.Visible = 'off';
            app.ImportBtnGp.SelectionChangedFcn = createCallbackFcn(app,@ImportBtnGpSelectionChanged,true);

            % ----------------------- Point-cloud surface -----------------------
            app.PointCloudRadioBtn = uiradiobutton(app.ImportBtnGp,'Text','Point Cloud','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval,app.ImportBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
                app.ButtonWidth,app.ButtonHeight]);
            app.PointCloudBtnGl = uigridlayout(app.SurfaceGl,[2,2]);
            app.PointCloudBtnGl.Layout.Row = [1,2];
            app.PointCloudBtnGl.Layout.Column = [2,3];
            app.PointCloudBtnGl.RowHeight = {'fit','1x'};
            app.PointCloudBtnGl.ColumnWidth = {'1x','fit'};
            app.PointCloudBtnGl.Visible = 'off';

            app.PointCloudLb = uilabel(app.PointCloudBtnGl,'Text','Surface path: ');
            app.PointCloudLb.Layout.Row = 1;
            app.PointCloudLb.Layout.Column = 1;
            app.PointCloudEf = uieditfield(app.PointCloudBtnGl,'text');
            app.PointCloudEf.Layout.Row = 2;
            app.PointCloudEf.Layout.Column = [1,2];

            app.PointCloudActivateBtn = uibutton(app.PointCloudBtnGl,'push','Text','Activate');
            app.PointCloudActivateBtn.Layout.Row = 1;
            app.PointCloudActivateBtn.Layout.Column = 2;
            app.PointCloudActivateBtn.ButtonPushedFcn = createCallbackFcn(app,@PointCloudActivateButtonPushed,true);

            % -------------------------------------------------------------------
            % ----------------------- 2D Geometry Surface -----------------------
            % -------------------------------------------------------------------
            app.Geometry2DBtnGp = uibuttongroup(app.SurfaceGl,'Title','','BorderType','line');
            app.Geometry2DBtnGp.Layout.Row = [2,3];
            app.Geometry2DBtnGp.Layout.Column = 1;
            app.Geometry2DBtnGp.Visible = 'off';
            app.Geometry2DBtnGp.SelectionChangedFcn = createCallbackFcn(app,@Geometry2DBtnGpSelectionChanged,true);

            % ----------------------- Rotating Paraboloid -----------------------
            app.ParaboloidRadioBtn = uiradiobutton(app.Geometry2DBtnGp,'Text','Rotating Paraboloid','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval, ...
                app.Geometry2DBtnGp.Position(4) - 2*app.ButtonHeight - app.ButtonRowInterval, ...
                app.ButtonWidth,2*app.ButtonHeight]);
            app.ParaboloidBtnGl = uigridlayout(app.SurfaceGl,[3,2]);
            app.ParaboloidBtnGl.Layout.Row = [1,2];
            app.ParaboloidBtnGl.Layout.Column = [2,3];
            app.ParaboloidBtnGl.RowHeight = {'fit','fit','fit'};
            app.ParaboloidBtnGl.ColumnWidth = {'fit','1x','fit'};
            app.ParaboloidBtnGl.Visible = 'off';

            app.ParaboloidFuncLb = uilabel(app.ParaboloidBtnGl,'Interpreter','latex');
            app.ParaboloidFuncLb.Layout.Row = 1;
            app.ParaboloidFuncLb.Layout.Column = [1,2];
            app.ParaboloidFuncLb.Text = '$Z = a\cdot r^{2} + b$';

            app.ParaboloidActivateBtn = uibutton(app.ParaboloidBtnGl,'push','Text','Activate');
            app.ParaboloidActivateBtn.Layout.Row = 1;
            app.ParaboloidActivateBtn.Layout.Column = 3;
            app.ParaboloidActivateBtn.ButtonPushedFcn = createCallbackFcn(app,@ParaboloidActivateBtnButtonPushed,true);

            app.ParaboloidALb = uilabel(app.ParaboloidBtnGl,'Text','a');
            app.ParaboloidALb.Layout.Row = 2;
            app.ParaboloidALb.Layout.Column = 1;
            app.ParaboloidASpin = uispinner(app.ParaboloidBtnGl,'Value',0.091/1000,'Limits',[-inf,inf]);
            app.ParaboloidASpin.Layout.Row = 2;
            app.ParaboloidASpin.Layout.Column = [2,3];

            app.ParaboloidBLb = uilabel(app.ParaboloidBtnGl,'Text','b');
            app.ParaboloidBLb.Layout.Row = 3;
            app.ParaboloidBLb.Layout.Column = 1;
            app.ParaboloidBSpin = uispinner(app.ParaboloidBtnGl,'Value',0,'Limits',[-inf,inf]);
            app.ParaboloidBSpin.Layout.Row = 3;
            app.ParaboloidBSpin.Layout.Column = [2,3];

            % ----------------------- Aspheric -----------------------
            app.AsphericRadioBtn = uiradiobutton(app.Geometry2DBtnGp,'Text','Aspheric','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval, ...
                app.Geometry2DBtnGp.Position(4) - 3*app.ButtonHeight - 2*app.ButtonRowInterval, ...
                app.ButtonWidth,app.ButtonHeight]);
            app.AsphericBtnGl = uigridlayout(app.SurfaceGl,[4,3]);
            app.AsphericBtnGl.Layout.Row = [1,2];
            app.AsphericBtnGl.Layout.Column = [2,3];
            app.AsphericBtnGl.RowHeight = {'fit','fit','fit','fit'};
            app.AsphericBtnGl.ColumnWidth = {'fit','1x','fit'};
            app.AsphericBtnGl.Visible = 'off';

            app.AsphericFuncLb = uilabel(app.AsphericBtnGl,'Interpreter','latex');
            app.AsphericFuncLb.Layout.Row = 1;
            app.AsphericFuncLb.Layout.Column = [1,2];
            app.AsphericFuncLb.Text = '$Z = \frac{c(r - r_{0})^2}{1 + \sqrt{1 - (1 + k)c^2(r - r_{0})^2}}$';

            app.AsphericActivateBtn = uibutton(app.AsphericBtnGl,'push','Text','Activate');
            app.AsphericActivateBtn.Layout.Row = 1;
            app.AsphericActivateBtn.Layout.Column = 3;
            app.AsphericActivateBtn.ButtonPushedFcn = createCallbackFcn(app,@AsphericActivateBtnButtonPushed,true);

            app.AsphericCurvatureLb = uilabel(app.AsphericBtnGl,'Text','Curvature (c)');
            app.AsphericCurvatureLb.Layout.Row = 2;
            app.AsphericCurvatureLb.Layout.Column = 1;
            app.AsphericCurvatureSpin = uispinner(app.AsphericBtnGl,'Value',0.69/1000,'Limits',[0,inf]);
            app.AsphericCurvatureSpin.Layout.Row = 2;
            app.AsphericCurvatureSpin.Layout.Column = 2;
            app.AsphericCurvatureSpin.ValueChangedFcn = createCallbackFcn(app,@AsphericCurvatureSpinValueChanged,true);

            app.AsphericRadiusLb = uilabel(app.AsphericBtnGl,'Text','Radius (1/c)');
            app.AsphericRadiusLb.Layout.Row = 3;
            app.AsphericRadiusLb.Layout.Column = 1;
            app.AsphericRadiusSpin = uispinner(app.AsphericBtnGl,'Value',1000/0.69,'Limits',[0,inf]);
            app.AsphericRadiusSpin.Layout.Row = 3;
            app.AsphericRadiusSpin.Layout.Column = 2;
            app.AsphericRadiusSpin.ValueChangedFcn = createCallbackFcn(app,@AsphericRadiusSpinValueChanged,true);
            app.AsphericRadiusUnitLb = uilabel(app.AsphericBtnGl,'Interpreter','latex');
            app.AsphericRadiusUnitLb.Layout.Row = 3;
            app.AsphericRadiusUnitLb.Layout.Column = 3;

            app.AsphericConicLb = uilabel(app.AsphericBtnGl,'Text','Conic (k)');
            app.AsphericConicLb.Layout.Row = 4;
            app.AsphericConicLb.Layout.Column = 1;
            app.AsphericConicSpin = uispinner(app.AsphericBtnGl,'Value',0,'Limits',[0,inf]);
            app.AsphericConicSpin.Layout.Row = 4;
            app.AsphericConicSpin.Layout.Column = 2;

            app.AsphericOffsetLb = uilabel(app.AsphericBtnGl,'Text','X Offset (x_0)');
            app.AsphericOffsetLb.Layout.Row = 5;
            app.AsphericOffsetLb.Layout.Column = 1;
            app.AsphericOffsetSpin = uispinner(app.AsphericBtnGl,'Value',0,'Limits',[-inf,inf]);
            app.AsphericOffsetSpin.Layout.Row = 5;
            app.AsphericOffsetSpin.Layout.Column = 2;
            app.AsphericOffsetUnitLb = uilabel(app.AsphericBtnGl,'Interpreter','latex');
            app.AsphericOffsetUnitLb.Layout.Row = 5;
            app.AsphericOffsetUnitLb.Layout.Column = 3;

            % -------------------------------------------------------------------
            % ----------------------- 3D Geometry Surface -----------------------
            % -------------------------------------------------------------------
            app.Geometry3DBtnGp = uibuttongroup(app.SurfaceGl,'Title','','BorderType','line');
            app.Geometry3DBtnGp.Layout.Row = [2,3];
            app.Geometry3DBtnGp.Layout.Column = 1;
            app.Geometry3DBtnGp.Visible = 'off';
            app.Geometry3DBtnGp.SelectionChangedFcn = createCallbackFcn(app,@Geometry3DBtnGpSelectionChanged,true);

            % ----------------------- Function-based surface -----------------------
            app.FuncsRadioBtn = uiradiobutton(app.Geometry3DBtnGp,'Text','Function-Based','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval, ...
                app.Geometry3DBtnGp.Position(4) - app.ButtonHeight - app.ButtonRowInterval, ...
                app.ButtonWidth,app.ButtonHeight]);
            app.FuncsBtnGl = uigridlayout(app.SurfaceGl,[2,2]);
            app.FuncsBtnGl.Layout.Row = [1,2];
            app.FuncsBtnGl.Layout.Column = [2,3];
            app.FuncsBtnGl.RowHeight = {'fit','1x'};
            app.FuncsBtnGl.ColumnWidth = {'1x','fit'};
            app.FuncsBtnGl.Visible = 'off';

            app.FuncsLb = uilabel(app.FuncsBtnGl,'Text','Function: ');
            app.FuncsLb.Layout.Row = 1;
            app.FuncsLb.Layout.Column = 1;
            app.FuncsEf = uieditfield(app.FuncsBtnGl,'text');
            app.FuncsEf.Layout.Row = 2;
            app.FuncsEf.Layout.Column = [1,2];
            app.FuncsEf.Value = 'C*sqrt(R.^2 - x.^2/A^2 - y.^2/B^2)';

            app.FuncsActivateBtn = uibutton(app.FuncsBtnGl,'push','Text','Activate');
            app.FuncsActivateBtn.Layout.Row = 1;
            app.FuncsActivateBtn.Layout.Column = 2;
            app.FuncsActivateBtn.ButtonPushedFcn = createCallbackFcn(app,@FuncsActivateButtonPushed,true);

            % ----------------------- Elliptical surface -----------------------
            app.EllipsoidRadioBtn = uiradiobutton(app.Geometry3DBtnGp,'Text','Ellipsoid','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval, ...
                app.Geometry3DBtnGp.Position(4) - 2*app.ButtonHeight - 2*app.ButtonRowInterval, ...
                app.ButtonWidth,app.ButtonHeight]);
            app.EllipsoidBtnGl = uigridlayout(app.SurfaceGl,[5,4]);
            app.EllipsoidBtnGl.Layout.Row = [1,2];
            app.EllipsoidBtnGl.Layout.Column = [2,3];
            app.EllipsoidBtnGl.RowHeight = {'fit','fit','fit','fit','fit'};
            app.EllipsoidBtnGl.ColumnWidth = {'fit','1x','1x','fit'};
            app.EllipsoidBtnGl.Visible = 'off';
            app.EllipsoidBtnGl.Scrollable = 'on';

            app.EllipsoidFuncLb = uilabel(app.EllipsoidBtnGl,'Interpreter','latex');
            app.EllipsoidFuncLb.Layout.Row = 1;
            app.EllipsoidFuncLb.Layout.Column = [1,3];
            app.EllipsoidFuncLb.Text = '$Z = \sqrt{C^{2}(D-\frac{x^{2}}{A^{2}}-\frac{y^{2}}{B^{2}})}$';

            app.EllipsoidActivateBtn = uibutton(app.EllipsoidBtnGl,'push','Text','Activate');
            app.EllipsoidActivateBtn.Layout.Row = 1;
            app.EllipsoidActivateBtn.Layout.Column = 4;
            app.EllipsoidActivateBtn.ButtonPushedFcn = createCallbackFcn(app,@EllipsoidActivateBtnButtonPushed,true);

            app.EllipsoidALb = uilabel(app.EllipsoidBtnGl,'Text','A');
            app.EllipsoidALb.Layout.Row = 2;
            app.EllipsoidALb.Layout.Column = 1;
            app.EllipsoidASpin = uispinner(app.EllipsoidBtnGl,'Value',3.5/2);
            app.EllipsoidASpin.Layout.Row = 2;
            app.EllipsoidASpin.Layout.Column = [2,4];

            app.EllipsoidBLb = uilabel(app.EllipsoidBtnGl,'Text','B');
            app.EllipsoidBLb.Layout.Row = 3;
            app.EllipsoidBLb.Layout.Column = 1;
            app.EllipsoidBSpin = uispinner(app.EllipsoidBtnGl,'Value',4/2);
            app.EllipsoidBSpin.Layout.Row = 3;
            app.EllipsoidBSpin.Layout.Column = [2,4];

            app.EllipsoidCLb = uilabel(app.EllipsoidBtnGl,'Text','C');
            app.EllipsoidCLb.Layout.Row = 4;
            app.EllipsoidCLb.Layout.Column = 1;
            app.EllipsoidCSpin = uispinner(app.EllipsoidBtnGl,'Value',5/2);
            app.EllipsoidCSpin.Layout.Row = 4;
            app.EllipsoidCSpin.Layout.Column = [2,4];

            app.EllipsoidDLb = uilabel(app.EllipsoidBtnGl,'Text','R');
            app.EllipsoidDLb.Layout.Row = 5;
            app.EllipsoidDLb.Layout.Column = 1;
            app.EllipsoidDSpin = uispinner(app.EllipsoidBtnGl,'Value',10/2*1000);
            app.EllipsoidDSpin.Layout.Row = 5;
            app.EllipsoidDSpin.Layout.Column = [2,4];

            % ----------------------- biconic surface -----------------------
            app.BiconicRadioBtn = uiradiobutton(app.Geometry3DBtnGp,'Text','Biconic','WordWrap','on', ...
                'Position',[app.ButtonColumnInterval, ...
                app.Geometry3DBtnGp.Position(4) - 3*app.ButtonHeight - 3*app.ButtonRowInterval, ...
                app.ButtonWidth,app.ButtonHeight]);
            app.BiconicBtnGl = uigridlayout(app.SurfaceGl,[5,4]);
            app.BiconicBtnGl.Layout.Row = [1,2];
            app.BiconicBtnGl.Layout.Column = [2,3];
            app.BiconicBtnGl.RowHeight = {'fit','fit','fit','fit','fit'};
            app.BiconicBtnGl.ColumnWidth = {'fit','1x','1x','fit'};
            app.BiconicBtnGl.Visible = 'off';
            app.BiconicBtnGl.Scrollable = 'on';

            % ----------------------- surface domain -----------------------
            app.SurfaceDomainPn = uipanel(app.SurfaceGl,'Title','Domain', ...
                'TitlePosition','centertop');
            app.SurfaceDomainPn.Layout.Row = 3;
            app.SurfaceDomainPn.Layout.Column = [2,3];
            app.SurfaceDomainPnGl = uigridlayout(app.SurfaceDomainPn,[3,4]);
            app.SurfaceDomainPnGl.RowHeight = {'fit','fit','fit'};
            app.SurfaceDomainPnGl.ColumnWidth = {'fit','2x','2x','1x'};
            app.SurfaceDomainDd = uidropdown(app.SurfaceDomainPnGl, ...
                'Items',{'XY','Polar'},'Value','Polar');
            app.SurfaceDomainDd.Layout.Row = 1;
            app.SurfaceDomainDd.Layout.Column = [1,4];
            app.SurfaceDomainDd.ValueChangedFcn = createCallbackFcn(app,@SurfaceDomainDdValueChanged,true);
            
            app.SurfaceXLb = uilabel(app.SurfaceDomainPnGl,'Text','X Range');
            app.SurfaceXLb.Layout.Row = 2;
            app.SurfaceXLb.Layout.Column = 1;
            app.SurfaceXLb.Visible = 'off';
            app.SurfaceX1Spin = uispinner(app.SurfaceDomainPnGl,'Value',-app.DomainDefault);
            app.SurfaceX1Spin.Layout.Row = 2;
            app.SurfaceX1Spin.Layout.Column = 2;
            app.SurfaceX1Spin.Visible = 'off';
            app.SurfaceX2Spin = uispinner(app.SurfaceDomainPnGl,'Value',app.DomainDefault);
            app.SurfaceX2Spin.Layout.Row = 2;
            app.SurfaceX2Spin.Layout.Column = 3;
            app.SurfaceX2Spin.Visible = 'off';
            app.SurfaceXUnitLb = uilabel(app.SurfaceDomainPnGl,'Interpreter','latex');
            app.SurfaceXUnitLb.Layout.Row = 2;
            app.SurfaceXUnitLb.Layout.Column = 4;
            app.SurfaceXUnitLb.Visible = 'off';

            app.SurfaceYLb = uilabel(app.SurfaceDomainPnGl,'Text','Y Range');
            app.SurfaceYLb.Layout.Row = 3;
            app.SurfaceYLb.Layout.Column = 1;
            app.SurfaceYLb.Visible = 'off';
            app.SurfaceY1Spin = uispinner(app.SurfaceDomainPnGl,'Value',-app.DomainDefault);
            app.SurfaceY1Spin.Layout.Row = 3;
            app.SurfaceY1Spin.Layout.Column = 2;
            app.SurfaceY1Spin.Visible = 'off';
            app.SurfaceY2Spin = uispinner(app.SurfaceDomainPnGl,'Value',app.DomainDefault);
            app.SurfaceY2Spin.Layout.Row = 3;
            app.SurfaceY2Spin.Layout.Column = 3;
            app.SurfaceY2Spin.Visible = 'off';
            app.SurfaceYUnitLb = uilabel(app.SurfaceDomainPnGl,'Interpreter','latex');
            app.SurfaceYUnitLb.Layout.Row = 3;
            app.SurfaceYUnitLb.Layout.Column = 4;
            app.SurfaceYUnitLb.Visible = 'off';

            app.SurfaceRLb = uilabel(app.SurfaceDomainPnGl,'Text','Y Range');
            app.SurfaceRLb.Layout.Row = 2;
            app.SurfaceRLb.Layout.Column = 2;
            app.SurfaceR1Spin = uispinner(app.SurfaceDomainPnGl,'Value',-app.DomainDefault);
            app.SurfaceR1Spin.Layout.Row = 2;
            app.SurfaceR1Spin.Layout.Column = 2;
            app.SurfaceR1Spin.ValueChangedFcn = createCallbackFcn(app,@SurfaceR1SpinValueChanged,true);
            app.SurfaceR2Spin = uispinner(app.SurfaceDomainPnGl,'Value',app.DomainDefault);
            app.SurfaceR2Spin.Layout.Row = 2;
            app.SurfaceR2Spin.Layout.Column = 3;
            app.SurfaceR2Spin.ValueChangedFcn = createCallbackFcn(app,@SurfaceR2SpinValueChanged,true);
            app.SurfaceRUnitLb = uilabel(app.SurfaceDomainPnGl,'Interpreter','latex');
            app.SurfaceRUnitLb.Layout.Row = 2;
            app.SurfaceRUnitLb.Layout.Column = 4;

            % ------------------------ state text area ------------------------
            app.SurfaceFuncsEf = uieditfield(app.SurfaceGl,'text', ...
                'Editable','off','BackgroundColor',[0.96,0.96,0.96]);
            app.SurfaceFuncsEf.Layout.Row = 4;
            app.SurfaceFuncsEf.Layout.Column = [1,2];
            app.SurfaceFuncsEf.Value = 'No surface is activated now.';

            app.SurfaceEnterBtn = uibutton(app.SurfaceGl,'push','Text','Enter');
            app.SurfaceEnterBtn.Layout.Row = 4;
            app.SurfaceEnterBtn.Layout.Column = 3;
            app.SurfaceEnterBtn.ButtonPushedFcn = createCallbackFcn(app,@SurfaceEnterButtonPushed,true);

            % enable the window
            app.SurfaceUI.Visible = 'on';
        end
    end

    methods (Access = public)
        function app = add_surface(varargin)
            % Create UIFigure and components
            createComponents(app);

            % Register the app with App Designer
            registerApp(app,app.SurfaceUI);

            % Execute the startup function
            runStartupFcn(app,@(app)startupFcn(app,varargin{:}));

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.SurfaceUI);
        end
    end
end