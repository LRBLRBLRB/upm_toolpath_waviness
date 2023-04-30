function ToolFileBtnPushed(app,event)
    if isempty(app.workspaceDir)
        uialert(app.UIFigure, ...
            'Workspace directory doesn''t exist, please choose one first!', ...
            'Alert Message','CloseFcn',createCallbackFcn(app,@WorkspaceDirBtnPushed,true));
        return
    end
    [fileName,dirName] = uigetfile({ ...
        '*.csv','Comma-Separated-Values files(*.csv)'; ...
        '*.xls;*.xlsx','Excel worksheet files(*.xls,*.xlsx)'; ...
        '*.mat','MAT files(*.mat)'; ...
        '*.txt','Text files(*.txt)'; ...
        '*.*','all files(*.*)'...
        }, ...
        'Select One Tool Tip Measurement Data', ...
        app.workspaceDir, ...
        'MultiSelect','off');
    app.toolPathName = fullfile(dirName,fileName);
    app.ToolFileEf.Value = app.toolPathName;
    [~,~,fileExt] = fileparts(app.toolPathName);
    switch fileExt
        case {'.mat'}
            %% toolData has been generated previously
            % if the tool interpolation result is loaded, then plot it
            % directly; else get the original data from the file
            app.toolData = load(app.toolPathName);
            if isfield(app.toolData,'toolBform')
                app.toolUnit = app.toolData.unit;
                app.ToolUnitDd.Value = app.toolUnit;
                plot(app.ToolDataAxes, ...
                    app.toolData.toolFit(2,:),app.toolData.toolFit(3,:), ...
                    '--.','MarkerSize',8,'Color',[0,0.447,0.741]);
                hold(app.ToolDataAxes,'on');
                plot(app.ToolDataAxes, ...
                    app.toolData.toolBform.coefs(2,:),app.toolData.toolBform.coefs(3,:), ...
                    'x','Color',[0.32,0.55,0.19],'MarkerSize',5);
                plot(app.ToolDataAxes, ...
                    app.toolData.toolEdgePt(2,:),app.toolData.toolEdgePt(3,:), ...
                    'Color',[0.635,0.078,0.184]);
                axis(app.ToolDataAxes,'equal');
                set(app.ToolDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.ToolDataAxes,['y(',app.toolUnit,')']);
                ylabel(app.ToolDataAxes,['z(',app.toolUnit,')']);
                legend(app.ToolDataAxes, ...
                    'Measured Pts','Control Pts','Fitting Pts','Location','best');
                app.CheckToolLamp.Color = 'g';
                app.Msg = 'tool Data has been loaded.';
                InfoTaValueChanged(app,true);
%                     app.ParamTabToolfit.HandleVisibility = 'off';
%                     app.ParamTabToolmod.HandleVisibility = 'off';
%                     app.ParamTabToolinterp.HandleVisibility = 'off';

                app.Msg = 'Tool proces finished. Continue to load surface.';
            else
                app.toolData = [];
            end
        case {'.csv'}
            %% tool data file that has not been processed by mmt software
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
            if size(app.toolOri,2) ~= 3 && size(app.toolOri,2) ~= 2
                app.toolOri = app.toolOri.data;
            end
            app.toolOri(:,3) = [];
            app.toolOri = sortrows(app.toolOri,1,'ascend');
            app.toolOri = app.toolOri';

            app.Msg = 'Please select the parameters and click ''Update''.';
        case '.txt'
            %% tool data file that has been processed in mmt software
            app.toolOri = importdata(app.toolPathName,' ',0);
            app.toolOri = sortrows(app.toolOri,1,'ascend');
            app.toolOri = app.toolOri';
            
            app.Msg = 'Please select the parameters and click ''Update''.';
    end
    InfoTaValueChanged(app,true);
end