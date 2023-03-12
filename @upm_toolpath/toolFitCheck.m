function toolFitCheck(app)
%TOOLFITCHECK to get 

    if isempty(app.toolDataFile)
        [fileName,dirName,fileType] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','all files(*.*)'}, ...
            'Select the fitted tool data file', ...
            app.workspaceDir, ...
            'MultiSelect','off');
        if ~fileType
            app.Msg = 'Tool tip isn''t fitted correctly!!!';
            InfoTaValueChanged(app,true);
            return;
        end
        app.toolDataFile = fullfile(dirName,fileName);
    end
    app.toolData = load(app.toolDataFile);
    app.RStepEf.Value = app.toolData.radius;
    app.Msg = 'Tool tip is interpolated successfully.';
    InfoTaValueChanged(app,true);
    app.CheckToolLamp.Color = 'g';
end