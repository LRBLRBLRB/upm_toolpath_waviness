function toolFitCheck(app,toolDataFile)
    if isempty(toolDataFile)
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
    else
        app.toolDataFile = toolDataFile;
    end
    app.toolData = load(app.toolDataFile);
    app.Msg = 'Tool tip is fitted successfully.';
    InfoTaValueChanged(app,true);
    app.CheckToolLamp.Color = 'g';
end