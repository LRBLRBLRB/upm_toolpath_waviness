function resetToolfitParams(app)
    app.ToolUnitDd.Value = app.toolUnitDefault;
    app.toolUnit = app.ToolUnitDd.Value;

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