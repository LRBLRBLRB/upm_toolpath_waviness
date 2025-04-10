% to reset the optimization tab paramaters
function resetOptimParams(app)
    app.CutDirectionDd.Value = app.cutDirectionDefault;
    app.cutDirection = app.CutDirectionDd.Value;
    app.StartDirectionDd.Value = app.startDirectionDefault;
    app.startDirection = app.StartDirectionDd.Value;
    app.AngularDiscreteDd.Value = app.angularDiscreteDefault;
    app.angularDiscrete = app.AngularDiscreteDd.Value;
    app.AimResEf.Value = app.aimResDefault;
    app.aimRes = app.AimResEf.Value;
    app.MaxIterSpin.Value = app.maxIterDefault;
    app.maxIter = app.MaxIterSpin.Value;
    app.RStepEf.Value = app.rStepDefault;
    app.rStep = app.RStepEf.Value;
    app.ArcLengthEf.Value = app.arcLengthDefault;
    app.arcLength = app.ArcLengthEf.Value;
    app.MaxAngPtDistEf.Value = app.maxAngPtDistDefault;
    app.maxAngPtDist = app.MaxAngPtDistEf.Value;
    app.AngularLengthEf.Value = app.angularLengthDefault;
    app.angularLength = app.AngularLengthEf.Value;
    app.RadialIncrementDd.Value = app.radialIncrementDefault;
    app.radialIncrement = app.RadialIncrementDd.Value;
    app.SpiralMethodDd.Value = app.spiralMethodDefault;
    app.spiralMethod = app.SpiralMethodDd.Value;
    app.ZAllowanceEf.Value = app.zAllowanceDefault;
    app.zAllowance = app.ZAllowanceEf.Value;
end