function resetSurfaceParams(app)
    clear surfType surfDomain surfPathName ...
        surfPt surfFuncs surfFx surfFy surfMesh;
    app.SurfaceDetailTa.Value = '';
    app.surfPlotSpar = app.surfPlotSparDefault;
    app.surfMesh = zeros(app.surfPlotSpar,app.surfPlotSpar,3);
end

