function surfacePlot(app,method)

app.surfMesh = zeros(app.surfPlotSpar,app.surfPlotSpar,3);
cla(app.SurfaceDataAxes,'reset');
% process the surface function or point cloud
switch app.surfType
    case app.ImportCell
        % sampling
        surfaceMeshGen(app);

        switch method
            case '2D'
                % plot the importing results
                surf(app.SurfaceDataAxes, ...
                    app.surfMesh(:,:,1),app.surfMesh(:,:,2),app.surfMesh(:,:,3), ...
                    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                view(app.SurfaceDataAxes,2);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            case '3D'
                % plot the importing results
                surf(app.SurfaceDataAxes, ...
                    app.surfMesh(:,:,1),app.surfMesh(:,:,2),app.surfMesh(:,:,3), ...
                    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                zlabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
        end
    case app.Geometry2DCell % 2D Geometry
        app.rMax = app.surfDomain(1,2);

        switch method
            case '2D'
                % sampling
                conx = linspace(app.surfDomain(1,1),app.surfDomain(1,2),app.surfPlotSpar);
                conz = app.surfFuncs(conx,zeros(1,length(conx)));

                % plot the importing results
                plot(app.SurfaceDataAxes,conx,conz);
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['r (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            case '3D'
                % sampling
                surfaceMeshGen(app);
        
                % plot the importing results
                surf(app.SurfaceDataAxes, ...
                    app.surfMesh(:,:,1),app.surfMesh(:,:,2),app.surfMesh(:,:,3), ...
                    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                zlabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
        end
    case app.Geometry3DCell % 3D Geometry
        % partial differential function
        syms x y;
        app.surfFx = matlabFunction(diff(app.surfFuncs,x));
        app.surfFy = matlabFunction(diff(app.surfFuncs,y));

        % sampling
        surfaceMeshGen(app);

        switch method
            case '2D'
                % plot the importing results
                surf(app.SurfaceDataAxes, ...
                    app.surfMesh(:,:,1),app.surfMesh(:,:,2),app.surfMesh(:,:,3), ...
                    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                view(app.SurfaceDataAxes,2);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
            case '3D'
                % plot the importing results
                surf(app.SurfaceDataAxes, ...
                    app.surfMesh(:,:,1),app.surfMesh(:,:,2),app.surfMesh(:,:,3), ...
                    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
                hold(app.SurfaceDataAxes,'on');
                grid(app.SurfaceDataAxes,'on');
                set(app.SurfaceDataAxes,'FontSize',app.fontSize, ...
                    'FontName',app.fontName);
                xlabel(app.SurfaceDataAxes,['x (',app.unit,')']);
                ylabel(app.SurfaceDataAxes,['y (',app.unit,')']);
                zlabel(app.SurfaceDataAxes,['z (',app.unit,')']);
                title(app.SurfaceDataAxes,'Surface Original Data','FontSize',14);
        end
end

end