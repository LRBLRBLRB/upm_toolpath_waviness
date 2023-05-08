function surfaceMeshGen(app)
% to generate the mesh points of the surface

% process the surface function or point cloud
switch app.surfType
    case app.ImportCell
        uialert(app.UIFigure,'There has not been methods fot import surface. ', ...
            'Warning','Icon','warning');

    case app.Geometry2DCell % 2D Geometry
        % partial differential function
        syms x y;
        app.surfFx = diff(app.surfFuncs,x);
        app.surfFy = diff(app.surfFuncs,y);

        % sampling
        conR = linspace(0,app.surfDomain(1,2),app.surfPlotSpar); % concentric radius vector
        conTheta = linspace(0,2*pi,app.surfPlotSpar);
        [rMesh,thetaMesh] = meshgrid(conR,conTheta);
        app.surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
        app.surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
        app.surfMesh(:,:,3) = app.surfFuncs(app.surfMesh(:,:,1),app.surfMesh(:,:,2));
        app.rMax = app.surfDomain(1,2);
    case app.Geometry3DCell % 3D Geometry
        % partial differential function
        syms x y;
        app.surfFx = diff(app.surfFuncs,x);
        app.surfFy = diff(app.surfFuncs,y);

        % sampling
        if size(app.surfDomain,1) == 1
            conR = linspace(0,app.surfDomain(1,2),app.surfPlotSpar); % concentric radius vector
            conTheta = linspace(0,2*pi,app.surfPlotSpar);
            [rMesh,thetaMesh] = meshgrid(conR,conTheta);
            app.surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
            app.surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
            app.rMax = app.surfDomain(1,2);
        else
            xx = linspace(app.surfDomain(1,1),app.surfDomain(1,2),app.surfPlotSpar);
            yy = linspace(app.surfDomain(2,1),app.surfDomain(2,2),app.surfPlotSpar);
            [app.surfMesh(:,:,1),app.surfMesh(:,:,2)] = meshgrid(xx,yy);
            app.rMax = norm(app.surfDomain(:,2) - app.surfDomain(:,1));
        end
        app.surfMesh(:,:,3) = app.surfFuncs(app.surfMesh(:,:,1),app.surfMesh(:,:,2));
end
