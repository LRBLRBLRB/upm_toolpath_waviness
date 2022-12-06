xlim = 1000*input('x range (unit: mm): \n');
if isempty(xlim)
    xlim = 1000*[-5,5];
elseif ~isequal(size(xlim),[1,2])
    error('Invalid x range. A 1*2 array is needed.');
end
ylim = 1000*input('y range (unit: mm): \n');
if isempty(ylim)
    ylim = 1000*[-5,5];
elseif ~isequal(size(ylim),[1,2])
    error('Invalid y range. A 1*2 array is needed.');
end
surfType = 'function';
if strcmp(surfType,'function') % analytic surface
    syms x y;
    surfSym = input("Type the function expression of the surface," + ...
        "with the independent variables of \'x\' and \'y\': \n");
    if isempty(surfSym) % default: biconic surface
        rx = 1000*10; % radius in x direction
        cx = 1/rx; % curvature in x direction
        ry = 1000*10; % radius in y direction
        cy = 1/ry; % curvature in y direction
        kx = -1; % conic in x direction
        ky = -1; % conic in y direction
        Ar = 0; % polynomials (same below)
        Ap = 0;
        Br = 0;
        Bp = 0;
        Cr = 0;
        Cp = 0;
        Dr = 0;
        Dp = 0;
        surfSym = (cx*x^2 + cy*y^2) ...
            /(1 + sqrt(1 - (1 + kx)*cx^2*x^2 - (1 + ky)*cy^2*y^2)) ...
            + Ar*((1 - Ap)*x^2 + (1 + Ap)*y^2)^2 ...
            + Br*((1 - Bp)*x^2 + (1 + Bp)*y^2)^3 ...
            + Cr*((1 - Cp)*x^2 + (1 + Cp)*y^2)^4 ...
            + Dr*((1 - Dp)*x^2 + (1 + Dp)*y^2)^5;
    elseif ~(contains(surfSym,'x') && contains(surfSym,'y'))
        error('invalid function expression: independent variables must be x & y.');
    end

    surfFunc = matlabFunction(surfSym,'vars',[x,y]);
    surfFx = matlabFunction(diff(surfFunc,x),'vars',[x,y]);
    surfFy = matlabFunction(diff(surfFunc,y),'vars',[x,y]);
    
    % surfMesh for drawing surf-mesh pictures
    R = sqrt(xlim(2)^2 + ylim(2)^2);
    reconR = linspace(0,R);
    reconTheta = linspace(0,2*pi);
    [surfMeshPolar(:,:,1),surfMeshPolar(:,:,2)] = meshgrid(reconR,reconTheta); % 这里的meshgrid，x和y的顺序是不是我想要的？
    surfMesh(:,:,1) = surfMeshPolar(:,:,1).*cos(surfMeshPolar(:,:,2));
    surfMesh(:,:,2) = surfMeshPolar(:,:,1).*sin(surfMeshPolar(:,:,2));
    surfMesh(:,:,3) = surfFunc(surfMeshPolar(:,:,1),surfMeshPolar(:,:,2));

    % discretization of the freeform surface (constant arc)
    R = sqrt(xlim(2)^2 + ylim(2)^2);
    r = (toolData.radius/2):(toolData.radius/2):R;
    arcLength = 30;
    maxAngPtDist = 6*pi/180;
    
    surfX = [];
    surfY = [];
    surfNorm = [];
    for ii = 1:length(r)
        conTheta = linspace(0,2*pi, ...
            ceil(2*pi/min(maxAngPtDist,arcLength/r(ii))) + 1);
        conTheta(end) = [];
        sparTheta(ii) = length(conTheta);
        surfX = [surfX,r(ii)*cos(conTheta)];
        surfY = [surfY,r(ii)*sin(conTheta)];
    end
    ptNum = length(surfX);
    surfPt(1:2,:) = [surfX;surfY];
    surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:));
    surfNorm(1,:) = surfFx(surfPt(1,:),surfPt(2,:));
    surfNorm(2,:) = surfFy(surfPt(1,:),surfPt(2,:));
    surfNorm(3,:) = -1*ones(1,ptNum);
    surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1));
    surfDirect = cutDirection(surfPt,[0;0;0]);
    clear surfX surfY;

      
else % point cloud input
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select the surface edge data file', ...
        'workspace\input_data\surface\ellipsoidAray.mat', ...
        'MultiSelect','off');
    surfName = fullfile(dirName,fileName);
    surfData = load(surfName);
    num = size(surfData.xyz,2);
    
    % B-spline interpolation of the freeform surface in polar coordinate
    R = sqrt(xlim(2)^2 + ylim(2)^2);
    reconR = linspace(0,R,ceil(sqrt(num))); % 重建所需的点数值得再考虑
    reconTheta = linspace(0,2*pi,ceil(sqrt(num)));
    [surfMeshPolar(:,:,1),surfMeshPolar(:,:,2)] = meshgrid(reconR,reconTheta); % 这里的meshgrid，x和y的顺序是不是我想要的？
    surfMeshPolar(:,:,3) = griddata(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3), ...
        surfMeshPolar(:,:,1),surfMeshPolar(:,:,2),'linear');
    % surfMesh for drawing surf-mesh pictures
    surfMesh(:,:,1) = surfMeshPolar(:,:,1).*cos(surfMeshPolar(:,:,2));
    surfMesh(:,:,2) = surfMeshPolar(:,:,1).*sin(surfMeshPolar(:,:,2));
    surfMesh(:,:,3) = surfFunc(surfMeshPolar(:,:,1),surfMeshPolar(:,:,2));
    % [surfCpts,U,V] = bSplineSurfCpts(surfData.mesh,3,3);
    u = 0:0.002:1;
    v = 0:0.002:1;
    % surfPts = zeros(length(u),length(v),dim);
    % surfPts = bSplineSurfPts(surfCpts,3,u,U,3,v,V);
    [surfPt,surfBform] = bsplineSurfPts_spapi(surfMeshPolar,3,3,u,v,'paramMethod','centripetal');
    
    figSurfInt = figure('Name','Freeform surface interpolation results');
    % pose = get(gcf,'Position');
    % set(gcf,"Position",[pose(1)-pose(3)/2,pose(2),2*pose(3),pose(4)]);
    plot3(surfData.xyz(:,1),surfData.xyz(:,2),surfData.xyz(:,3), ...
        '.','Color',[0,0.45,0.74]); 
    hold on;
    plot3(reshape(surfCpts(:,:,1),[],1), ...
        reshape(surfCpts(:,:,2),[],1), ...
        reshape(surfCpts(:,:,3),[],1), ...
        'x','Color',[0.85,0.33,0.10]);
    surf(surfPt(:,:,1),surfPt(:,:,2),surfPt(:,:,3), ...
        'FaceColor','interp','EdgeColor','none');
    cb = colorbar('Location','eastoutside');
    set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
    xlabel(['x (',unit,')']);
    ylabel(['y (',unit,')']);
    zlabel(['z (',unit,')']);
    grid on
    legend('Measured Pts','Interpolation Cpts','Interpolation Results', ...
        'Orientation','vertical','Location','best');

    % calculate the derivative and cutting direction of the surface
    surfBform = fnder(surfBform,{1,0});
    % ???????????????????????????????????????????????????????????????????????????????????????????????????

    % discretization of the freeform surface (constant arc)
    r = (toolData.radius/2):(toolData.radius/2):R;
    arcLength = 30;
    maxAngPtDist = 6*pi/180;

    surfX = [];
    surfY = [];
    surfNorm = [];
    for ii = 1:length(r)
        conTheta = linspace(0,2*pi, ...
            ceil(2*pi/min(maxAngPtDist,arcLength/r(ii))) + 1);
        conTheta(end) = [];
        sparTheta(ii) = length(conTheta);
        surfX = [surfX,r(ii)*cos(conTheta)];
        surfY = [surfY,r(ii)*sin(conTheta)];
    end
    ptNum = length(surfX);
    surfPt(1:2,:) = [surfX;surfY];
    surfPt(3,:) = surfFunc(surfPt(1,:),surfPt(2,:));
    surfNorm(1,:) = surfFx(surfPt(1,:),surfPt(2,:));
    surfNorm(2,:) = surfFy(surfPt(1,:),surfPt(2,:));
    surfNorm(3,:) = -1*ones(1,ptNum);
    surfNorm = -1*(surfNorm./vecnorm(surfNorm,2,1));
    surfDirect = cutDirection(surfPt,[0;0;0]);
    clear surfX surfY;
end
