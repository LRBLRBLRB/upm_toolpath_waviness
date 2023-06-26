% The file is aimed to deal with the data measured by the Zygo white light interferometer, 
% and regenerate the measured surface

close all;
clear;
clc;
addpath(genpath('funcs'));
if ~(exist('workspaceDir','var'))
    workspaceDir = fullfile('..','workspace','20230510','MEASURE');
    surfParamUnit = 'mm';
    unit = 'um';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    unitList = {'m','mm','um','nm'};
    presUnit = find(strcmp(unitList,surfParamUnit),1);
    aimUnit = find(strcmp(unitList,unit),1);
end

% surface data import
if ~(exist('surfFunc','var'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A = tand(20)/(2*2000);
    c = 0.69/(1000^(aimUnit - presUnit));
    syms x y;
    surfFunc = @(c,x,y) c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-500,500;-500,500];
    surfDomain = 1.2*surfDomain;
    rMax = max(surfDomain(1,2),surfDomain(2,2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % sampling density
    spar = 1000;
    conR = linspace(0,rMax,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = surfFunc(c,surfMesh(:,:,1),surfMesh(:,:,2));
end

% optim selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FitType = 'trans'; % trans trans-rot rigid rigid+c

surfFitOpt.Algorithm = 'trust-region-reflective';
surfFitOpt.Display = 'iter-detailed';
surfFitOpt.UseParallel = true;
surfFitOpt.MaxIterations = 1000;
surfFitOpt.MaxFunctionEvaluations = 5000;
surfFitOpt.FunctionTolerance = 1e-12;
surfFitOpt.StepTolerance = 1e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3D surface data load
[surfFileName,surfDirName] = uigetfile({ ...
    '*.mat','MATLAB data files(*.mat)'; ...
    '*.datx','data-files(*.datx)'; ...
    '*.xyz','xyz-files(*.xyz)'; ...
    '*.txt','text files(*.txt)'; ...
    '*.*','all files(*.*)'}, ...
    'Select one zygo data file', ...
    fullfile(workspaceDir), ...
    'MultiSelect','off');

surfFilePath = fullfile(surfDirName,surfFileName);
[~,~,surfFileType] = fileparts(surfFilePath);

% extract the original data
switch surfFileType
    case '.mat'
        % case for the mat data exported from datx files
        surfData0 = zygo_datx(surfFilePath,unit);
    case '.datx'
        % case for the .datx data
        return;
    case '.xyz'
        % case for the .xyz data
%         surfData00 = loadXYZ(surfFilePath,1000,1000,8.70983e-07,'um','um');
%         surfData0(:,:,1) = surfData00.X;
%         surfData0(:,:,2) = surfData00.Y;
%         surfData0(:,:,3) = surfData00.Z;
        surfData0 = zygo_xyz(surfFilePath,unit);
    case '.txt'
        % case for the txt data files
        surfData0 = load(surfFilePath);
    otherwise
        return;
end

%% plot the original result
figure('Name','1 Zygo original result');
figPos = get(gcf,'Position');
set(gcf,'Position',[figPos(1)-0.5*figPos(3),figPos(2),2*figPos(3),figPos(4)]);
tiledlayout(1,2);
nexttile;
s = surf( ...
    surfData0(:,:,1),surfData0(:,:,2),surfData0(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
% surf( ...
%     surfData1(:,:,1),surfData1(:,:,2),surfData1(:,:,3), ...
%     'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
nexttile;
waitBar = waitbar(0,'Figure Plotting ...','Name','CNC Data Plot', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(waitBar,'canceling',0);
videoSpeed = 50;
for ii = 1:size(surfData0,1)/videoSpeed
    for jj = 1:size(surfData0,2)/videoSpeed
        % Check for clicked Cancel button
        if getappdata(waitBar,'canceling')
            break;
        end
        displayData = ((ii - 1)*size(surfData0,1)/videoSpeed + jj) ...
            /(size(surfData0,1)*size(surfData0,2)/videoSpeed/videoSpeed); % Calculate percentage
        waitbar(displayData,waitBar,['Figure Plotting ... ', ...
            num2str(roundn(displayData*100,-2),'%.2f'),'%']); % Progress bar dynamic display
        plot3(surfData0(videoSpeed*ii,videoSpeed*jj,1), ...
            surfData0(videoSpeed*ii,videoSpeed*jj,2), ...
            surfData0(videoSpeed*ii,videoSpeed*jj,3),'.', ...
            'Color',[0,0.4470,0.7410],'MarkerSize',0.5); hold on;
    end
end

set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
delete(waitBar);   

questOpt.Interpreter = 'tex';
questOpt.Default = 'OK & Continue';
msgfig = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s}', ...
    'Surface was loaded successfully!\n'],textFontSize,textFontType), ...
    sprintf('The parameters are listed below:'), ...
    sprintf('1. Unit: %s',unit), ...
    sprintf('2. Surface fitting type: %s',FitType), ...
    sprintf('3. Surface fitting algorithm: %s',surfFitOpt.Algorithm), ...
    sprintf('   Max. iterations: %i',surfFitOpt.MaxIterations), ...
    sprintf('   Max. function evaluations: %i',surfFitOpt.MaxFunctionEvaluations), ...
    sprintf('   Function tolerance: %d',surfFitOpt.FunctionTolerance), ...
    sprintf('   Step x tolerance: %d',surfFitOpt.StepTolerance), ...
    'Ready to continue?'}, ...
    'Surface Generation','OK & Continue','Cancel & quit',questOpt);
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    return;
end

%% point cloud registration
% % rough registration
% tmpZ = 0 - min(surfData0(:,:,3),[],'all');
% surfData0(:,:,3) = surfData0(:,:,3) + tmpZ;
% surfDataPt0 = reshape(surfData0,[],3);
% surfMeshPt = reshape(surfMesh,[],3);
% 
% % using zyq's icp function
% surfDataPC0 = pointCloud(surfDataPt);
% surfMeshPC0 = pointCloud(surfMeshPt);
% gridStep = (surfDataPC0.XLimits(2) - surfDataPC0.XLimits(1))/25; % down sampling
% surfDataPC = pcdownsample(surfDataPC0,'gridAverage',gridStep);
% surfMeshPC = pcdownsample(surfMeshPC0,'gridAverage',gridStep);
% [surfRot,surfTran] = icp_yq(surfDataPC.Location,surfMeshPC.Location);
% surfDataPt = (surfRot*(surfDataPt.') + surfTran).';
% 
% 
% using computer vision toolbox
% % out liers removing
% surfDataPt = rmoutliers(surfDataPt,'median');
% % low-pass filter
% % surfDataPt
% 
% figure; plot3(surfDataPt(:,1),surfDataPt(:,2),surfDataPt(:,3),'.', ...
%     'Color',[0.8500 0.3250 0.0980],'MarkerSize',1);
% 
% surfDataPC0 = pointCloud(surfDataPt);
% surfMeshPC = pointCloud(surfMeshPt);
% 
% 
% gridStep1 = (surfDataPC0.XLimits(2) - surfDataPC0.XLimits(1))/25; % down sampling
% [surfHomo1,surfDataPC1,rmse1] = pcregisterndt(surfDataPC0,surfMeshPC,gridStep1);
% gridStep2 = (surfDataPC1.XLimits(2) - surfDataPC1.XLimits(1))/200; % down sampling
% surfDataPC2 = pcdownsample(surfDataPC1,'gridAverage',gridStep2);
% surfMeshPC = pcdownsample(surfMeshPC,'gridAverage',gridStep2);
% [surfHomo2,surfDataPC3,rmse2] = pcregistericp(surfDataPC2,surfMeshPC);
% surfDataPt = reshape(surfData,[],3);
% surfDataPt = (surfDataPt*surfHomo1.Rotation + surfHomo1.Translation) ...
%     *surfHomo2.Rotation + surfHomo2.Translation;
% 
% % manual registration
% surfRot = eul2rotm([0,pi/900,0],'ZYX');
% surfTran = [0;0;0];
% surfDataPt1 = (surfRot*(surfDataPt0') + surfTran)';
% 
% % get the registered measured surface
% surfData1 = reshape(surfDataPt1,size(surfData0,1),size(surfData0,2),3);
% 
% figure('Name','2 measured surface and designed surface');
% surf(surfData1(:,:,1),surfData1(:,:,2),surfData1(:,:,3), ...
%     'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','FaceAlpha',0.1);
% hold on;
% surf(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
%     'FaceColor',[0 0.4470 0.7410],'EdgeColor','none','FaceAlpha',0.1);
% xlabel(['x (',unit,')']);
% ylabel(['y (',unit,')']);
% zlabel(['z (',unit,')']);

%% separate the flat region from the surface, and remove the tilt
% % surf data convertion to image data
% surfSize = size(surfData0);
% surfImg = zeros(surfSize(2),surfSize(1));
% surfImgColorMap = colormap(parula(256));
% surfZ = reshape(surfData0(:,:,3),[],1);
% maxZ = max(surfZ,[],'all');
% minZ = min(surfZ,[],'all');
% surfZInv = linspace(minZ,maxZ,256 + 1);
% surfZsort = sort(surfZ,'descend');
% 
% for col = 1:surfSize(1)
%     for row = 1:surfSize(2)
%         if ~isnan(surfZ(row + (col - 1)*surfSize(1)))
%             surfImg(row,col) = find(surfZ(row + (col - 1)*surfSize(1)) <= surfZInv,1);
%         end
%     end
% end

% to get the flat region and fit it to the z = 0 plane
surfBW = imbinarize(rescale(surfData0(:,:,3),0,1),0.99);
fig2 = figure('Name','2 surface located','WindowState','maximized');
tiled2 = tiledlayout(2,2,'TileSpacing','compact');
ax21 = nexttile(1,[2,1]);
imshow(surfBW,'Parent',ax21);
axis(ax21,'xy');
axis(ax21,'on');
xlabel(ax21,['x (',unit,')']);
ylabel(ax21,['y (',unit,')']);
[fitCenterGp,fitRadiusGp] = imfindcircles(surfBW,[100,1000], ...
    'ObjectPolarity','dark'); % use circular Hough transform to find the largest circle
if ~isempty(fitRadiusGp)
    [fitRadius,fitInd] = max(fitRadiusGp);
    fitCenter = fitCenterGp(fitInd,:);
    viscircles(fitCenter,fitRadius,'Color',[0,0.4450,0.7410]);
end
questOpt.Interpreter = 'tex';
questOpt.Default = 'Re-select';
msgfig1 = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
    'Auto surface selection finished successfully!'],textFontSize,textFontType), ...
    'Re-select or not?'}, ...
    'Concentric tool path Generation','Re-select','Cancel & continue','OK',questOpt);
msgfig2 = 'Re-select';
ifFitPlane = true;
while strcmp(msgfig2,'Re-select')
    if strcmp(msgfig1,'Cancel & continue')
        planeMesh1 = surfData0;
        planeData1 = reshape(planeMesh1,[],3);
        ifFitPlane = false;
        break;
    elseif strcmp(msgfig1,'Re-select')
        msgfig1 = 'OK';
        hold(ax21,'off');
        % imshow(surfImg,surfImgColorMap);
        [fitCenter,fitRadius] = drawCircle(ax21,surfData0, ...
            'textFontSize',textFontSize,'textFontType',textFontType,'unit',unit);
    end

    % to get the indices of the plane
    % isSurf = inROI(roi,surfData0(:,:,1),surfData0(:,:,2));
    isPlane = (surfData0(:,:,1) - fitCenter(1)).^2 + (surfData0(:,:,2) - fitCenter(1)).^2 > fitRadius.^2;
    % plane extraction
    planeMesh1 = surfData0;
    planeMesh1(~(repmat(isPlane,[1,1,3]))) = nan;
    planeData1 = reshape(planeMesh1,[],3);
    ax22 = nexttile(tiled2,2);
    hold(ax22,'off');
    surf(ax22,planeMesh1(:,:,1),planeMesh1(:,:,2),planeMesh1(:,:,3), ...
        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
    colormap(parula(256));
    colorbar(ax22,'eastoutside');
    clim([min(planeMesh1(:,:,3),[],'all'),max(planeMesh1(:,:,3),[],'all')]);
    xlabel(ax22,['x (',unit,')']);
    ylabel(ax22,['y (',unit,')']);
    zlabel(ax22,['z (',unit,')']);
    % plot3(planeData(:,1),planeData(:,2),planeData(:,3),'.','Color',[0,0.4450,0.7410]);
    % pause();
    msgfig2 = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
        'Surface selection finished successfully!'],textFontSize,textFontType), ...
        'Re-selet or not?'}, ...
        'Concentric tool path Generation','Re-select','OK','Cancel & continue',questOpt);
    if strcmp(msgfig2,'Cancel & continue')
        ifFitPlane = false;
        break;
    end
end

% to fit the plane and to remove the tilt
% a.*x + b.*y + c.*z + d = 0;
if ifFitPlane
    [planeNorm,~] = fitPlane(planeData1,'normalization',true);
    planeRot = vecRot(planeNorm,[0;0;1]);
else
    planeRot = eye(3);
end

planeData2 = (planeRot*planeData1')'; % scatters projected to the xOy plane
meshSize = size(planeMesh1);
planeMesh2 = reshape(planeData2,meshSize);

ax23 = nexttile(tiled2,4);
surf(ax23,planeMesh2(:,:,1),planeMesh2(:,:,2),planeMesh2(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
colormap(parula(256));
colorbar(ax23,'eastoutside');
clim(ax23,[min(planeMesh2(:,:,3),[],'all'),max(planeMesh2(:,:,3),[],'all')]);
xlabel(ax23,['x (',unit,')']);
ylabel(ax23,['y (',unit,')']);
zlabel(ax23,['z (',unit,')']);

%% extract the theoretical surface
% z = A.*((x - x0).^2 + (y - y0).^2)./2 + C + z0;

% surface extraction (the same as plane extraction)
fig3 = figure('Name','3 Surface Extraction & fitting','WindowState','maximized');
tiled3 = tiledlayout(2,2,'TileSpacing','compact');
ax31 = nexttile(tiled3,1,[2,1]);

msgfig2 = 'Re-select';
while strcmp(msgfig2,'Re-select')
    % imshow(surfImg,surfImgColorMap);
    [fitCenter,fitRadius] = drawCircle(ax31,surfData0, ...
        'textFontSize',textFontSize,'textFontType',textFontType);

    % to get the indices of the plane
    % isSurf = inROI(roi,surfData0(:,:,1),surfData0(:,:,2));
    isSurface = (surfData0(:,:,1) - fitCenter(1)).^2 + (surfData0(:,:,2) - fitCenter(1)).^2 < fitRadius.^2;
    % surface extraction
    surfMesh1 = surfData0;
    surfMesh1(~(repmat(isSurface,[1,1,3]))) = nan;
    surfData1 = reshape(surfMesh1,[],3);

    ax32 = nexttile(tiled3,2);
    surf(ax32,surfMesh1(:,:,1),surfMesh1(:,:,2),surfMesh1(:,:,3), ...
        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
    colormap(ax32,turbo(256));
    colorbar(ax32,'eastoutside');
    clim(ax23,[min(surfData0(:,:,3),[],'all'),max(surfData0(:,:,3),[],'all')]);
    set(ax32,'FontSize',textFontSize,'FontName',textFontType);
    xlabel(ax32,['x (',unit,')']);
    ylabel(ax32,['y (',unit,')']);
    zlabel(ax32,['z (',unit,')']);
    % plot3(planeData(:,1),planeData(:,2),planeData(:,3),'.','Color',[0,0.4450,0.7410]);
    % pause();
    msgfig2 = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
        'Surface selection finished successfully!'],textFontSize,textFontType), ...
        'Re-selet or not?'}, ...
        'Concentric tool path Generation','Re-select','OK',questOpt);
    hold(ax31,'off');
    hold(ax32,'off');
end

%% least-square nonlinear fitting of the surface
% coeffMat = [surfDataPt1(:,1).^2 + surfDataPt1(:,2).^2,ones(size(surfDataPt1,1),1)];
% reMat = surfDataPt1(:,3);
% fitResult = (coeffMat'*coeffMat)\(coeffMat'*reMat);

% lsqcurvefit (x means the offset of the surface)
surfData1(:,3) = surfData1(:,3) - min(surfData1(:,3));
surfDataFit = surfData1(all(~isnan(surfData1),2),:); % remove all the nans in surfData1
surfFitOpts = optimoptions('lsqnonlin', ...
    'Algorithm',surfFitOpt.Algorithm, ...
    'Display',surfFitOpt.Display,'UseParallel',surfFitOpt.UseParallel, ...
    'MaxIterations',surfFitOpt.MaxIterations, ...
    'MaxFunctionEvaluations',surfFitOpt.MaxFunctionEvaluations, ...
    'FunctionTolerance',surfFitOpt.FunctionTolerance, ...
    'StepTolerance',surfFitOpt.StepTolerance);
switch FitType
    case 'trans'
        f = @(x,data) surfFunc(c,data(:,1) - x(1),data(:,2) - x(2)) + x(3);
        surfFitOpts = optimoptions('lsqcurvefit',surfFitOpts);
        surfFitResult = lsqcurvefit(f,[0;0;min(surfDataFit(:,3))], ...
            surfDataFit(:,1:2),surfDataFit(:,3),[],[],surfFitOpts);
        surfData2 = surfData1 - meshgrid(surfFitResult,1:size(surfData1,1));
        surfTheo(:,1:2) = surfData2(:,1:2); % fit surface
        surfTheo(:,3) = f([0;0;0],surfTheo(:,1:2));
    case 'trans-rot'
    case 'rigid'
        f = @(H) surfHomogeneous(surfFunc,surfDataFit',c,H(1:3),H(4:5));
        surfFitOpts = optimoptions('lsqnonlin',surfFitOpts);
        surfFitResult = lsqnonlin(f,[0;0;0;0;0], ...
            [-100;-100;-100;-5;-5],[100;100;100;5;5],surfFitOpts);
        surfData2 = (rotx(surfFitResult(4))*roty(surfFitResult(5))*(surfData1') ...
            + ndgrid(surfFitResult(1:3),1:size(surfData1,1)))';
        surfTheo(:,1:2) = surfData2(:,1:2); % fit surface
        surfTheo(:,3) = surfFunc(c,surfTheo(:,1),surfTheo(:,2));
    case 'rigid+c'
        f = @(H) surfHomogeneous(surfFunc,surfDataFit',H(6),H(1:3),H(4:5));
        surfFitOpts = optimoptions('lsqnonlin',surfFitOpts);
        surfFitResult = lsqnonlin(f,[0;0;-1*min(surfDataFit(:,3));0;0;c],[],[],surfFitOpts);
        surfData2 = (rotx(surfFitResult(4))*roty(surfFitResult(5))*(surfData1') ...
            + ndgrid(surfFitResult(1:3),1:size(surfData1,1)))';
        surfTheo(:,1:2) = surfData2(:,1:2); % fit surface
        surfTheo(:,3) = surfFunc(surfFitResult(6),surfTheo(:,1),surfTheo(:,2));
end

% plot the result
ax33 = nexttile(tiled3,4);
surfMesh2 = reshape(surfData2,meshSize); % surface mesh data whose offset has been removed
surf(ax33,surfMesh2(:,:,1),surfMesh2(:,:,2),surfMesh2(:,:,3), ...
        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
colormap(ax33,parula(256));
colorbar(ax33,'eastoutside');
clim(ax33,[min(surfMesh2(:,:,3),[],'all'),max(surfMesh2(:,:,3),[],'all')]);
hold(ax33,'on');
surfTheoMesh = reshape(surfTheo,meshSize);
surf(ax33,surfTheoMesh(:,:,1),surfTheoMesh(:,:,2),surfTheoMesh(:,:,3), ...
        'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha',0.3,'LineStyle','none');
set(ax33,'FontSize',textFontSize,'FontName',textFontType);
xlabel(ax33,['x (',unit,')']);
ylabel(ax33,['y (',unit,')']);
zlabel(ax33,['z (',unit,')']);

%% surface error calculation
deltaZ(:,:,1:2) = surfMesh2(:,:,1:2);
deltaZ(:,:,3) = surfMesh2(:,:,3) - surfTheoMesh(:,:,3);
% deltaZFit = deltaZ(~isnan(deltaZ));
% Sa = mean(abs(deltaZFit));
% Sz = max(deltaZFit,[],'all') - min(deltaZFit,[],'all');
% Sq = sqrt(mean(deltaZFit.^2)); % rms

figure('Name','4 - 3D surface error');
tiledlayout(2,1);
nexttile;
surf(deltaZ(:,:,1),deltaZ(:,:,2),deltaZ(:,:,3),'EdgeColor','none');
% axis equal;
hold('on');
colormap(turbo(256));
colorbar('eastoutside');
clim([min(deltaZ(:,:,3),[],'all'),max(deltaZ(:,:,3),[],'all')]);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['\Deltaz (',unit,')']);

% cut
[lineAng,lineData] = viewError(deltaZ,textFontSize + 2,textFontType,unit);
lineRot = rotz(lineAng);
minX = min(deltaZ(:,:,1),[],'all');
maxX = max(deltaZ(:,:,1),[],'all');
minZ = min(deltaZ(:,:,3),[],'all');
maxZ = max(deltaZ(:,:,3),[],'all');
linePts = [1.2*minX,1.2*maxX,1.2*maxX,1.2*minX;0,0,0,0; ...
    1.2*minZ,1.2*minZ,1.2*maxZ,1.2*maxZ];
linePts = lineRot*linePts;
fill3(linePts(1,:),linePts(2,:),linePts(3,:), ...
    [0.6350 0.0780 0.1840],'EdgeColor','none','FaceAlpha',0.3);
line([1.2*minX;1.2*maxX],[0;0],[1.2*maxZ;1.2*maxZ], ...
    'Color',[0.6350 0.0780 0.1840],'LineWidth',3, ...
    'Marker','.','MarkerSize',18);

nexttile;
plot(lineData(:,1),lineData(:,2));
set(gca,'FontSize',textFontSize,'FontName',textFontType);
grid on;
xlabel(['r (',unit,')']);
ylabel(['\Deltaz (',unit,')']);

% s7_extract;

%%
% rmpath(genpath('funcs'));

% function of surface homogeneous transformation parameters fitting
function deltaZ = surfHomogeneous(f,data,c,x,ang)
data1 = rotx(ang(1))*roty(ang(2))*data + ndgrid(x,1:size(data,2));
deltaZ = data1(3,:) - f(c,data1(1,:),data1(2,:));
end