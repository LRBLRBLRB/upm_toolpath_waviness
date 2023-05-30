% The file is aimed to deal with the data measured by the Zygo white light interferometer, 
% and regenerate the measured surface

% close all;
clear;
clc;
addpath(genpath('funcs'));
if ~(exist('workspaceDir','var'))
    workspaceDir = fullfile('..','workspace','20230510','MEASURE');
    surfUnit = 'um';
    unit = 'um';
    textFontSize = 12;
    textFontType = 'Times New Roman';
    unitList = {'m','mm','um','nm'};
    presUnit = find(strcmp(unitList,surfUnit),1);
    aimUnit = find(strcmp(unitList,unit),1);
end

%% surface data import
if ~(exist('surfFunc','var'))
    % A = tand(20)/(2*2000);
    c = 0.18/(1000^(aimUnit - presUnit));
    syms x y;
    surfSym = c.*(x.^2 + y.^2)./(1 + sqrt(1 - c.^2.*(x.^2 + y.^2)));
    surfFunc = matlabFunction(surfSym);
    surfFx = diff(surfFunc,x);
    surfFy = diff(surfFunc,y);
    surfDomain = [-1000,1000;-1000,1000];
    surfDomain = 1.05*surfDomain;
    rMax = max(surfDomain(1,2),surfDomain(2,2));
    % sampling density
    spar = 1000;
    conR = linspace(0,rMax,spar); % concentric radius vector
    conTheta = linspace(0,2*pi,spar);
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = surfFunc(surfMesh(:,:,1),surfMesh(:,:,2));
end

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
    'Re-selet or not?'}, ...
    'Concentric tool path Generation','Re-select','OK',questOpt);
msgfig2 = 'Re-select';
while strcmp(msgfig2,'Re-select')
    if strcmp(msgfig1,'Re-select')
        hold(ax21,'off');
        % imshow(surfImg,surfImgColorMap);
        contourf(ax21,surfData0(:,:,1),surfData0(:,:,2),surfData0(:,:,3),16,'--','LineWidth',0.1); hold on;
        colormap(parula(256));
        colorbar(ax21,'southoutside');
        axis(ax21,'equal');
        set(ax21,'FontSize',textFontSize,'FontName',textFontType);
        xLim = get(ax21,'XLim');
        yLim = get(ax21,'YLim');
        set(ax21,'XLim',sqrt(2)*xLim,'YLim',sqrt(2)*yLim);
        xlabel(ax21,['x (',unit,')']);
        ylabel(ax21,['y (',unit,')']);
        fprintf(['Please draw a circle to include all the surface: \n' ...
            '(Double-click to return to the main procedure)\n\n']);
        roi = drawcircle(ax21,'Color',[0.8500 0.3250 0.0980]);
        wait(roi);
        fitCenter = roi.Center;
        fitRadius = roi.Radius;
    end

    % to get the indices of the plane
    % isSurf = inROI(roi,surfData0(:,:,1),surfData0(:,:,2));
    isPlane = (surfData0(:,:,1) - fitCenter(1)).^2 + (surfData0(:,:,2) - fitCenter(1)).^2 > fitRadius.^2;
    % plane extraction
    planeMesh1 = surfData0;
    planeMesh1(~(repmat(isPlane,[1,1,3]))) = nan;
    planeData1 = reshape(planeMesh1,[],3);
    ax22 = nexttile(2);
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
        'Concentric tool path Generation','Re-select','OK',questOpt);
end

% to fit the plane and to remove the tilt
% a.*x + b.*y + c.*z + d = 0;
[planeNorm,~] = fitPlane(planeData1,'normalization',true);
planeRot = vecRot(planeNorm,[0;0;1]);

planeData2 = (planeRot*planeData1')'; % scatters projected to the xOy plane
planeMeshSize = size(planeMesh1);
planeMesh2 = reshape(planeData2,planeMeshSize);

ax23 = nexttile(4);
surf(ax23,planeMesh2(:,:,1),planeMesh2(:,:,2),planeMesh2(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
colormap(parula(256));
colorbar(ax23,'eastoutside');
clim([min(planeMesh2(:,:,3),[],'all'),max(planeMesh2(:,:,3),[],'all')]);
xlabel(ax23,['x (',unit,')']);
ylabel(ax23,['y (',unit,')']);
zlabel(ax23,['z (',unit,')']);

%% fit the theoretical surface
% surface extraction (the same as plane extraction)





surfMesh1 = surfData0;
surfMesh1((repmat(isPlane,[1,1,3]))) = nan;
surfData1 = reshape(surfMesh1,[],3);

% z = A.*((x - x0).^2 + (y - y0).^2)./2 + C + z0;

% least-square result of linear equations
% coeffMat = [surfDataPt1(:,1).^2 + surfDataPt1(:,2).^2,ones(size(surfDataPt1,1),1)];
% reMat = surfDataPt1(:,3);
% fitResult = (coeffMat'*coeffMat)\(coeffMat'*reMat);

% lsqcurvefit
f = @(x,xdata) c.*((xdata(:,1) - x(1)).^2 + (xdata(:,2) - x(2)).^2) ...
    ./(1 + sqrt(1 - c.^2.*((xdata(:,1) - x(1)).^2 + (xdata(:,2) - x(2)).^2))) + x(3);
x = lsqcurvefit(f,[0;0;0],surfData1(:,1:2),surfData1(:,3));

%% surface error calculation
deltaZ = surfData1(:,:,3) - surfFunc(surfData1(:,:,1),surfData1(:,:,2));
Sa = mean(abs(deltaZ));
Sz = max(deltaZ,[],'all') - min(deltaZ,[],'all');
Sq = sqrt(mean(deltaZ.^2));

figure('Name','3 - surface error');
surf(surfData1(:,:,1),surfData1(:,:,2),deltaZ,'EdgeColor','none');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['\Deltaz (',unit,')']);

%%
% rmpath(genpath('funcs'));