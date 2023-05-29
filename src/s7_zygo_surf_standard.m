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

% plot the original result
figure('Name','1 Zygo original result');
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
for ii = 1:size(surfData0,1)
    for jj = 1:size(surfData0,2)
        % Check for clicked Cancel button
        if getappdata(waitBar,'canceling')
            break;
        end
        displayData = num2str(roundn(((ii - 1)*size(surfData0,1) + jj) ...
            /(size(surfData0,1)*size(surfData0,2))*100,-2)); % Calculate percentage
        displayStr = ['Figure Plotting ... ',displayData,'%']; % Show Calculate State
        waitbar(((ii - 1)*size(surfData0,1) + jj) ...
            /(size(surfData0,1)*size(surfData0,2)),waitBar,displayStr); % Progress bar dynamic display
        plot3(surfData0(ii,jj,1),surfData0(ii,jj,2),surfData0(ii,jj,3),'.', ...
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

%% rigid transform
% surf data convertion to image data
surfSize = size(surfData0);
surfImgInd = zeros(surfSize(2),surfSize(1));
surfImgColorMap = colormap(parula(512));
surfZ = reshape(surfData0(:,:,3),[],1);
[~,surfZInd] = sort(surfZ,'descend');
for col = 1:surfSize(1)
    for row = 1:surfSize(2)
        surfImgInd(row,col) = surfZInd(row + (col - 1)*surfSize(1));
    end
end
maxZ = max(surfImgInd,[],'all');
minZ = min(surfImgInd,[],'all');
surfZInv = linspace(maxZ,minZ,512 + 1);

surfImgColorInd = surfImgInd;
for ii = 1:surfSize(2)
    for jj = 1:surfSize(1)
        if ~isnan(surfImgInd(ii,jj))
            surfImgColorInd(ii,jj) = find(surfImgInd(ii,jj) >= surfZInv,1);
        end
    end
end

% to get the flat region and fit it to the z = 0 plane
figure('Name','2 surface located');
% hPlot1 = contourf(surfData0(:,:,1),surfData0(:,:,2),surfData0(:,:,3)); hold on;
% view(0,90);
imshow(surfImgColorInd,surfImgColorMap);
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);

% isContinue = checkextractcir();
% while ~isContinue
%     fprintf('Please select the surface data with selecting the center and radius of a circle: \n');
%     disableDefaultInteractivity(gca);
%     k = waitforbuttonpress; % 等待鼠标按下
%     point1 = get(gca,'CurrentPoint'); % 鼠标按下了
%     finalRect = rbbox; %
%     point2 = get(gca,'CurrentPoint'); % 鼠标松开了
%     circCenter = point1(1,1:2); % the center of the drawed circle
%     circEdge = point2(1,1:2); % the edge point of it
%     circRadius = norm(circEdge - circCenter); % the radius of it
%     boxDel = [circCenter(1) - circRadius,circCenter(2) - circRadius,2*circRadius,2*circRadius];
%     hplot = rectangle('Position',boxDel,'Curvature',[1,1],'FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980]); %在原图上显示截取的区域
%     [isContinue] = checkextractcir(hplot);
%     delete(hPlot1);
%     hPlot1 = contourf(surfData0(:,:,1),surfData0(:,:,2),surfData0(:,:,3)); hold on;
% end


% roi = images.roi.Circle('StripeColor',[0.8500 0.3250 0.0980],'SelectedColor',[0.6350 0.0780 0.1840]);
% fprintf(['Please draw a circle to include all the surface: \n' ...
%     '(Double-click to return to the main procedure)\n']);
% roi = drawcircle('Color',[0.8500 0.3250 0.0980]);
% wait(roi);
% isSurf = inROI(roi,surfData0(:,:,1),surfData0(:,:,2));

%% fit the theoretical surface
% z = A.*((x - x0).^2 + (y - y0).^2)./2 + C + z0;
surfDataPt1 = reshape(surfData0,[],3);

% least-square result of linear equations
% coeffMat = [surfDataPt1(:,1).^2 + surfDataPt1(:,2).^2,ones(size(surfDataPt1,1),1)];
% reMat = surfDataPt1(:,3);
% fitResult = (coeffMat'*coeffMat)\(coeffMat'*reMat);

% lsqcurvefit
f = @(x,xdata) x(1)/2*(xdata(:,1).^2 + xdata(:,2).^2) + x(2);
x = lsqcurvefit(f,c,surfDataPt1(:,1:2),surfDataPt1(:,3));

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