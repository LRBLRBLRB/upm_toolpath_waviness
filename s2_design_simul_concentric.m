% simulation of the designed tool path of a concentric surface
% Step one: tool and surface data import
% Step two: tool path calculation
% Step three: residual height calculation and correction
% Step four: simulation of the machining surface
close all;
clear;
clc;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType;
unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

msgOpts.Default = 'Cancel and quit';
msgOpts.Interpreter = 'tex';
% msgOpts.modal = 'non-modal';
profile on
parObj = gcp;

%% concentric surface generation / import
% [fileName,dirName] = uigetfile({ ...
%     '*.mat','MAT-files(*.mat)'; ...
%     '*,*','all files(*.*)'}, ...
%     'Select one tool edge data file', ...
%     'output_data\tool\tooltheo.mat', ...
%     'MultiSelect','off');
% toolName = fullfile(dirName,fileName);
toolName = 'output_data\tool\toolTheo_3D.mat';
toolData = load(toolName);

default = false;
if default
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','all files(*.*)'}, ...
        'Select the surface edge data file', ...
        'input_data\surface\ellipsoidAray.mat', ...
        'MultiSelect','off');
    surfName = fullfile(dirName,fileName);
    load(surfName);
% % % % % % % % % % % % % % % % % % % % % % % % %     sparTheta
% % % % % % % % % % % % % % % % % % % % % % % % %     surfMesh
% % % % % % % % % % % % % % % % % % % % % % % % %     surfNorm
else % ellipsoid
    R = 10/2*1000;
    A = 3.5/2;
    B = 4/2;
    C = 5/2;
    % sampling density
    r = [0,R/4]; % concentric radius range
    sparTheta = 101;
    surfCenter = [0,0,sqrt(C^2*(R.^2-r(2).^2))]; % concentric circle center
    conR = (toolData.radius/8):(toolData.radius/2):r(2); % concentric radius vector
    densR = length(conR);
    conTheta = linspace(0,2*pi,sparTheta);
    
    [rMesh,thetaMesh] = meshgrid(conR,conTheta);
    surfMesh(:,:,1) = A*rMesh.*cos(thetaMesh);
    surfMesh(:,:,2) = B*rMesh.*sin(thetaMesh);
    surfMesh(:,:,3) = sqrt(C^2*(R.^2-rMesh.^2));
    %  y = reshape(surfMesh(:,:,2),[],1);
    % z = reshape(surfMesh(:,:,3),[],1);
    % surfXYZ = [x,y,z];
    
    % calculate the normal vector of the analytic surface
    % syms X Y z (X,Y);
    % z (X,Y) = sqrt(C^2*(R^2 - X.^2/A^2 - Y.^2/B^2));
    % ZDX = diff(Z,X);
    % ZDY = diff(Z,Y);
    % surfNorm(:,1) = eval(subs(ZDX,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,2) = eval(subs(ZDY,{X,Y},{surfXYZ(:,1),surfXYZ(:,2)}));
    % surfNorm(:,3) = -ones(densTheta*densR,1);
    [surfNorm(:,:,1),surfNorm(:,:,2),surfNorm(:,:,3)] = surfnorm( ...
        surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3));
%     save('input_data/surface/ellipsoidAray.mat', ...
%         "surfMesh","surfNorm","surfCenter");
end

%% plot the freeform surface
figure('Name','original xyz scatters of the surface (sparsely)');
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
tiledlayout(1,2);
nexttile;
rSpar = 1;
theSpar = 1;
plot3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    '.','Color',[0,0.45,0.74]);
hold on;
legend('Original Points','Location','northeast');
quiver3( ...
    surfMesh(1:theSpar:end,1:rSpar:end,1), ...
    surfMesh(1:theSpar:end,1:rSpar:end,2), ...
    surfMesh(1:theSpar:end,1:rSpar:end,3), ...
    surfNorm(1:theSpar:end,1:rSpar:end,1), ...
    surfNorm(1:theSpar:end,1:rSpar:end,2), ...
    surfNorm(1:theSpar:end,1:rSpar:end,3), ...
    'AutoScale','on','Color',[0.85,0.33,0.10],'DisplayName','Normal Vectors');
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
axis equal; grid on;
% title({'Radially & circunferentially sparse','by 50 and 2 times, respectively'}, ...
%     'FontSize',textFontSize,'FontName',textFontType);
nexttile;
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
hold on; axis equal;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% cb2 = colorbar;
msgfig = questdlg({'Surface was generated successfully!', ...
    'Ready to continue?'}, ...
    'Surface Generation','OK & continue','Cancel & quit','OK & continue');
if strcmp(msgfig,'Cancel & quit') || isempty(msgfig)
    msgbox('Exit for the program','Exit','error','modal');
    uiwait(msgbox);
    return;
end
surfPt = transpose(reshape(surfMesh,[],3));
surfNorm = transpose(reshape(surfNorm,[],3));
surfDirect = cutDirection(surfPt,surfCenter);

%% Calculation of Tool Path & Spindle Drection
ptNum = size(surfPt,2);
toolQuat = zeros(ptNum,4);
toolVec = zeros(3,ptNum);
toolPathPt = zeros(3,ptNum);
toolCutDirect = zeros(3,ptNum);
toolNormDirect = zeros(3,ptNum);
toolContactU = zeros(1,ptNum);
isCollision = false(1,ptNum);
tic
parfor ii = 1:ptNum
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),surfDirect(:,ii),[0;0;1]);
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
end
tToolpath = toc;
fprintf('The time spent in the tool path generating process is %fs.\n',tToolpath);

% still need procedures to deal with the invalid points, which will cause
% interference between the tool edge and the designed surface. 

figure('Name','tool center position & tool normal vector');
plotSpar = 1;
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','Color',[0.6350,0.0780,0.1840]);
hold on;
quiver3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    toolCutDirect(1,1:plotSpar:end), ...
    toolCutDirect(2,1:plotSpar:end), ...
    toolCutDirect(3,1:plotSpar:end), ...
    'AutoScale','on','Color',[0,0.4470,0.7410]);
quiver3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    toolNormDirect(1,1:plotSpar:end), ...
    toolNormDirect(2,1:plotSpar:end), ...
    toolNormDirect(3,1:plotSpar:end), ...
    'AutoScale','on','Color',[0.85,0.33,0.10]);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',0.6,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = 'Height (mm)';
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
legend('tool center point','tool cutting direction', ...
    'tool spindle direction','','Location','northeast');
msgfig = msgbox('Tool Path was calculated successfully!', ...
    'Tool Path Simulation','warn','non-modal');
uiwait(msgfig);

%% Calculation of Residual Height & Cutting Surface
toolSp = toolData.toolBform;
toolRadius = toolData.radius;
resNum = ptNum - sparTheta;
res = zeros(1,2*resNum);
peakPt = zeros(3,2*resNum);
uLim = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of each toolpath
% uLimTmp = [zeros(1,ptNum);ones(1,ptNum)]; % the interval of the projective toolpath
angle = atan2(toolPathPt(2,:),toolPathPt(1,:));

tic
% outer side of each point in the tool path
parfor ii = 1:resNum
    % 如果是沿同一个极径的，就可以直接不用投影；否则还是需要这样子找
    nLoop = floor((ii - 1)/sparTheta) + 1;
    angleN = angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1));
    % ind2(ii) remains the index of angle nearest to angle(ii) within 
    % those which is larger than the angle(ii) and in angleN
    if isempty(angleN(angleN - angle(ii) >= 0))
        % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
        tmp = angleN - angle(ii) + 2*pi >= 0;
    else
        tmp = angleN - angle(ii) >= 0;
    end
    ind2 = sparTheta*nLoop + find(angleN == min(angleN(tmp)));
    % ind3(ii) remains the index of angle nearest to angle(ii) within 
    % those which is smaller than the angle(ii) and in angleN
    if isempty(angleN(angleN - angle(ii) < 0))
        tmp = angleN - angle(ii) - 2*pi < 0;
    else
        tmp = angleN - angle(ii) < 0;
    end
    ind3 = sparTheta*nLoop + find(angleN == max(angleN(tmp)));
    [res(ii),peakPt(:,ii),uLim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp,toolRadius, ...
        uLim(:,ii),ii,ind2,ind3);
end
% inner side of each point on the tool path
parfor ii = (sparTheta + 1):ptNum
    % 如果是沿同一个极径的，就可以直接不用投影；否则还是需要这样子找
    nLoop = floor((ii - 1)/sparTheta) - 1;
    angleN = angle(sparTheta*nLoop + 1:sparTheta*(nLoop + 1));
    % ind2(ii) remains the index of angle nearest to angle(ii) within 
    % those which is larger than the angle(ii) and in angleN
    if isempty(angleN(angleN - angle(ii) >= 0))
        % to avoid that angle(ii) is cloesd to -pi, and smaller than each elements
        tmp = angleN - angle(ii) + 2*pi >= 0;
    else
        tmp = angleN - angle(ii) >= 0;
    end
    ind2 = sparTheta*nLoop + find(angleN == min(angleN(tmp)));
    % ind3(ii) remains the index of angle nearest to angle(ii) within 
    % those which is smaller than the angle(ii) and in angleN
    if isempty(angleN(angleN - angle(ii) < 0))
        tmp = angleN - angle(ii) - 2*pi < 0;
    else
        tmp = angleN - angle(ii) < 0;
    end
    ind3 = sparTheta*nLoop + find(angleN == max(angleN(tmp)));
    [res(ii + resNum),peakPt(:,ii + resNum),uLim(:,ii)] = residual3D( ...
        toolPathPt,toolNormDirect,toolCutDirect,toolContactU,toolSp,toolRadius, ...
        uLim(:,ii),ii,ind2,ind3);
end
tRes = toc;
fprintf('The time spent in the residual height calculation process is %fs.\n',tRes);

% post-processing of the residual height data
% ( this process is in the visualization part now)

% for ii = 1:ptNum
%     uLim(1,ii) = max([uLim(1,ii),uLimTmp(1,ind2 == ii)]);
%     uLim(2,ii) = min([uLim(2,ii),uLimTmp(2,ind2 == ii)]);
% end
clear angle ind2 ind3 uLimTmp;

%%
figure('Name','residual height calculation');
plotSpar = 1;
tic
plot3(toolPathPt(1,1:plotSpar:end), ...
    toolPathPt(2,1:plotSpar:end), ...
    toolPathPt(3,1:plotSpar:end), ...
    '.','MarkerSize',6,'Color',[0,0.4470,0.7410]);
hold on;
% quiver3(toolPathPt(1,1:plotSpar:end), ...
%     toolPathPt(2,1:plotSpar:end), ...
%     toolPathPt(3,1:plotSpar:end), ...
%     toolCutDirect(1,1:plotSpar:end), ...
%     toolCutDirect(2,1:plotSpar:end), ...
%     toolCutDirect(3,1:plotSpar:end), ...
%     'AutoScale','on','Color',[0.6350,0.0780,0.1840]);
% quiver3(toolPathPt(1,1:plotSpar:end), ...
%     toolPathPt(2,1:plotSpar:end), ...
%     toolPathPt(3,1:plotSpar:end), ...
%     toolNormDirect(1,1:plotSpar:end), ...
%     toolNormDirect(2,1:plotSpar:end), ...
%     toolNormDirect(3,1:plotSpar:end), ...
%     'AutoScale','on','Color',[0.85,0.33,0.10]);
surf( ...
    surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3), ...
    'FaceColor','flat','FaceAlpha',1,'LineStyle','none');
colormap('summer');
cb = colorbar;
cb.Label.String = ['Height (',unit,')'];
for ii = 1:ptNum
    toolSp = toolData.toolBform;
    toolCoefs = toolData.toolBform.coefs;
    toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
    Q = fnval(toolSp,uLim(1,ii):0.01:uLim(2,ii));
    plot3(Q(1,:),Q(2,:),Q(3,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5); hold on;
end
axis equal; grid on;
set(gca,'FontSize',textFontSize,'FontName',textFontType,'ZDir','reverse');
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);
% legend('tool center point','tool cutting direction', ...
%     'tool spindle direction','','tool edge','Location','northeast');
legend('tool center point','','tool edge','Location','northeast');
tPlot = toc;
fprintf('The time spent in the residual height plotting process is %fs.\n',tPlot);


%% Visualization & Simulation
while true
    msgfig = questdlg({'Residual Height was calculated successfully!', ...
        'Ready for residual height visualization or machining simulation?'}, ...
        'Tool Path Simulation','Residual height','Machining simulation', ...
        'Cancel and quit',msgOpts);
    switch msgfig
        case 'Cancel and quit'
            % delete(parObj);
            profile off
            % profile viewer
            % profsave(profile("info"),"profile_data");
            % rmpath(genpath('funcs'));
            return;
        case 'Residual height'
            tic
            plotNum = 1000;
            xPlot = linspace(min(peakPt(1,:)),max(peakPt(1,:)),plotNum);
            yPlot = linspace(min(peakPt(2,:)),max(peakPt(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [resUnique,peakPtUnique] = groupsummary(res',peakPt(1:2,:)',@max);
            resMesh = griddata(peakPtUnique{1},peakPtUnique{2},resUnique,xMesh,yMesh);
            figure('Name','Residual height');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            surf(xMesh,yMesh,resMesh,'EdgeColor','interp'); hold on;
            colormap("parula");
            grid on;
            plot3(peakPt(1,:),peakPt(2,:),res,'o', ...
                'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
            % cb1 = colorbar;
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['residual height (',unit,')']);
            legend('residual height in each peakPt','residual height map', ...
                'Location','best');
            nexttile;
            contourf(xMesh,yMesh,resMesh); hold on;
            colormap("parula");
            axis equal; grid on;
            cb2 = colorbar;
            cb2.Label.String = ['Residual Height (',unit,')'];
            cb2.Layout.Tile = 'east';
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            tRes = toc;
            fprintf('The time spent in the residual map process is %fs.\n',tRes);
        case 'Machining simulation'
            stepLength = 0.01;
            nLoop = ceil(ptNum/sparTheta);
            uLimRound = round(uLim,2);
            toolPathList = [];
            tic
            for ii = 1:nLoop % each loop
                Q = cell(sparTheta,1);
                for jj = 1:sparTheta
                    toolSp = toolData.toolBform;
                    toolSp.coefs = quat2rotm(toolQuat((ii-1)*nLoop+jj,:))*toolCoefs + toolVec(:,(ii-1)*nLoop+jj);
                    tmp = fnval(toolSp,uLimRound(1,(ii-1)*nLoop+jj):stepLength:uLimRound(2,(ii-1)*nLoop+jj));
                    % Q{jj} = tmp;
                    toolPathList = [toolPathList,tmp];
                end
%                 for u = 0:stepLength:1
%                     for jj = 1:sparTheta
%                         if u >= uLimRound(1,(ii-1)*nLoop+jj) && u <= uLimRound(2,(ii-1)*nLoop+jj)
%                             tmp = Q{jj}(:,round((u - uLimRound(1,(ii-1)*nLoop+jj))/stepLength + 1));
%                             toolPathList = [toolPathList,tmp];
%                         end
%                     end
%                 end
            end
            % plot3(toolPathList(1,:),toolPathList(2,:),toolPathList(3,:),'.','Color',[0,0.4450,0.7410]);
            plotNum = 1000;
            xPlot = linspace(min(toolPathList(1,:)),max(toolPathList(1,:)),plotNum);
            yPlot = linspace(min(toolPathList(2,:)),max(toolPathList(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [toolPathZUnique,toolPathXYUnique] = groupsummary(toolPathList(3,:)',toolPathList(1:2,:)',@max);
            zMesh = griddata(toolPathXYUnique{1},toolPathXYUnique{2},toolPathZUnique,xMesh,yMesh);
            % calculate the error based on the designed surface
            z0Mesh = griddata(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3),xMesh,yMesh);
            % 这里有更好的仿真方式：每个点都计算到设计曲面的距离，而不是沿z方向的距离！！！！
            figure('Name','Machining surface simulation');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            mesh(xMesh,yMesh,zMesh,'EdgeColor','interp'); hold on;
            colormap("jet");
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['z (',unit,')']);
            grid on;
            nexttile;
            mesh(xMesh,yMesh,zMesh - z0Mesh,'EdgeColor','interp'); hold on;
            colormap("jet");
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['z error (',unit,')']);
            grid on;
            tSimul = toc;
            fprintf('The time spent in the simulation calculation process is %fs.\n',tSimul);
    end
end
