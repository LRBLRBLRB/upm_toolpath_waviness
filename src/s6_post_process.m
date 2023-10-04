% This programme aims to export the toolpath file, for the continuing post
% processing procedure in the IMSpost software to get the nc file.
close all;
%     clear;

unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

addpath(genpath('funcs'));
if ~(exist('spiralPath','var') && exist('spiralNorm','var'))
    workspaceDir = uigetdir(fullfile('..','workspace'),'select the workspace directory:');
    % the spiral path does not exist in the workspace
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','All Files(*.*)'}, ...
        'Select one tool path data file', ...
        fullfile(workspaceDir,'spiralpath.mat'), ...
        'MultiSelect','off');
    toolPathPath = fullfile(dirName,fileName);
    load(toolPathPath);
%     processData = load(toolPathPath);
%     spiralAngle = processData.spiralAngle;
%     spiralPath = processData.spiralPath;
%     spiralNorm = processData.spiralNorm;
%     spiralCut = processData.spiralCut;
%     spiralQuat = processData.spiralQuat;
%     surfMesh = processData.surfMesh;
%     approxMethod = processData.approxMethod;
end

%% generate the 5-axis tool path from original data
% toolpath should be saved in "x y z i j k" format

spiralAngle1 = 180/pi*spiralAngle; % rad -> deg
spiralPath1 = 0.001*spiralPath; % um -> mm
postType = 2;

%% post-processing to export the axial position
%         axisC = atan2(2*(spiralQuat(:,2).*spiralQuat(:,3) - spiralQuat(:,1).*spiralQuat(:,4)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,3).^2) - 1);
%         axisB = atan2(-2*(spiralQuat(:,2).*spiralQuat(:,4) + spiralQuat(:,1).*spiralQuat(:,3)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,4).^2) - 1);
%         axisZ = (spiralPath1(3,:).')./cos(axisB);
%         axisX = axisZ.*sin(axisB).*cos(axisC) - spiralPath1(1,:).';
%         axisY = axisZ.*sin(axisB).*sin(axisC) - spiralPath1(2,:).';

axisC = (wrapTo360(spiralAngle1));
axisZ = spiralPath1(3,:);
axisX = sign(spiralPath1(1,1))*vecnorm(spiralPath1(1:2,:),2,1);

% zero point of the C axis
if axisC(1) ~= 0
    axisC = axisC - axisC(1);
    axisC = wrapTo360(axisC);
end

% zero point of z axis
if axisZ(end) ~= 0
    axisZ = axisZ - axisZ(end);
end

% direction correction
%         if strcmp(startDirection,'X Plus')
% axisX = -1.*axisX;
% axisC = wrapTo360(-1.*axisC);
axisC(find(abs(axisC - 360) < 1e-3)) = 0;
%         end

% cnc header parameter
% app = post_process;

%% simulation of the tool path 
waitBar2 = waitbar(0,'Figure Plotting ...','Name','CNC Data Plot', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(waitBar2,'canceling',0);

figure; 
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);

ax1 = nexttile;
surf(ax1, ...
    surfMesh(:,:,1)*0.001,surfMesh(:,:,3)*0.001,-1*surfMesh(:,:,2)*0.001, ...
    'FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
hold(ax1,'on');
colormap(ax1,'summer');
xlabel(ax1,'x'); ylabel(ax1,'z'); zlabel(ax1,'y');
view(ax1,-255,15);

ax2 = nexttile;
surf(ax2, ...
    surfMesh(:,:,1)*0.001,surfMesh(:,:,3)*0.001,surfMesh(:,:,2)*0.001, ...
    'FaceColor','flat','FaceAlpha',0.3,'LineStyle','none');
hold(ax2,'on');
colormap(ax2,'summer');
xlabel(ax2,'x'); ylabel(ax2,'z'); zlabel(ax2,'y');
view(ax2,-255,15);

% tool
toolSp = toolData.toolBform;
if strcmp(startDirection,'X Plus')
    toolSp.coefs = toolSp.coefs*0.001 + [0.3;0;0.5];
else
    toolSp.coefs = toolSp.coefs*0.001 + [-0.3;0;0.5];
end
toolSpPt = fnval(toolSp,0:0.01:1);
TEdge = plot3(ax2,toolSpPt(2,:),toolSpPt(3,:),toolSpPt(1,:),'Color',[0,0.4470,0.7410]);
TFill = patch(ax2,'XData',toolSpPt(2,:),'YData',toolSpPt(3,:),'ZData',toolSpPt(1,:),...
    'EdgeColor','none','FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.3);

cncNum = length(axisC);
X = axisX.*cosd(axisC);
Y = axisX.*sind(axisC);
Z = axisZ;

Pt1 = scatter3(ax1,spiralPath1(1,1),spiralPath1(3,1),-1*spiralPath1(2,1),12, ...
    'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]);
Pt2 = scatter3(ax2,X(1),Z(1),Y(1),12, ...
    'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]);
for ii = 2:cncNum
    % Check for clicked Cancel button
    if getappdata(waitBar2,'canceling')
        break;
    end

    delete(Pt1);
    delete(Pt2);
    Pt1 = scatter3(ax1,spiralPath1(1,ii),spiralPath1(3,ii),-1*spiralPath1(2,ii),12, ...
        'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]);
    line(ax1,[spiralPath1(1,ii - 1),spiralPath1(1,ii)], ...
        [spiralPath1(3,ii - 1),spiralPath1(3,ii)], ...
        [-1*spiralPath1(2,ii - 1),-1*spiralPath1(2,ii)], ...
        'Color',[0,0.4470,0.7410],'LineWidth',0.1);
    Pt2 = scatter3(ax2,X(ii),Z(ii),Y(ii),6, ...
        'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]);
    line(ax2,[X(ii - 1),X(ii)],[Z(ii - 1),Z(ii)],[Y(ii - 1),Y(ii)], ...
        'Color',[0,0.4470,0.7410],'LineWidth',0.1);

    displayData = num2str(roundn(ii/cncNum*100,-2)); % Calculate percentage
    displayStr = ['Figure Plotting ... (',num2str(ii),'/', ...
        num2str(cncNum),') ',displayData,'%']; % Show Calculate State
    waitbar(ii/cncNum,waitBar2,displayStr); % Progress bar dynamic display
end

delete(waitBar2); 


%% put in the .nc file
spiralncFolderName = getlastfoldername(workspaceDir);
[ncFname,ncPath,ncInd] = uiputfile( ...
    {'*.nc','Numerical control files(*.nc)';'*.*','All files'}, ...
    'Enter the file to save the CNC code',fullfile( ...
    workspaceDir,[spiralncFolderName,'-spiralPath-',approxMethod,datestr(now,'yyyymmddTHHMMSS'),'.nc']));
if ~ncFname
    msgbox('No CNC file saved.','Message','warn','non-modal');
    return;
end
ncFile = fullfile(ncPath,ncFname);
loop.num = 1;
loop.offset = 0;
loop.step = 0;
writecnc_STS(ncFile,'G55','T0303',axisC,axisX,axisZ,loop);

%%
% rmpath(genpath('funcs'));