% This programme aims to export the toolpath file, for the continuing post
% processing procedure in the IMSpost software to get the nc file.
close all;
%     clear;

unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

addpath(genpath('funcs'));
if ~(exist('spiralPath','var') && exist('spiralNorm','var'))
    workspaceDir = uigetdir(fullfile('..','workspace','20230504 D906'),'select the workspace directory:');
    % the spiral path does not exist in the workspace
    [fileName,dirName] = uigetfile({ ...
        '*.mat','MAT-files(*.mat)'; ...
        '*,*','All Files(*.*)'}, ...
        'Select one tool path data file', ...
        fullfile(workspaceDir,'spiralpath.mat'), ...
        'MultiSelect','off');
    toolPathPath = fullfile(dirName,fileName);
    processData = load(toolPathPath);
    % toolName = 'output_data\tool\toolTheo_3D.mat';
    spiralAngle = processData.spiralAngle;
    spiralPath = processData.spiralPath;
    spiralNorm = processData.spiralNorm;
    spiralCut = processData.spiralCut;
    spiralQuat = processData.spiralQuat;
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
axisX = -1.*axisX;
axisC = wrapTo360(-1.*axisC);
axisC(find(abs(axisC - 360) < 1e-3)) = 0;
%         end

% cnc header parameter
% app = post_process;

%% simulation of the tool path 
waitBar2 = waitbar(0,'Figure Plotting ...','Name','CNC Data Plot', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(waitBar2,'canceling',0);

figure;
% rSpar = linspace(0,surfRange,cncNum/10);
% plot(rSpar,surfFunc(rSpar),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on;
cncNum = length(axisC);
x = axisX.*cosd(axisC);
y = axisX.*sind(axisC);
Pt = scatter3(x(1),y(1),axisZ(1),6, ...
    'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]); hold on;
for ii = 2:cncNum
    % Check for clicked Cancel button
    if getappdata(waitBar2,'canceling')
        break;
    end

    delete(Pt);
    Pt = scatter3(x(ii),y(ii),axisZ(ii),6, ...
        'MarkerEdgeColor','flat','MarkerFaceColor',[0.8500 0.3250 0.0980]);
    line([x(ii - 1),x(ii)],[y(ii - 1),y(ii)],[axisZ(ii - 1),axisZ(ii)], ...
        'Color',[0,0.4470,0.7410],'LineWidth',0.1); hold on;

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
    msgbox(sprintf('\\fontname{%s}\\fontsize{%d} No CNC file saved.', ...
        textFontType,textFontSize),'Message','warn',msgMode);
end
ncFile = fullfile(ncPath,ncFname);
loop.num = 1;
loop.offset = 0;
loop.step = 0;
writecnc_STS(ncFile,'G55','T0303',axisC,axisX,axisZ,loop);

%%
% rmpath(genpath('funcs'));