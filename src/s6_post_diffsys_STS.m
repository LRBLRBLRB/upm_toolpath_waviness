% to load the diffsys-exported cnc file and modified the head and tail of
% it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Thanks to the Nanocamfile_read_211122.m, written by Doctor Junnan Chen, 
%  I wrote the post processing DIFFSYS Jobfile 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc
addpath(genpath('funcs'));

%% modification of process data
unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

msgMode.WindowStyle = 'non-modal';
msgMode.Interpreter = 'tex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tool
toolRadius = 0.18746;

% surface
c = 0.69;
surfRange = 0.6;
z0 = 0;
syms x;
surfSym = c.*x.^2./(1 + sqrt(1 - c.^2.*x.^2)) + z0;
surfFunc = matlabFunction(surfSym);

% cnc
loop.num = 1;
loop.offset = 0;
loop.step = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the CNC file
workspaceDir = uigetdir(fullfile(['D:\OneDrive - sjtu.edu.cn\Research\Projects' ...
    '\202111-考虑刀具几何的路径规划\experiment\非球面加工\20230510\cnc']), ...
    'Select the workspace directory:');
[jobFileName,jobDirName] = uigetfile({ ...
    '*.pgm','DIFFSYS Jobfile(*.pgm)'; ...
    '*.nc','CNC-files(*.nc)'; ...
    '*.*','All Files(*.*)'}, ...
    'Select one CNC data file', ...
    fullfile(workspaceDir,'spiralpath.nc'), ...
    'MultiSelect','off');
if ~jobFileName
    msgbox(sprintf('\\fontname{%s}\\fontsize{%d} No CNC file saved.', ...
        textFontType,textFontSize),'Message','warn',msgMode);
end
jobPath = fullfile(jobDirName,jobFileName);

fprintf('Importing data...\r\n');
jobID = fopen(jobPath,'r');

tic
maxCounter = 0;
while ~feof(jobID)
    % read 10000 str, calculate the number of \n (char(10)=\n)
    % '*char' indicates 1 str per read, *indicates output str
    maxCounter = maxCounter + sum(fread(jobID,10000,'*char') == newline);
end
fclose(jobID);
fprintf('File row number is %d.\r\n',maxCounter);
toc

tic
jobID = fopen(jobPath,'r');
waitBar1 = waitbar(0,'Data Importing ...','Name','CNC Data Load', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(waitBar1,'canceling',0);
cncLine = zeros(maxCounter,1);
cncLine = string(cncLine);
for ii = 1 : maxCounter
    % Check for clicked Cancel button
    if getappdata(waitBar1,'canceling')
        break;
    end
    if mod(ii,fix(maxCounter/100)) == 0
        displayData = num2str(roundn(ii/maxCounter*100,-2));
        % Calculate percentage
        displayStr = ['Import Progress: ',displayData,'%'];
        % Show Calculate State
        waitbar(ii/maxCounter,waitBar1,displayStr);
        % Progress bar dynamic display
    end
    str = fgetl(jobID);
    cncLine(ii) = string(str);
end

waitbar(1,waitBar1,'Data Imported ! ');
pause(1);
delete(waitBar1);   % Close Progress bar window
fclose(jobID);
fprintf('Data imported.\r\n')
toc

%% delete the original header & read the main part
semicolonInd = startsWith(cncLine,';');
cncLine(semicolonInd) = [];
cncNum = length(cncLine);

axisC = zeros(cncNum,1);
axisX = zeros(cncNum,1);
axisZ = zeros(cncNum,1);
for ii = 1:cncNum
    data = sscanf(cncLine(ii),'C%fX%fZ%f');
    axisC(ii) = data(1);
    axisX(ii) = data(2);
    axisZ(ii) = data(3);
end

% diffsys - nanocam convertion
axisZ = axisZ - axisZ(end);
axisX = -1.*axisX;
axisC = wrapTo360(-1.*axisC);
axisC(find(abs(axisC - 360) < 1e-3)) = 0;

%% simulate the actual surface
% waitBar2 = waitbar(0,'Figure Plotting ...','Name','CNC Data Plot', ...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(waitBar2,'canceling',0);
% 
% figure;
% rSpar = linspace(0,axisX(1),cncNum/10);
% plot(rSpar,surfFunc(rSpar),'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on;
% plotNum = 200;
% plotRange = 1:cncNum;
% plotRange = plotRange(mod(1:cncNum,fix(cncNum/plotNum)) == 1);
% for ii = plotRange
%     % Check for clicked Cancel button
%     if getappdata(waitBar2,'canceling')
%         break;
%     end
%     displayData = num2str(roundn(ii/(2*cncNum)*100,-2)); % Calculate percentage
%     displayStr = ['Figure Plotting ... ',displayData,'%']; % Show Calculate State
%     waitbar(ii/cncNum,waitBar2,displayStr); % Progress bar dynamic display
% 
%     rectangle('Position',[axisX(ii)-toolRadius,axisZ(ii)-toolRadius, ...
%         2*toolRadius,2*toolRadius],'Curvature',[1,1],'EdgeColor',[0 0.4470 0.7410],'LineWidth',0.2);
% end
% 
% for ii = plotRange
%     % Check for clicked Cancel button
%     if getappdata(waitBar2,'canceling')
%         break;
%     end
%     displayData = num2str(roundn((ii + cncNum)/(2*cncNum)*100,-2)); % Calculate percentage
%     displayStr = ['Figure Plotting ...',displayData,'%']; % Show Calculate State
%     waitbar(ii/cncNum,waitBar2,displayStr); % Progress bar dynamic display
% 
%     scatter(axisX(ii),axisZ(ii),'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerEdgeAlpha',1);
%     % plot(axisX(ii),axisZ(ii),'Color',[0 0.4470 0.7410],'LineWidth',2);
% end
% set(gca,'FontSize',textFontSize,'FontName',textFontType);
% xlabel(['r (',unit,')']);
% ylabel(['z (',unit,')']);
% 
% delete(waitBar2);   

%% export the right file
% [configFileName,configDirName] = uigetfile({ ...
%     '*.pgm','Text file(*.txt)';
%     '*,*','All Files(*.*)'}, ...
%     'Select one CNC configuration file', ...
%     fullfile(src,'DIFFSYS_2D.nc'), ...
%     'MultiSelect','off');
% configPath = fullfile(configDirName,configFileName);
% configFid = fopen(configPath,'r');
% 
% fclose(configFid);

% cncFileName = getlastfoldername(jobPath);
[cncFile,cncDir] = uiputfile({'*.nc','Numerical control files(*.nc)'; ...
    '*.*','All files'},'Enter the file to save the CNC code',fullfile( ...
    jobDirName,[jobFileName,'(',datestr(now,'yyyymmddTHHMMSS'),').nc']));
cncPath = fullfile(cncDir,cncFile);
if ~cncFile
    msgbox(sprintf('\\fontname{%s}\\fontsize{%d} No CNC file saved.', ...
        textFontType,textFontSize),'Message','warn',msgMode);
end

writecnc_STS(cncPath,'G55','T0303',axisC,axisX,axisZ,loop);
