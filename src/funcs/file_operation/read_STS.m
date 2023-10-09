function [axisC,axisX,axisZ] = read_STS(workspaceDir,cncFormat,unit)
%IMPORTNC read the 3-axes CL points from cnc file, and only support the
%nanotech 650FG V2

%% load the cnc file
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

fprintf('Importing data...\n');
jobID = fopen(jobPath,'r');

tic
maxCounter = 0;
while ~feof(jobID)
    % read 10000 str, calculate the number of \n (char(10)=\n)
    % '*char' indicates 1 str per read, *indicates output str
    maxCounter = maxCounter + sum(fread(jobID,10000,'*char') == newline);
end
fclose(jobID);
fprintf('File row number is %d, with time consumption of %fs.\n',maxCounter,toc);

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
fprintf('Data imported, with time consumption of %fs\n',toc);

%% delete the original header & read the main part
semicolonInd = startsWith(cncLine,';');
cncLine(semicolonInd) = [];
cncNum = length(cncLine);

data = zeros(cncNum,3);
for ii = 1:cncNum
    data(ii,:) = sscanf(cncLine(ii),cncFormat);
end

%% diffsys - nanocam convertion

unitList = {'m','mm','\mum','nm'};
presUnit = find(strcmp(unitList,'mm'),1);
aimUnit = find(strcmp(unitList,unit),1);
1000^(aimUnit - presUnit)*data;


axisC = data(:,1);
axisX = data(:,2);
axisZ = data(:,3);

axisZ = axisZ - axisZ(end);
axisX = -1.*axisX;
axisC = wrapTo360(-1.*axisC);
axisC(find(abs(axisC - 360) < 1e-3)) = 0;

end