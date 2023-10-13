function [cncData,jobFileName,jobDirName] = read_my_STS(workspaceDir,cncFormat)
%IMPORTNC read the 3-axes CL points from cnc file, which only supports the
%nanotech 650FG V2 and the cnc files that have been rewritten.

%% load the cnc file
[jobFileName,jobDirName] = uigetfile({ ...
    '*.pgm','DIFFSYS Jobfile(*.pgm)'; ...
    '*.nc','CNC-files(*.nc)'; ...
    '*.*','All Files(*.*)'}, ...
    'Select one CNC data file', ...
    fullfile(workspaceDir,'spiralpath.nc'), ...
    'MultiSelect','off');
if ~jobFileName
    msgbox('No CNC file saved.','Message','warn','non-modal');
end
jobPath = fullfile(jobDirName,jobFileName);

fprintf('Importing data...\n');
jobFid = fopen(jobPath,'r');

% get rid of the header of the nc file
numHeader = 1;
while ~feof(jobFid)
    tmpLine = fgetl(jobFid);
    if strncmp(tmpLine,'#105',4)
        feedVel = sscanf(tmpLine,'#105=%d%.s');
        continue;
    end
    if strncmp(tmpLine,'#201',4)
        spindleVel = sscanf(tmpLine,'#201=%d%.s');
        continue;
    end
    % if the line begins with %d%d or -%d, then break
    if strcmp(tmpLine,'( CUTTING BLOCK )')
        break;
    end
    numHeader = numHeader + 1;
end

% load the X Z data
cncData = zeros(3,0);
tmpLine = fgetl(jobFid);
while ~feof(jobFid)
    tmpLine = fgetl(jobFid);
    % if the line begins with %d%d or -%d, then break
    if strcmp(tmpLine,'(linking block)') || strcmp(tmpLine,'(LINKING BLOCK)')
        break;
    end
    cncData = [cncData,sscanf(tmpLine,cncFormat)];
end
fclose(jobFid);

end