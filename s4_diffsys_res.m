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

% workspaceDir = 'workspace\20220925-contrast\nagayama_concentric';
% workspaceDir = 'workspace\20221020-tooltip\tooltip fitting result';
workspaceDir = 'workspace\20230417';
diaryFile = fullfile('workspace\diary',['diary',datestr(now,'yyyymmddTHHMMSS')]);
% diary diaryFile;
% diary on;

%% load nc file
[fileName,dirName] = uigetfile({ ...
    '*.nc;.pgm','CNC-files(*.nc,*pgm)'; ...
    '*,*','all files(*.*)'}, ...
    'Select one cnc file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');

ncName = fullfile(dirName,fileName);
fid = fopen(ncName,'r');
% get rid of the header of the csv file
numHeader = 1;
while ~feof(fid)
    tmpLine = fgetl(fid);
    if strncmp(tmpLine,'#105',4)
        cutVel = sscanf(tmpLine,'#105=%d%.s');
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
ncData = zeros(2,0);
while ~feof(fid)
    tmpLine = fgetl(fid);
    % if the line begins with %d%d or -%d, then break
    if strcmp(tmpLine,'(linking block)')
        break;
    end
    ncData = [ncData,sscanf(tmpLine,'X%fZ%f')];
end
% ncData = fscanf(fid,'X%fZ%f');
fclose(fid);

%% 2-axes convertion to cartesian
nData = size(ncData,2);
spiralPath = [nData(1);0;nData(3)];
for ii = 2:nData
    nData(1,ii)
end





