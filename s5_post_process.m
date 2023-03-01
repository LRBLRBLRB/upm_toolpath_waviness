% This programme aims to export the toolpath file, for the continuing post
% processing procedure in the IMSpost software to get the nc file.
isAPP = false;
if isAPP
    workspaceDir = app.workspace;
    spiralPath = app.spiralPath;
    spiralNorm = app.spiralNorm;
else
    close all;
    addpath(genpath('funcs'));
    workspaceDir = 'workspace\20220925-contrast\nagayama_concentric';
    if ~(exist('spiralPath','var') && exist('spiralNorm','var'))
        % the spiral path does not exist in the workspace
        [fileName,dirName] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','All Files(*.*)'}, ...
            'Select one tool path data file', ...
            fullfile(workspaceDir,'tooltheo.mat'), ...
            'MultiSelect','off');
        toolPathPath = fullfile(dirName,fileName);
        % toolName = 'output_data\tool\toolTheo_3D.mat';
        spiralPath = load(toolPathPath,'spiralPath');
        spiralNorm = load(toolPathPath,'spiralNorm');
    end
end

%% case 1: export apt file

[aptFname,aptPath,aptInd] = uiputfile( ...
    {'*.apt','Automatically Programmed Tool(*.apt)';'*.*','All Files'}, ...
    'Enter the file to save the tool path','toolPath.apt');
aptFile = fullfile(aptPath,aptFname);
aptFid = fopen(aptFile,'w');
for ii = 1:size(spiralPath,2)
    fprintf(aptFid,'GOTO/%f,%f,%f,%f,%f,%f\n', ...
        spiralPath(1,ii),spiralPath(2,ii),spiralPath(3,ii), ...
        spiralNorm(1,ii),spiralNorm(2,ii),spiralNorm(3,ii));
end

%% case 2: post-processing to export the axial position

%%
% rmpath(genpath('funcs'));