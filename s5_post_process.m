% This programme aims to export the toolpath file, for the continuing post
% processing procedure in the IMSpost software to get the nc file.
isAPP = false;
if isAPP
    workspaceDir = app.workspace;
    if app.isSpiralFile
        toolPathPath = fullfile(app.dirName,app.fileName);
        % toolName = 'output_data\tool\toolTheo_3D.mat';
        spiralPath = load(toolPathPath,'spiralPath');
        spiralNorm = load(toolPathPath,'spiralNorm');
        spiralCut = load(toolPathPath,'spiralCut');
        spiralQuat = load(toolPathPath,'spiralQuat');
    else
        spiralPath = app.spiralPath;
        spiralNorm = app.spiralNorm;
        spiralCut = app.spiralCut;
        spiralQuat = app.spiralQuat;
    end
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
        spiralCut = load(toolPathPath,'spiralCut');
        spiralQuat = load(toolPathPath,'spiralQuat');
    end
end

%% generate the 5-axis tool path from original data
% toolpath should be saved in "x y z i j k" format



switch postType
    case 1
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
    case 2
        %% case 2: post-processing to export the axial position
        p
end

%%
% rmpath(genpath('funcs'));