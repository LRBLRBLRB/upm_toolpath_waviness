% to show the tool measuring data
%% 2D curve results
% close all;
clear; clc;
isAPP = false;
addpath(genpath('funcs'));

% global variables
% global textFontSize textFontType unit fitMethod paramMethod;
% workspaceDir = '../workspace/20221020-tooltip';
workspaceDir = '../workspace/20230424';
unit = '\mum';
textFontSize = 12;
textFontType = 'Times New Roman';

toolUnit = 'mm';

%% tool measurement file loading
while true
    [fileName,dirName] = uigetfile({ ...
        '*.csv','Comma-Separated Values-files(*.csv)'; ...
        '*.txt','text-files(*.txt)'; ...
        '*.*','all files(*.*)'...
        }, ...
        'Select one tool tip measurement data', ...
        workspaceDir, ...
        'MultiSelect','off');
    pathName = fullfile(dirName,fileName);
    if ~fileName
        fprintf('No tool data file has been loaded.\n');
        return;
    end
    workspaceDir = dirName;
    
    [~,~,fileExt] = fileparts(pathName);
    switch fileExt
        case {'.csv'}
            % get rid of the header of the csv file
            numHeader = 0;
            tooltipFile = fopen(pathName);
            while ~feof(tooltipFile)
                tmpLine = fgets(tooltipFile);
                % if the line begins with %d%d or -%d, then break
                if ~isnan(str2double(tmpLine(1:2)))
                    break;
                end
                numHeader = numHeader + 1;
            end
            fclose(tooltipFile);
            toolOri = importdata(pathName,',',numHeader);
            if size(toolOri,2) ~= 3 && size(toolOri,2) ~= 2
                toolOri = toolOri.data;
            end
            toolOri(:,3) = [];
%             toolOri = sortrows(toolOri,1,'ascend');
            toolOri = toolOri';
        case {'.txt'}
            % tool data file that has been processed in mmt software
            toolOri = importdata(pathName,' ',0);
%             toolOri = sortrows(toolOri,1,'ascend');
            toolOri = toolOri';
    end

    % unit convertion
    unitList = {'m','mm','\mum','nm'};
    presUnit = find(strcmp(unitList,toolUnit),1);
    aimUnit = find(strcmp(unitList,unit),1);
    toolOri = 1000^(aimUnit - presUnit)*toolOri;
    
    %% plot the importing result.
    dispInd = strfind(dirName,'workspace');
    if isempty(dispInd)
        dispFileName = fileName;
    else
        dispFileName = fullfile(dirName(dispInd + 10:end),fileName);
    end
    fig1 = figure('Name',dispFileName);
    ax1 = plot(toolOri(1,:),toolOri(2,:),'.','MarkerSize',2);
    hold on;
    grid on;
    xlabel(['x (',unit,')']);
    ylabel(['y (',unit,')']);

    msgfig = questdlg({'Tool data loaded successfully!', ...
        'Another tool data to load?'}, ...
        'Exit','Yes','Quit','Yes');
    if strcmp(msgfig,'Quit') || isempty(msgfig)
        return;
    end
end

