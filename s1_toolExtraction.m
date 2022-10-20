% to deal with the tool tip measurement data
% and to get the 3D point cloud of the tool tip

if true
    close all;
    clear; clc;
    addpath(genpath('funcs'));
    
    % global variables
    % global textFontSize textFontType unit fitMethod paramMethod;
    workspaceDir = 'workspace/20221020-tooltip';
    fitMethod = 'Levenberg-Marquardt';
    paramMethod = 'centripetal';
    unit = '\mum';
    textFontSize = 12;
    textFontType = 'Times New Roman';
end

%% 2D results
fileName = fullfile(workspaceDir,"tooltip result/20221019-strategy-2+40-5.csv");
data = load(fileName);
figure('Name','Original tool data');
plot(data(:,1),data(:,2),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);

% outliners removed

rmoutliers(data,2,"median");
% 如果多边形的直线段和圆弧段重合部分的数据差别比较大，说明目前的测量模式有问题。



% line approx and angle calculation





% tool tip arc extraction



%% 3D results
fileName = fullfile(workspaceDir,"tooltip result/20221019-strategy-3+40-3.csv");
% get rid of the header of the csv file
numHeader = 0;
tooltipFile = fopen(fileName);
while ~feof(tooltipFile)
    tmpLine = fgets(tooltipFile);
    if ~isnan(str2double(tmpLine(1:2)))
        break;
    end
    numHeader = numHeader + 1;
end
fclose(tooltipFile);
data = importdata(fileName,',',numHeader);
data = data.data;
figure('Name','Original tool data');
plot3(data(:,1),data(:,2),data(:,3),'.','MarkerSize',2);
hold on;
grid on;
xlabel(['x (',unit,')']);
ylabel(['y (',unit,')']);
zlabel(['z (',unit,')']);

figure('Name','Original 2D data');
















%%

% rmpath(genpath('.')