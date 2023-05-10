% 切削仿真：考虑刀具几何误差，仿真计算加工自由曲面的结果

close all;
clear; clc;

%% get the code from the nc file
% dirName = 'tool_path\';
% fileName = 'Lens array_210721_RT0.2_tilt0_rot0_SF_s6a12d0.66.nc';
if ~exist('workspaceDir','var')
    workspaceDir = fullfile('..','workspace','';
end
[fileName,dirName] = uigetfile({ ...
    '*.nc','nc-file(*.nc)'; ...
    '*.txt','txt-file(*.txt)'; ...
    '*.*','All the files(*.*)' ...
    },'Choose a file',workspaceDir);
pathName = fullfile(dirName,fileName);
ncFile = fopen(pathName);

% get each line of the nc code (using file I/O)
% tic
% lineNum = 1;
% while ~feof(ncFile)
%     tmpLine = fgets(ncFile); % extract the line (fgets faster than fgetl)
%     if contains(tmpLine,'(') % delete useless lines
%         continue;
%     end
%     ncLine{lineNum} = tmpLine;
%     lineNum = lineNum + 1;
% end
% fclose(ncFile);
% tFile = toc
% 
% axesPos = zeros(1,5); % columns mean X, Y, Z, B, C coor'
% tic
% for ii = 2:lineNum
%     lineLength = length(ncLine{ii});
%     % template of the position of motion axes
%     pat = digitsPattern + "." + digitsPattern | digitsPattern;
%     axesPos(ii,1) = getPos(axesPos(ii-1,1), ncLine{ii}, "X", pat);
%     axesPos(ii,2) = getPos(axesPos(ii-1,2), ncLine{ii}, "Y", pat);
%     axesPos(ii,3) = getPos(axesPos(ii-1,3), ncLine{ii}, "Z", pat);
%     axesPos(ii,4) = getPos(axesPos(ii-1,4), ncLine{ii}, "B", pat);
%     axesPos(ii,5) = getPos(axesPos(ii-1,5), ncLine{ii}, "C", pat);
% end
% tExt = toc

% get the matrix of each point (using readtable)
axesPos = importNCFile(pathName,[53 2214856]);
len = size(axesPos,1);
axesPos(1,5) = 0.217; % Z position
for ii = 2:len
    if isnan(axesPos(ii,3))
        axesPos(ii,5) = axesPos(ii-1,5);
    else
        axesPos(ii,5) = axesPos(ii,3);
    end
end
axesPos(:,4) = zeros(len,1); % Y position
axesPos(:,3) = axesPos(:,2);
axesPos(:,2) = axesPos(:,1);
axesPos(:,1) = zeros(len,1);

%% Nanotech 650FG kinematics
syms B C X Y Z;
TB = [cos(B),0,-sin(B),0; 0,1,0,0; sin(B),0,cos(B),0; 0,0,0,1];
TC = [cos(C),-sin(C),0,0; sin(C),cos(C),0,0; 0,0,1,0; 0,0,0,1];
TX = [1,0,0,X; 0,1,0,0; 0,0,1,0; 0,0,0,1];
TY = [1,0,0,0; 0,1,0,Y; 0,0,1,0; 0,0,0,1];
TZ = [1,0,0,0; 0,1,0,0; 0,0,1,Z; 0,0,0,1];
T = (TB*TZ)\(TC*TX*TZ);
endPosion = zeros(len,3);
endPose = zeros(len,3);
parfor ii = 1:len
    Ti = subs(T,[B,C,X,Y,Z],axesPos(ii,:));
    endPosion(ii,:) = Ti(1:3,4);
    tmp = logm(Ti(1:3,1:3));
    endPose(ii,:) = [tmp(3,2),tmp(1,3),tmp(2,1)];
end

% surf(endPos(:,1),endPos(:,2),endPos(:,3));

%% simulation situation one: 




%% simulation situation two: planned trajectory
% equidistant curves

%% functions
function axesPos = getPos(axesForm,ncLine,Cha,pat)
tmp = str2double(extract(extract(ncLine, Cha + pat),pat));
if isempty(tmp)
    axesPos = axesForm;
else
    axesPos = tmp;
end
end

