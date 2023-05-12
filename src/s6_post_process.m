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
%     clear;
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
end

feedRate = 0.001;


%% generate the 5-axis tool path from original data
% toolpath should be saved in "x y z i j k" format

spiralAngle1 = 180/pi*spiralAngle; % rad -> deg
spiralPath1 = 0.001*spiralPath; % um -> mm
postType = 2;
switch postType
    case 1
        %% case 1: export apt file
        [aptFname,aptPath,aptInd] = uiputfile( ...
            {'*.apt','Automatically Programmed Tool(*.apt)';'*.*','All Files'}, ...
            'Enter the file to save the tool path','toolPath.apt');
        aptFile = fullfile(aptPath,aptFname);
        aptFid = fopen(aptFile,'w');
        for ii = 1:size(spiralPath1,2)
            fprintf(aptFid,'GOTO/%f,%f,%f,%f,%f,%f\n', ...
                spiralPath1(1,ii),spiralPath1(2,ii),spiralPath1(3,ii), ...
                spiralNorm(1,ii),spiralNorm(2,ii),spiralNorm(3,ii));
        end
    case 2
        %% case 2: post-processing to export the axial position
%         axisC = atan2(2*(spiralQuat(:,2).*spiralQuat(:,3) - spiralQuat(:,1).*spiralQuat(:,4)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,3).^2) - 1);
%         axisB = atan2(-2*(spiralQuat(:,2).*spiralQuat(:,4) + spiralQuat(:,1).*spiralQuat(:,3)), ...
%             2*(spiralQuat(:,1).^2 + spiralQuat(:,4).^2) - 1);
%         axisZ = (spiralPath1(3,:).')./cos(axisB);
%         axisX = axisZ.*sin(axisB).*cos(axisC) - spiralPath1(1,:).';
%         axisY = axisZ.*sin(axisB).*sin(axisC) - spiralPath1(2,:).';

        axisC = (wrapTo360(spiralAngle1));
        axisZ = spiralPath1(3,:);
        axisX = vecnorm(spiralPath1(1:2,:),2,1);

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
        if ~isrev
            axisX = -1.*axisX;
            axisC = wrapTo360(-1.*axisC);
            axisC(find(abs(axisC - 360) < 1e-3)) = 0;
        end

        % cnc header parameter
        % app = post_process;

        % put in the .nc file
        spiralncFolderName = getlastfoldername(workspaceDir);
        [ncFname,ncPath,ncInd] = uiputfile( ...
            {'*.nc','Numerical control files(*.nc)';'*.*','All files'}, ...
            'Enter the file to save the CNC code',fullfile( ...
            workspaceDir,[spiralncFolderName,'spiralPath',datestr(now,'yyyymmddTHHMMSS'),'.nc']));
        ncFile = fullfile(ncPath,ncFname);
        ncFid = fopen(ncFile,'w');
        fprintf(ncFid,'#100 = 1\t\t( TOTAL PASSES OF Z)\n');
        fprintf(ncFid,'#103 = 0.0\t\t( DEPTH OF CUT Z )\n');
        fprintf(ncFid,'#555 = 0.001\t\t( FEED RATE )\n');

        fprintf(ncFid,'\nG52 G63 G71 G103 G40 G18 G90\n');
        fprintf(ncFid,'\n( Working Coordinates )\n%s\n%s\n\n','G55','T0303');

        fprintf(ncFid,'G18 G40 G94\n');
        fprintf(ncFid,'M78\n');
        fprintf(ncFid,'G94 Z20 F500\n\n');

        fprintf(ncFid,'M26.1\n\n');

        fprintf(ncFid,'#5 = 0\t\t( COUNT VARIABLE )\n');
        fprintf(ncFid,'#8 = 0\t\t( CUT OFFSET Z )\n');
        fprintf(ncFid,'WHILE[#5 LT #100] DO 2\n\n');

        fprintf(ncFid,'#8 = #8 - #103\t\t( Z AXIS CUTTING )\n');
        fprintf(ncFid,'G52 X0 Y0 Z[#8]\n');
        
        fprintf(ncFid,'( LINK BLOCK )\n');
        fprintf(ncFid,'G94 X%f\n',axisX(1));
        fprintf(ncFid,'G93 F0.001 C%f\n\n',axisC(1));
        
        fprintf(ncFid,'G94 Z%d F500\n',5 + floor(axisZ(1)));
        fprintf(ncFid,'Z%f F200\n',ceil(axisZ(1)));
        fprintf(ncFid,'Z%f F20\n',axisZ(1) + 0.1);
        fprintf(ncFid,'Z%f F1\n\n',axisZ(1));

        fprintf(ncFid,'( CUTTING BLOCK )\nG93 F[#555]\n');
        for ii = 1:length(axisC)
            fprintf(ncFid,'C%f X%f Z%f\n',axisC(ii),axisX(ii),axisZ(ii));
%             fprintf(ncFid,'%f %f %f\n',axisC(ii),axisX(ii),axisZ(ii));
%             fprintf(ncFid,'GOTO/%f,%f,%f,%f,%f,%f\n', ...
%                 spiralPath1(1,ii),spiralPath1(2,ii),spiralPath1(3,ii), ...
%                 spiralNorm(1,ii),spiralNorm(2,ii),spiralNorm(3,ii));
        end

        fprintf(ncFid,'(linking block)\n');
        fprintf(ncFid,'G94 Z5 F200\n\n');
        fprintf(ncFid,'#5 = #5 + 1\n');
        fprintf(ncFid,'END 2\n\n');
        fprintf(ncFid,'G94 Z20 F800\n');
        fprintf(ncFid,'M29\n');
        fprintf(ncFid,'M30\n');

        fclose(ncFid);
end

%%
% rmpath(genpath('funcs'));