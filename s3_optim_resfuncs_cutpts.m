% to generate the function of residual height 
% related variables: spindle direction, cutting width and arc step.
close all;
clear;
clc;
addpath(genpath('funcs'));
t0 = tic;

% global variables
% global textFontSize textFontType;
workspaceDir = 'workspace/20220925-contrast/nagayama_concentric';
unit = '\mum';
textFontSize = 14;
textFontType = 'Times New Roman';
% profile on
parObj = gcp('nocreate');

[fileName,dirName] = uigetfile({ ...
    '*.mat','MAT-files(*.mat)'; ...
    '*,*','all files(*.*)'}, ...
    'Select one tool edge data file', ...
    fullfile(workspaceDir,'tooltheo.mat'), ...
    'MultiSelect','off');
toolName = fullfile(dirName,fileName);
% toolName = 'output_data\tool\toolTheo_3D.mat';
toolData = load(toolName);


%% the influence of cut width
tool2D = toolData;
tool2D.toolBform.coefs(1,:) = [];
tool2D.toolBform.dim = 2;
r = tool2D.radius;
tool2D.center = [0;0];
tool2D.cutDirect = [];
tool2D.toolDirect = [1;0];
tool2D.toolEdgeNorm = [0;1];
tool2D.toolEdgePt = [];
tool2D.toolFit = [];

tWidth0 = tic;
% widthAvi = VideoWriter("debug/videos/cutWidth.mp4",'MPEG-4');
% widthAvi.FrameRate = 60;
% open(widthAvi);
% figure('Name','Video of cut width - residual height funciton');

travNum = 1000;
cutWidth = linspace(0,r,travNum);
slope = linspace(-1,1,travNum);
vec1 = [0;-1];
vec2 = [0;-1];
c1 = [0;0];
s1 = tool2D;
[quat1,trans1,u1] = toolPos(tool2D,c1,related to slope,vec1);
s1.coefs = quat2rotm(quat1)*s1.coefs + trans1;
res = zeros(travNum,travNum);
parfor ii = 1:travNum*travNum
    kk = ceil(ii/travNum);
    jj = ii - (kk - 1)*travNum;
    c2 = [cutWidth(kk);cutWidth(kk)*slope(jj)];
    s2 = tool2D;

    [quat2,trans2,u2,isCollision] = toolPos(tool2D,c2,related to slope,vec2);
    if isCollision == false
        s2.coefs = quat2rotm(quat2)*s2.coefs + trans2(1:2);
    end
    [res(kk,jj),~] = residual2D_numeric(s1,s2,1e-3, ...
        c1,c2,vec1,vec2,r,'DSearchn');
%     scatter(c1,0,'MarkerEdgeColor',[0,0.4470,0.7410]); hold on;
%     scatter(c2(ii),0,'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
%     scatter(peakPt(1),peakPt(2),'MarkerFaceColor',[0.4940,0.1840,0.5560]);
%     pt1 = fnval(s1,0:1e-3:1);
%     % pt1 = fnval(s1,ind1:1e-3:1);
%     % pt1Shade = fnval(s1,0:1e-3:ind1);
%     pt2 = fnval(s2,0:1e-3:1);
%     % pt2 = fnval(s2,0:1e-3:ind2);
%     % pt2Shade = fnval(s2,ind2:1e-3:1);
%     plot(pt1(1,:),pt1(2,:),'Color',[0,0.4470,0.7410]);
%     % plot(pt1Shade(1,:),pt1Shade(2,:),'Color',[0,0.4470,0.7410]*0.5);
%     plot(pt2(1,:),pt2(2,:),'Color',[0.8500,0.3250,0.0980]);
%     % plot(pt2Shade(1,:),pt2Shade(2,:),'Color',[0.8500,0.3250,0.0980]*0.5);
%     grid on;
%     axis equal;
%     xlabel(['x (',unit,')']);
%     ylabel(['y (',unit,')']);
%     set(gca,'FontSize',textFontSize,'FontName',textFontType);
%     hold off;
%     drawnow;
%     frame = getframe(gcf);
%     writeVideo(widthAvi,frame);
%     disp(ii);
end
% close(widthAvi);

figure('Name','Residual height function');
[widthRes(:,:,1),widthRes(:,:,2)] = meshgrid(cutWidth,slope);
widthRes(:,:,3) = res;
surf(widthRes(:,:,1),widthRes(:,:,2),widthRes(:,:,3),'EdgeColor','none');
xlabel(['Cut Width (',unit,')']);
ylabel(['Residual Height (',unit,')']);
widthResPath = fullfile(workspaceDir,'widthRes.mat');
save(widthResPath,"widthRes");
tWidth = toc(tWidth0);
fprintf("The time spent in the spindle-direction-traversing process is %fs.\n",tWidth);



%%
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% delete(parObj);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));