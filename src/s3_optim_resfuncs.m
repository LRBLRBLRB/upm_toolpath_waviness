% to generate the function of residual height 
% related variables: spindle direction, cutting width and arc step.
close all;
clear;
clc;
addpath(genpath('funcs'));
t0 = tic;

% global variables
% global textFontSize textFontType;
workspaceDir = fullfile('..','workspace','/20220925-contrast/nagayama_concentric';
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
tool2D = toolData.toolBform;
tool2D.coefs(1,:) = [];
tool2D.dim = 2;
r = toolData.radius;

tWidth0 = tic;
% widthAvi = VideoWriter("debug/videos/cutWidth.mp4",'MPEG-4');
% widthAvi.FrameRate = 60;
% open(widthAvi);
% figure('Name','Video of cut width - residual height funciton');

travNum = 1000;
c1 = 0;
c2 = linspace(0.5*r,2*r,travNum);
vec1 = [0;-1];
vec2 = [0;-1];
R1 = rotz(pi);
s1 = tool2D;
s1.coefs = R1(1:2,1:2)*s1.coefs;
R2 = rotz(pi);
R2 = R2(1:2,1:2);
res = zeros(1,travNum);
parfor ii = 1:travNum
    s2 = tool2D;
    s2.coefs = R2*s2.coefs + [c2(ii);0];
    [res(ii),peakPt,ind1,ind2] = residual2D_numeric(s1,s2,1e-3, ...
        [c1;0],[c2(ii);0],vec1,vec2,r,'DSearchn');
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

figure('Name','Cut width - residual height function');
plot(c2(1,:),res);
xlabel(['Cut Width (',unit,')']);
ylabel(['Residual Height (',unit,')']);
widthRes(1,:) = c2(1,:);
widthRes(2,:) = res;
widthResPath = fullfile(workspaceDir,'widthRes.mat');
save(widthResPath,"widthRes");
tWidth = toc(tWidth0);
fprintf("The time spent in the spindle-direction-traversing process is %fs.\n",tWidth);

%% the influence of spindle direction
tool2D = toolData.toolBform;
tool2D.coefs(1,:) = [];
tool2D.dim = 2;
r = toolData.radius;

tSpin0 = tic;
% spindleAvi = VideoWriter("debug/videos/spindleDir.mp4",'MPEG-4');
% spindleAvi.FrameRate = 60;
% open(spindleAvi);
% figure('Name','Video of Spindle Direction - Residual Height Funciton');

travNum = 1000;
c1 = 0;
c2 = r*ones(1,travNum);
vec1 = [0;-1];
vec2Ref = [0;-1];
R1 = rotz(pi);
s1 = tool2D;
s1.coefs = R1(1:2,1:2)*s1.coefs;
res = zeros(1,travNum);
theta = linspace(pi - 10*pi/180,pi + 10*pi/180,travNum);
for ii = 1:travNum
    s2 = tool2D;
    R2 = rotz(theta(ii));
    R2 = R2(1:2,1:2);
    s2.coefs = R2*s2.coefs + [c2(ii);0];
    vec2 = R2*vec2Ref;
    [res(ii),peakPt,ind1,ind2] = residual2D_numeric(s1,s2,1e-3, ...
        [c1;0],[c2(ii);0],vec1,vec2,r,'DSearchn');
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
%     writeVideo(spindleAvi,frame);
%     disp(ii);
end
% close(spindleAvi);

figure('Name','Spindel Direction - Residual Height Function');
plot((theta - pi)/pi*180,res);
xlabel(['Spindle Angle (',unit,')']);
ylabel(['Residual Height (',unit,')']);
spinRes(1,:) = theta - pi;
spinRes(2,:) = res;
spinResPath = fullfile(workspaceDir,'spinRes.mat');
save(spinResPath,"spinRes");
tSpin = toc(tSpin0);
fprintf("The time spent in the spindle-direction-traversing process is %fs.\n",tSpin);

%% the influence of both spindle direction and cut width
tool2D = toolData.toolBform;
tool2D.coefs(1,:) = [];
tool2D.dim = 2;
r = toolData.radius;

tBoth0 = tic;
% bothAvi = VideoWriter("debug/videos/both.mp4",'MPEG-4');
% bothAvi.FrameRate = 60;
% open(bothAvi);
% figure('Name','Video of Residual Height Funciton');

travNum = 100;
c1 = 0;
c2 = linspace(0.5*r,2*r,travNum);
vec1 = [0;-1];
vec2Ref = [0;-1];
R1 = rotz(pi);
s1 = tool2D;
s1.coefs = R1(1:2,1:2)*s1.coefs;
res = zeros(travNum,travNum);
theta = linspace(pi - 10*pi/180,pi + 10*pi/180,travNum);
for ii = 1:travNum
    for jj = 1:travNum
%         ii = 80;
%         jj = 100;
        s2 = tool2D;
        R2 = rotz(theta(jj));
        R2 = R2(1:2,1:2);
        s2.coefs = R2*s2.coefs + [c2(ii);0];
        vec2 = R2*vec2Ref;
        [res(ii,jj),peakPt,ind1,ind2] = residual2D_numeric(s1,s2,1e-3, ...
            [c1;0],[c2(ii);0],vec1,vec2,r,'DSearchn');
%         scatter(c1,0,'MarkerEdgeColor',[0,0.4470,0.7410]); hold on;
%         scatter(c2(ii),0,'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
%         scatter(peakPt(1),peakPt(2),'MarkerFaceColor',[0.4940,0.1840,0.5560]);
%         pt1 = fnval(s1,0:1e-3:1);
%         % pt1 = fnval(s1,ind1:1e-3:1);
%         % pt1Shade = fnval(s1,0:1e-3:ind1);
%         pt2 = fnval(s2,0:1e-3:1);
%         % pt2 = fnval(s2,0:1e-3:ind2);
%         % pt2Shade = fnval(s2,ind2:1e-3:1);
%         plot(pt1(1,:),pt1(2,:),'Color',[0,0.4470,0.7410]);
%         % plot(pt1Shade(1,:),pt1Shade(2,:),'Color',[0,0.4470,0.7410]*0.5);
%         plot(pt2(1,:),pt2(2,:),'Color',[0.8500,0.3250,0.0980]);
%         % plot(pt2Shade(1,:),pt2Shade(2,:),'Color',[0.8500,0.3250,0.0980]*0.5);
%         grid on;
%         axis equal;
%         xlabel(['x (',unit,')']);
%         ylabel(['y (',unit,')']);
%         set(gca,'FontSize',textFontSize,'FontName',textFontType);
%         hold off;
%         drawnow;
%         frame = getframe(gcf);
%         writeVideo(bothAvi,frame);
%         disp((ii - 1)*travNum + jj);
    end
end
% close(bothAvi);

figure('Name','Residual height function');
[thetaMesh,c2Mesh] = meshgrid(theta,c2);
surf(c2Mesh,(thetaMesh - pi)/pi*180,res, ...
    'FaceColor','interp','EdgeColor','none'); hold on;
cb = colorbar;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['Cut Width (',unit,')']);
ylabel('Spindle Angle (\circ)');
zlabel(['Residual Height (',unit,')']);
bothRes(:,:,1) = c2Mesh;
bothRes(:,:,2) = thetaMesh - pi;
bothRes(:,:,3) = res;
bothResPath = fullfile(workspaceDir,'bothRes.mat');
save(bothResPath,"bothRes");
tBoth = toc(tBoth0);
fprintf("The time spent in the both-traversing process is %fs.\n",tBoth);

%%
tTol = toc(t0);
fprintf("The time spent in the whole process is %fs.\n",tTol);
% delete(parObj);
% profile viewer
% profsave(profile("info"),"profile_data");
% rmpath(genpath('funcs'));