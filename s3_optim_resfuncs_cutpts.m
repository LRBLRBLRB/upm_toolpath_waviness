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
res = zeros(1,travNum);
parfor ii = 1:travNum
    s2 = tool2D;
    s2.coefs = R2*s2.coefs + [c2(ii);0];
    [toolQuat(ii,:),toolVec(:,ii),toolContactU(ii),isCollision(ii)] = toolPos( ...
        toolData,surfPt(:,ii),surfNorm(:,ii),surfDirect(:,ii),[0;0;1]);
    if isCollision(ii) == false
        toolPathPt(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.center + toolVec(:,ii);
        toolCutDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.cutDirect;
        toolNormDirect(:,ii) = quat2rotm(toolQuat(ii,:))*toolData.toolEdgeNorm;
    end
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
% tool2D = toolData.toolBform;
% tool2D.coefs(1,:) = [];
% tool2D.dim = 2;
% r = toolData.radius;
% openAngle = toolData.openAngle;
% 
% tSpin0 = tic;
% % spindleAvi = VideoWriter("debug/videos/spindleDir.mp4",'MPEG-4');
% % spindleAvi.FrameRate = 60;
% % open(spindleAvi);
% % figure('Name','Video of Spindle Direction - Residual Height Funciton');
% 
% travNum = 1000;
% p1 = [0;0];
% p2 = [r*ones(1,travNum);zeros(1,travNum)];
% res = zeros(1,travNum);
% theta = linspace(pi - openAngle*pi/180/2,pi + openAngle*pi/180/2,travNum);
% for ii = 1:travNum
%     s1 = tool2D;
%     R1 = rotz(theta(ii));
%     s1.coefs = R1(1:2,1:2)*s1.coefs;
%     s2 = tool2D;
%     R2 = rotz(theta(ii));
%     s2.coefs = R2(1:2,1:2)*s2.coefs + [p2(ii);0];
%     vec2 = R2*vec2Ref;
%     [res(ii),peakPt,ind1,ind2] = residual2D_numeric(s1,s2,1e-3,p1,p2,'DSearchn');
% %     scatter(c1,0,'MarkerEdgeColor',[0,0.4470,0.7410]); hold on;
% %     scatter(c2(ii),0,'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
% %     scatter(peakPt(1),peakPt(2),'MarkerFaceColor',[0.4940,0.1840,0.5560]);
% %     pt1 = fnval(s1,0:1e-3:1);
% %     % pt1 = fnval(s1,ind1:1e-3:1);
% %     % pt1Shade = fnval(s1,0:1e-3:ind1);
% %     pt2 = fnval(s2,0:1e-3:1);
% %     % pt2 = fnval(s2,0:1e-3:ind2);
% %     % pt2Shade = fnval(s2,ind2:1e-3:1);
% %     plot(pt1(1,:),pt1(2,:),'Color',[0,0.4470,0.7410]);
% %     % plot(pt1Shade(1,:),pt1Shade(2,:),'Color',[0,0.4470,0.7410]*0.5);
% %     plot(pt2(1,:),pt2(2,:),'Color',[0.8500,0.3250,0.0980]);
% %     % plot(pt2Shade(1,:),pt2Shade(2,:),'Color',[0.8500,0.3250,0.0980]*0.5);
% %     grid on;
% %     axis equal;
% %     xlabel(['x (',unit,')']);
% %     ylabel(['y (',unit,')']);
% %     set(gca,'FontSize',textFontSize,'FontName',textFontType);
% %     hold off;
% %     drawnow;
% %     frame = getframe(gcf);
% %     writeVideo(spindleAvi,frame);
% %     disp(ii);
% end
% % close(spindleAvi);
% 
% figure('Name','Spindel Direction - Residual Height Function');
% plot((theta - pi)/pi*180,res);
% xlabel(['Spindle Angle (',unit,')']);
% ylabel(['Residual Height (',unit,')']);
% spinRes(1,:) = theta - pi;
% spinRes(2,:) = res;
% spinResPath = fullfile(workspaceDir,'spinRes.mat');
% save(spinResPath,"spinRes");
% tSpin = toc(tSpin0);
% fprintf("The time spent in the spindle-direction-traversing process is %fs.\n",tSpin);

%% the influence of both spindle direction and cut width
tool2D = toolData.toolBform;
tool2D.coefs(1,:) = [];
tool2D.dim = 2;
r = toolData.radius;
openAngle = toolData.openAngle;

tBoth0 = tic;
% bothAvi = VideoWriter("debug/videos/both.mp4",'MPEG-4');
% bothAvi.FrameRate = 60;
% open(bothAvi);
% figure('Name','Video of Residual Height Funciton');

travNum = 100;
p1 = [0;0];
p2 = [linspace(0.5*r,2*r,travNum);zeros(1,travNum)];
res = zeros(travNum,travNum,travNum);
theta2 = linspace(-openAngle*pi/180/2,openAngle*pi/180/2,travNum);
phi = linspace(pi-pi,pi+pi,travNum);

for kk = 1:travNum
    s1 = tool2D;
    R1 = rotz(phi(kk));
    s1.coefs = R1(1:2,1:2)*s1.coefs + p1;
    for ii = 1:travNum
        for jj = 1:travNum
            s2 = tool2D;
            R2 = rotz(phi(kk) + theta2(jj));
            R2 = R2(1:2,1:2);
            s2.coefs = R2*s2.coefs + [c2(ii);0];
            vec2 = R2*vec2Ref;
            [res(ii,jj,kk),peakPt,ind1,ind2] = residual2D_numeric(s1,s2,1e-3, ...
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
end
% close(bothAvi);

figure('Name','Residual height function');
[xMesh,yMesh] = meshgrid(c2,theta2);
surf(xMesh,(yMesh - pi)/pi*180,res,'FaceColor','interp'); hold on;
cb = colorbar;
set(gca,'FontSize',textFontSize,'FontName',textFontType);
xlabel(['Cut Width (',unit,')']);
ylabel('Spindle Angle (\circ)');
zlabel(['Residual Height (',unit,')']);
bothRes(:,:,1) = xMesh;
bothRes(:,:,2) = yMesh - pi;
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