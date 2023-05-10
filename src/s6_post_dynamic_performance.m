% test whether the dynamic performances is good or not

isAPP = false;
if isAPP
else
    workspaceDir = uigetdir('../workspace/20230504 D906', ...
        'Select the Workspace Directory');
    if workspaceDir
        msgfig = msgbox('New tool path will be selected!','Exit','help','non-modal');
        uiwait(msgfig);
        [fileName,dirName] = uigetfile({ ...
            '*.mat','MAT-files(*.mat)'; ...
            '*,*','All Files(*.*)'}, ...
            'Select one tool path data file', ...
            fullfile(workspaceDir,'spiralpath.mat'), ...
            'MultiSelect','off');
        toolPathPath = fullfile(dirName,fileName);
        load(toolPathPath);
    end
end

%% concentric acceleration check
acce = diff(diff(curvePathPt,1,2),1,2);
maxAcce = max(acce);
figure;
plot(curvePathPt(1,2:end - 1),acce(3,:)./acce(1,:));
xlabel(['r (',unit,')']);
ylabel(['a (',unit,'^{-1})']);

%% spiral check
spiralZ = spiralPath(3,:);
spiralR = vecnorm(spiralPath(1:2,:),2,1);

dz = diff(spiralZ);
dr = diff(spiralR);
max(dz)
max(dr)

figure;
tiledlayout(2,1);
nexttile;
yyaxis left;
plot(spiralAngle(1:end),spiralZ,'.-','Color',[0,0.4470,0.7410]);
hold on;
ylabel('z');
yyaxis right;
bar(spiralAngle(2:end),dz,'grouped','EdgeColor','none','FaceColor',"flat",'FaceAlpha',0.3);
nexttile;
yyaxis left;
plot(spiralAngle(1:end),spiralR,'.-');
ylabel('r');
yyaxis right;
bar(spiralAngle(2:end),dr,'grouped','EdgeColor','none','FaceColor',"flat",'FaceAlpha',0.3);
xlabel('\theta');

ddz = diff(diff(spiralZ));
ddr = diff(diff(spiralR));

figure;
tiledlayout(2,1);
nexttile;
plot(ddz);
hold on;
ylabel(['a_{z} (',unit,'/s^2)']);
nexttile;
plot(ddr);
ylabel(['a_{r} (',unit,'/s^2)']);