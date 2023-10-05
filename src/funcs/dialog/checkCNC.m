function checkCNC(axisC,axisX,axisZ,textFontSize,textFontType)
%CHECKCNC 此处显示有关此函数的摘要
%   此处显示详细说明
minZ = min(axisZ);
diffC = diff(axisC);
if diffC(1) > 0
    CDir = '+';
else
    CDir = '-';
end
diffX = diff(axisX);
if diffX(1) > 0
    XDir = '+';
else
    XDir = '-';
end
absDiffX = diff(abs(axisX));
isX = ~all(absDiffX(1).*absDiffX > 0);
if isX
    isX = 'true';
else
    isX = 'false';
end

figure('WindowState','maximized');
tiledlayout(1,2);
pos = get(gcf,'position');
set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
nexttile(1);
X = axisX.*cosd(axisC);
Y = axisX.*sind(axisC);
Z = axisZ;
plot3(X,Y,Z,'Color',[0,0.4470,0.7410],'LineStyle',':', ...
    'LineWidth',0.05,'Marker','.','MarkerSize',2);
xlabel('x'); ylabel('y'); zlabel('z');
nexttile(2);
plot(abs(axisX),'LineStyle',':','LineWidth',0.05, ...
    'Marker','.','MarkerSize',2);
ylabel('$\vert{x}\vert$','Interpreter','latex', ...
    'FontSize',textFontSize,'FontName',textFontType, ...
    'FontSmoothing','on','Rotation',0);

opts.WindowStyle = 'non-modal';
opts.Interpreter = 'tex';
warndlg({sprintf(['\\fontsize{%d}\\fontname{%s}', ...
    'Surface was generated successfully!\n'],textFontSize,textFontType), ...
    sprintf('\\fontsize{%d}Minimum of Z-axis is:    \\fontsize{%d}%fmm', ...
    textFontSize,textFontSize + 2,minZ), ...
    sprintf('\\fontsize{%d}Direction of C-axis is:  \\fontsize{%d}%s', ...
    textFontSize,textFontSize + 2,CDir), ...
    sprintf('\\fontsize{%d}Direction of X-axis is:  \\fontsize{%d}%s', ...
    textFontSize,textFontSize + 2,XDir), ...
    sprintf('\\fontsize{%d}Any Wrong loop?          \\fontsize{%d}%s', ...
    textFontSize,textFontSize + 2,isX) ...
    },'CNC code check',opts);
end