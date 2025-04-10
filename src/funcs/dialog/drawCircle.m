function [fitCenter,fitRadius] = drawCircle(ax,surfData,options)
%SELECTCIRCLE 此处显示有关此函数的摘要
%   此处显示详细说明

arguments
    ax
    surfData double
    options.textFontSize uint8 = 14
    options.textFontType char = 'Times New Roman'
    options.unit char = '\mum'
    options.xlim = []
    options.ylim = []
    options.level = 16
end

% imshow(surfImg,surfImgColorMap);
contourf(ax,surfData(:,:,1),surfData(:,:,2),surfData(:,:,3),options.level, ...
    '--','LineWidth',0.1); hold on;
colormap(ax,jet(256));
colorbar(ax,'southoutside');
axis(ax,'equal');
set(ax,'FontSize',options.textFontSize,'FontName',options.textFontType);
if isempty(options.xlim)
    options.xlim = sqrt(2)*get(ax,'XLim');
    options.ylim = sqrt(2)*get(ax,'YLim');
end
set(ax,'XLim',options.xlim,'YLim',options.ylim);
xlabel(ax,['x (',options.unit,')']);
ylabel(ax,['y (',options.unit,')']);
fprintf(['Please draw a circle to include all the surface: \n' ...
    '(Double-click to return to the main procedure)\n\n']);
roi = drawcircle(ax,'Color',[0.8500 0.3250 0.0980]);
wait(roi);
fitCenter = roi.Center;
fitRadius = roi.Radius;

end