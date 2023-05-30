function [fitCenter,fitRadius] = drawCircle(ax21,surfData,options)
%SELECTCIRCLE 此处显示有关此函数的摘要
%   此处显示详细说明

arguments
    ax21
    surfData double
    options.textFontSize uint8 = 14
    options.textFontType char = 'Times New Roman'
    options.unit char = '\mum'
end

% imshow(surfImg,surfImgColorMap);
contourf(ax21,surfData(:,:,1),surfData(:,:,2),surfData(:,:,3),16,'--','LineWidth',0.1); hold on;
colormap(parula(256));
colorbar(ax21,'southoutside');
axis(ax21,'equal');
set(ax21,'FontSize',options.textFontSize,'FontName',options.textFontType);
xLim = get(ax21,'XLim');
yLim = get(ax21,'YLim');
set(ax21,'XLim',sqrt(2)*xLim,'YLim',sqrt(2)*yLim);
xlabel(ax21,['x (',options.unit,')']);
ylabel(ax21,['y (',options.unit,')']);
fprintf(['Please draw a circle to include all the surface: \n' ...
    '(Double-click to return to the main procedure)\n\n']);
roi = drawcircle(ax21,'Color',[0.8500 0.3250 0.0980]);
wait(roi);
fitCenter = roi.Center;
fitRadius = roi.Radius;

end