function [ax23,outputArg2] = selectCircle(varargin)
%SELECTCIRCLE 此处显示有关此函数的摘要
%   此处显示详细说明

if isempty(varargin)
    fig2 = figure('Name','2 surface located','WindowState','maximized');
    tiled2 = tiledlayout(2,2,'TileSpacing','compact');
    ax21 = nexttile(1,[2,1]);
end

questOpt.Interpreter = 'tex';
questOpt.Default = 'Re-select';

msgfig2 = 'Re-select';
while strcmp(msgfig2,'Re-select')
    if strcmp(msgfig1,'Re-select')
        hold(ax21,'off');
        % imshow(surfImg,surfImgColorMap);
        contourf(ax21,surfData0(:,:,1),surfData0(:,:,2),surfData0(:,:,3),16,'--','LineWidth',0.1); hold on;
        colormap(parula(256));
        colorbar(ax21,'southoutside');
        axis(ax21,'equal');
        set(ax21,'FontSize',textFontSize,'FontName',textFontType);
        xLim = get(ax21,'XLim');
        yLim = get(ax21,'YLim');
        set(ax21,'XLim',sqrt(2)*xLim,'YLim',sqrt(2)*yLim);
        xlabel(ax21,['x (',unit,')']);
        ylabel(ax21,['y (',unit,')']);
        fprintf(['Please draw a circle to include all the surface: \n' ...
            '(Double-click to return to the main procedure)\n\n']);
        roi = drawcircle(ax21,'Color',[0.8500 0.3250 0.0980]);
        wait(roi);
        fitCenter = roi.Center;
        fitRadius = roi.Radius;
    end

    % to get the indices of the plane
    % isSurf = inROI(roi,surfData0(:,:,1),surfData0(:,:,2));
    isPlane = (surfData0(:,:,1) - fitCenter(1)).^2 + (surfData0(:,:,2) - fitCenter(1)).^2 > fitRadius.^2;
    % plane extraction
    planeMesh1 = surfData0;
    planeMesh1(~(repmat(isPlane,[1,1,3]))) = nan;
    ax22 = nexttile(2);
    hold(ax22,'off');
    surf(ax22,planeMesh1(:,:,1),planeMesh1(:,:,2),planeMesh1(:,:,3), ...
        'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
    colormap(parula(256));
    colorbar(ax22,'eastoutside');
    clim([min(planeMesh1(:,:,3),[],'all'),max(planeMesh1(:,:,3),[],'all')]);
    xlabel(ax22,['x (',unit,')']);
    ylabel(ax22,['y (',unit,')']);
    zlabel(ax22,['z (',unit,')']);
    % plot3(planeData(:,1),planeData(:,2),planeData(:,3),'.','Color',[0,0.4450,0.7410]);
    % pause();
    msgfig2 = questdlg({sprintf(['\\fontsize{%d}\\fontname{%s} ', ...
        'Surface selection finished successfully!'],textFontSize,textFontType), ...
        'Re-selet or not?'}, ...
        'Concentric tool path Generation','Re-select','OK',questOpt);
end

end

