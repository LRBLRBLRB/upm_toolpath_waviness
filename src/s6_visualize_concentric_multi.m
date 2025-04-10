% function visualize_process(toolPathPt,toolNormDirect,toolCutDirect,toolSp, ...
%     toolQuat,toolVec,uLim,peakPt,res)
% to plot the simulation or optimization results

while true
    msgfig = questdlg({'Residual Height was calculated successfully!', ...
        'Ready for residual height visualization or machining simulation?'}, ...
        'Concentric Tool Path Simulation','Residual height','Concentric machining simulation', ...
        'Save and quit','Save and quit');
    switch msgfig
        case 'Save and quit'
            pathFolderName = getlastfoldername(workspaceDir);
            [pathFileName,pathDirName,pathFileType] = uiputfile({ ...
                '*.mat','MAT-file(*.mat)'; ...
                '*.txt','text-file(.txt)';...
                '*.*','all file(*.*)';...
                }, ...
                'Select the directory and filename to save the surface concentric tool path', ...
                fullfile(workspaceDir,[pathFolderName,'-toolPath',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
            pathName = fullfile(pathDirName,pathFileName);
            switch pathFileType
                case 0
                    msgfig = msgbox("No tool path saved","Warning","warn","non-modal");
                    uiwait(msgfig);
                case 1
                    Comments = cell2mat(inputdlg( ...
                        'Enter Comment of the tool path:', ...
                        'Saving Comments', ...
                        [5 60], ...
                        string(datestr(now))));
                    save(pathName,"Comments", ...
                        "curvePathPt","curveQuat","curveContactU","curvePt","curveRes", ...
                        "curvePeakPt","curveInterPt","curveULim", ...
                        "toolPathAngle","toolPathPt","toolQuat","toolContactU", ...
                        "surfPt","res","peakPt","interPt","uLim", ...
                        "toolNormDirect","toolCutDirect","toolNAccum","accumPtNum","toolRAccum");
            end
            clear msgfig tRes0 tSimul0;
            return;
        case 'Residual height'
            tRes0 = tic;
            plotNum = 1000;
            resLine = [res(1,:),res(2,:)];
            peakPtLine = [peakPt(1:3,:),peakPt(4:6,:)];
            xPlot = linspace(min(peakPtLine(1,:)),max(peakPtLine(1,:)),plotNum);
            yPlot = linspace(min(peakPtLine(2,:)),max(peakPtLine(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [resUnique,peakPtUnique] = groupsummary(resLine',peakPtLine(1:2,:)',@max);
            resMaxInd = find(resUnique == max(resUnique));
            resUnique(resMaxInd) = [];
            peakPtUnique{1}(resMaxInd) = [];
            peakPtUnique{2}(resMaxInd) = [];
            resMesh = griddata(peakPtUnique{1},peakPtUnique{2},resUnique,xMesh,yMesh);
            figure('Name','Residual height map');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            surf(xMesh,yMesh,resMesh,'EdgeColor','interp'); hold on;
            colormap("parula");
            grid on;
            % plot3(peakPtUnique{1},peakPtUnique{2},resUnique,'o', ...
            %     'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
            % cb1 = colorbar;
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['residual height (',unit,')']);
            % legend('residual height in each peakPt','residual height map', ...
            %     'Location','best');
            nexttile;
            contourf(xMesh,yMesh,resMesh,'LineStyle',':'); hold on;
            colormap("turbo");
            axis equal; grid on;
            cb2 = colorbar;
            cb2.Label.String = ['Residual Height (',unit,')'];
            cb2.Layout.Tile = 'east';
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            tRes = toc(tRes0);
            fprintf('The time spent in the residual map process is %fs.\n',tRes);
        case 'Concentric machining simulation'
            stepLength = abs(log10(abs(curvePlotSpar)));
            toolPathList = [];
            tSimul0 = tic;
            toolCoefs = toolData.toolBform.coefs;
            % figure;
            parfor ii = 1:accumPtNum(end) % each tool path point
                toolSp1 = toolData.toolBform;
                toolSp1.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolPathPt(:,ii);
                uLimRound = round(uLim{ii},stepLength);
                for jj = 1:size(uLim{ii},2)
                    tmp = fnval(toolSp1,uLimRound(1,jj):curvePlotSpar:uLimRound(2,jj));
                    % Q{jj} = tmp;
                    toolPathList = [toolPathList,tmp];
                end
            end
            plotNum = 1000;
            xPlot = linspace(min(toolPathList(1,:)),max(toolPathList(1,:)),plotNum);
            yPlot = linspace(min(toolPathList(2,:)),max(toolPathList(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [toolPathZUnique,toolPathXYUnique] = groupsummary(toolPathList(3,:)',toolPathList(1:2,:)',@mean);
            zMesh = griddata(toolPathXYUnique{1},toolPathXYUnique{2},toolPathZUnique,xMesh,yMesh);
            % calculate the error based on the designed surface
            z0Mesh = griddata(surfMesh(:,:,1),surfMesh(:,:,2),surfMesh(:,:,3),xMesh,yMesh);
            % 这里有更好的仿真方式：每个点都计算到设计曲面的距离，而不是沿z方向的距离！！！！
            figure('Name','Machining surface simulation');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            mesh(xMesh,yMesh,zMesh,'EdgeColor','interp'); hold on;
            colormap("jet");
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['z (',unit,')']);
            grid on;
            nexttile;
            mesh(xMesh,yMesh,zMesh - z0Mesh,'EdgeColor','interp'); hold on;
            colormap("jet");
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['z error (',unit,')']);
            grid on;
            tSimul = toc(tSimul0);
            fprintf('The time spent in the simulation calculation process is %fs.\n',tSimul);
        otherwise
            msgfig = msgbox("No tool path saved","Warning","warn","non-modal");
            uiwait(msgfig);
            return;
    end
end