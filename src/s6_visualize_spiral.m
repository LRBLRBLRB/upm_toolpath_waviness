% function visualize_process(toolPathPt,toolNormDirect,toolCutDirect,toolSp, ...
%     toolQuat,toolVec,uLim,peakPt,res)
% to plot the simulation or optimization the spiral path results

while true
    msgfig = questdlg({'Sprial tool path was calculated successfully!', ...
        'Ready for tool path error visualization & machining simulation?'}, ...
        'Spiral Tool Path Simulation','Spiral path error','Machining simulation', ...
        'Save and quit','Save and quit');
    switch msgfig
        case 'Save and quit'
            [spiralPathFileName,spiralPathDirName,spiralPathFileType] = uiputfile({ ...
                '*.mat','MAT-file(*.mat)'; ...
                '*.txt','text-file(.txt)';...
                '*.*','all file(*.*)';...
                }, ...
                'Select the directory and filename to save the surface tool path', ...
                fullfile(workspaceDir,['spiralPath',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
            spiralPathName = fullfile(spiralPathDirName,spiralPathFileName);
            switch spiralPathFileType
                case 0
                    msgfig = msgbox("No spiral tool path saved","Warning","warn","non-modal");
                    uiwait(msgfig);
                case 1
                    Comments = cell2mat(inputdlg( ...
                        'Enter Comment of the spiral tool path:', ...
                        'Saving Comments', ...
                        [5 60], ...
                        string(datestr(now))));
                    save(spiralPathName,"Comments","spiralAngle", ...
                        "spiralPath","spiralQuat","spiralNorm","spiralCut", ...
                        "spiralULim","spiralRes",'spiralPeakPt', ...
                        'spiralInterPtIn','spiralInterPtOut');
            end
            return;
        case 'Spiral path error'
            tResError0 = tic;
            plotNum = 1000;
            spiralResLine = [spiralRes(1,:),spiralRes(2,:)];
            spiralPeakPtLine = [spiralPeakPt(1:3,:),spiralPeakPt(4:6,:)];
            spiralResMaxInd = find(spiralResLine == 5*aimRes);
            spiralResLine(spiralResMaxInd) = [];
            spiralPeakPtLine(:,spiralResMaxInd) = [];
            xPlot = linspace(min(spiralPeakPtLine(1,:)),max(spiralPeakPtLine(1,:)),plotNum);
            yPlot = linspace(min(spiralPeakPtLine(2,:)),max(spiralPeakPtLine(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [spiralResUnique,spiralPeakPtUnique] = groupsummary(spiralResLine',spiralPeakPtLine(1:2,:)',@max);
            resMesh = griddata(spiralPeakPtUnique{1},spiralPeakPtUnique{2},spiralResUnique,xMesh,yMesh);
            figure('Name','Residual height map');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            surf(xMesh,yMesh,resMesh,'EdgeColor','interp'); hold on;
            colormap("parula");
            grid on;
            % plot3(spiralPeakPtUnique{1},spiralPeakPtUnique{2},spiralResUnique,'o', ...
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
            tResError = toc(tResError0);
            fprintf('The time spent in the residual map process is %fs.\n',tResError);
        case 'Machining simulation'
            stepLength = 0.01;
            uLimRound = round(spiralULim,2);
            spiralPathList = [];
            tSimul0 = tic;
            parfor ii = 1:spiralPtNum % each tool path point
                toolSp = toolData.toolBform;
                toolSp.coefs = quat2rotm(spiralQuat(ii,:))*toolCoefs + spiralPath(:,ii);
                tmp = fnval(toolSp,uLimRound(1,ii):stepLength:uLimRound(2,ii));
                % Q{jj} = tmp;
                spiralPathList = [spiralPathList,tmp];
            end
            plotNum = 1000;
            xPlot = linspace(min(spiralPathList(1,:)),max(spiralPathList(1,:)),plotNum);
            yPlot = linspace(min(spiralPathList(2,:)),max(spiralPathList(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [spiralPathZUnique,spiralPathXYUnique] = groupsummary(spiralPathList(3,:)',spiralPathList(1:2,:)',@max);
            zMesh = griddata(spiralPathXYUnique{1},spiralPathXYUnique{2},spiralPathZUnique,xMesh,yMesh);
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