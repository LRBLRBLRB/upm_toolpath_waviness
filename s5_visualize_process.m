while true
    msgfig = questdlg({'Residual Height was calculated successfully!', ...
        'Ready for residual height visualization or machining simulation?'}, ...
        'Tool Path Simulation','Residual height','Machining simulation', ...
        'Save and quit','Save and quit');
    switch msgfig
        case 'Save and quit'
            [pathFileName,pathDirName,pathFileType] = uiputfile({ ...
                '*.mat','MAT-file(*.mat)'; ...
                '*.txt','text-file(.txt)';...
                '*.*','all file(*.*)';...
                }, ...
                'Select the directory and filename to save the surface tool path', ...
                fullfile(workspaceDir,['toolPath',datestr(now,'yyyymmddTHHMMSS'),'.mat']));
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
                    save(pathName,"Comments","toolPathPt","toolNormDirect","toolCutDirect", ...
                        "toolSp","toolQuat","toolVec",'uLim',"res");
            end
            % delete(parObj);
            profile off
            tTol = toc(t0);
            fprintf("The time spent in the whole process is %fs.\n",tTol);
            % profile viewer
            % profsave(profile("info"),"profile_data");
            % rmpath(genpath('funcs'));
            return; 
        case 'Residual height'
            tRes0 = tic;
            plotNum = 1000;
            xPlot = linspace(min(peakPt(1,:)),max(peakPt(1,:)),plotNum);
            yPlot = linspace(min(peakPt(2,:)),max(peakPt(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [resUnique,peakPtUnique] = groupsummary(res',peakPt(1:2,:)',@max);
            resMesh = griddata(peakPtUnique{1},peakPtUnique{2},resUnique,xMesh,yMesh);
            figure('Name','Residual height map');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1)+pos(4)/2-pos(4),pos(2),2*pos(3),pos(4)]);
            tiledlayout(1,2);
            nexttile;
            surf(xMesh,yMesh,resMesh,'EdgeColor','interp'); hold on;
            colormap("parula");
            grid on;
            plot3(peakPt(1,:),peakPt(2,:),res,'o', ...
                'MarkerEdgeColor',[0.8500,0.3250,0.0980]);
            % cb1 = colorbar;
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            zlabel(['residual height (',unit,')']);
            legend('residual height in each peakPt','residual height map', ...
                'Location','best');
            nexttile;
            contourf(xMesh,yMesh,resMesh); hold on;
            colormap("parula");
            axis equal; grid on;
            cb2 = colorbar;
            cb2.Label.String = ['Residual Height (',unit,')'];
            cb2.Layout.Tile = 'east';
            set(gca,'FontSize',textFontSize,'FontName',textFontType);
            xlabel(['x (',unit,')']);
            ylabel(['y (',unit,')']);
            tRes = toc(tRes0);
            fprintf('The time spent in the residual map process is %fs.\n',tRes);
        case 'Machining simulation'
            stepLength = 0.01;
            nLoop = ceil(ptNum/sparTheta);
            uLimRound = round(uLim,2);
            toolPathList = [];
            tSimul0 = tic;
            for ii = 1:nLoop % each loop
                Q = cell(sparTheta,1);
                for jj = 1:sparTheta
                    toolSp = toolData.toolBform;
                    toolSp.coefs = quat2rotm(toolQuat((ii-1)*nLoop+jj,:))*toolCoefs + toolVec(:,(ii-1)*nLoop+jj);
                    tmp = fnval(toolSp,uLimRound(1,(ii-1)*nLoop+jj):stepLength:uLimRound(2,(ii-1)*nLoop+jj));
                    % Q{jj} = tmp;
                    toolPathList = [toolPathList,tmp];
                end
%                 for u = 0:stepLength:1
%                     for jj = 1:sparTheta
%                         if u >= uLimRound(1,(ii-1)*nLoop+jj) && u <= uLimRound(2,(ii-1)*nLoop+jj)
%                             tmp = Q{jj}(:,round((u - uLimRound(1,(ii-1)*nLoop+jj))/stepLength + 1));
%                             toolPathList = [toolPathList,tmp];
%                         end
%                     end
%                 end
            end
            % plot3(toolPathList(1,:),toolPathList(2,:),toolPathList(3,:),'.','Color',[0,0.4450,0.7410]);
            plotNum = 1000;
            xPlot = linspace(min(toolPathList(1,:)),max(toolPathList(1,:)),plotNum);
            yPlot = linspace(min(toolPathList(2,:)),max(toolPathList(2,:)),plotNum);
            [xMesh,yMesh] = meshgrid(xPlot,yPlot);
            % elliminate the smaller residual height at the same peak
            [toolPathZUnique,toolPathXYUnique] = groupsummary(toolPathList(3,:)',toolPathList(1:2,:)',@max);
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
    end
end