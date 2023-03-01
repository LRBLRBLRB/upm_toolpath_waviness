% function visualize_process(toolPathPt,toolNormDirect,toolCutDirect,toolSp, ...
%     toolQuat,toolVec,uLim,peakPt,res)
% to plot the simulation or optimization results

while true
    msgfig = questdlg({'Sprial tool path was calculated successfully!', ...
        'Ready for tool path error visualization & machining simulation?'}, ...
        'Tool Path Simulation','Tool path error','Machining simulation', ...
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
                        "toolSp","toolQuat",'uLim',"res");
            end
            return;
        case 'Tool path error'
            
        case 'Machining simulation'
            stepLength = 0.01;
            uLimRound = round(uLim,2);
            toolPathList = [];
            tSimul0 = tic;
            for ii = accumPtNum(1):accumPtNum(end) % each tool path point
                toolSp = toolData.toolBform;
                toolSp.coefs = quat2rotm(toolQuat(ii,:))*toolCoefs + toolVec(:,ii);
                tmp = fnval(toolSp,uLimRound(1,ii):stepLength:uLimRound(2,ii));
                % Q{jj} = tmp;
                toolPathList = [toolPathList,tmp];
            end
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
        otherwise
            msgfig = msgbox("No tool path saved","Warning","warn","non-modal");
            uiwait(msgfig);
            return;
    end
end