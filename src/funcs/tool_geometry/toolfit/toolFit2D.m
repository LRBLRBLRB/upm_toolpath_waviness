function [circ2D,scatterDst,RMSE,varargout] = toolFit2D(scatterOri, ...
    arcRansacMaxDist,lineFitMaxDist,options)
% usage: [circ2D,scatterDst,RMSE] = toolFit2D(scatterOri)
%
% solve the edge sharpness of a arc turning tool: 
% calculate the radius of the tool tip arc, 
% as well as transforming the tool coordinate
%
% Inputs: 
%   scatterMat: original tool edge point from measuring (2,n)
%   options: 
%       fitMethod: the method of curve fitting
% Outputs: 
%   circ2D: the struct that includes all below
%       c: center of the arc (2,1)
%       r: radius of the arc (1,1)
%       ang: the open angle of the tool (1,1)
%   scatterDst: tool edge points with pose adjustment (2,n)
%   RMSE: the square-mean-root error of curve fitting
% method:
%   least square method by normal equation solving

arguments
    scatterOri (2,:) {mustBeFinite}
    arcRansacMaxDist {mustBePositive} = 1
    lineFitMaxDist {mustBePositive} = 1
    options.toolFitType {mustBeMember(options.toolFitType, ...
        {'onlyArc','arcRansac','lineArc'})} = 'onlyArc'
    options.arcFitMethod {mustBeMember(options.arcFitMethod, ...
        {'gradient-decent','normal-equation','levenberg-marquardt'})} ...
        = 'levenberg-marquardt'
    options.lineFitMethod {mustBeMember(options.lineFitMethod, ...
        {'polyfit','ransac'})} = 'polyfit'
    options.arcFitdisplayType {mustBeMember(options.arcFitdisplayType, ...
        {'off','none','iter','iter-detailed','final','final-detailed'})} = 'final'
end

%% tool tip fitting
switch options.toolFitType
    case 'onlyArc' 
        % ------------tool tip remains only the arc, without line edges------------
    case 'arcRansac'
        % ------------tool tip with line edges, and use ransac to fit the arc directly------------
        % ransac
        sampleSz = 3; % number of points to sample per trial
        % maxDist = 0.003; % max allowable distance for inliers
        
        fitLineFcn = @(pts) arcFit2D(pts','displayType','off');  % fit function
        evalLineFcn = ...   % distance evaluation function
          @(mdl, pts) abs(vecnorm(pts - (mdl.center)',2,2) - mdl.radius);
                
        [~,inlierInd] = ransac(scatterOri',fitLineFcn,evalLineFcn, ...
          sampleSz,arcRansacMaxDist);
        scatterOri = scatterOri(:,inlierInd);
        varargout{1} = [];
    case 'lineArc'
        % get the start and end point of the line
        fprintf(['Please select the points of the two edges of the tool tip: \n' ...
            'The first two is the interval of the left one, ' ...
            'while the other two is that of the right.\n']);
        % we have three methods to get the points
        % [lineX,lineY] = ginput(2);
        % [lineX,lineY] = getpts(fig1);
        % ax1.Children

        isContinue = 0;
        while ~isContinue
            lineMid = ginput(4);
            fprintf('\nSuccessfully selected!\n\n')
            leftPts = scatterOri(:,(scatterOri(1,:) > lineMid(1,1)) & (scatterOri(1,:) < lineMid(2,1)));
            rightPts = scatterOri(:,(scatterOri(1,:) > lineMid(3,1)) & (scatterOri(1,:) < lineMid(4,1)));
            if strcmp(options.lineFitMethod,'ransac')
                % ------------ransac line fitting to remove outlieres & lsc arc fitting------------
                % ransac line fitting
                sampleSz = 2; % number of points to sample per trial
                % maxDist = 0.00001; % max allowable distance for inliers
                fitLineFcn = @(pts) polyfit(pts(:,1),pts(:,2),1); % fit function using polyfit
                evalLineFcn = ...   % distance evaluation function
                  @(mdl, pts) sum((pts(:, 2) - polyval(mdl, pts(:,1))).^2,2);
    
                [leftPoly,leftInlierIdx] = ransac(leftPts',fitLineFcn,evalLineFcn, ...
                  sampleSz,lineFitMaxDist);
                [rightPoly,rightInlierIdx] = ransac(rightPts',fitLineFcn,evalLineFcn, ...
                  sampleSz,lineFitMaxDist);
    
                % leftDist = pt2Line(scatterOri,leftPoly);
                % leftInlierIdx = leftDist < lineFitMaxDist;
                % rightDist = pt2Line(scatterOri,rightPoly);
                % rightInlierIdx = rightDist < lineFitMaxDist;
                % plot the line fitting process
                figure('Name','Line Fitting of the Tool Tip Arc');
                plot(leftPts(1,leftInlierIdx),leftPoly(1)*leftPts(1,leftInlierIdx) + leftPoly(2), ...
                    '.','Color',[0.8500    0.3250    0.0980],'MarkerSize',8); hold on;
                plot(rightPts(1,rightInlierIdx),rightPoly(1)*rightPts(1,rightInlierIdx) + rightPoly(2), ...
                    '.','Color',[0.9290    0.6940    0.1250],'MarkerSize',8);
                plot(scatterOri(1,:),scatterOri(2,:),'LineWidth',0.5,'Color',[0    0.4470    0.7410]);
                yLim = get(gca,'YLim');
                line([lineMid(1,1),lineMid(1,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(2,1),lineMid(2,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(3,1),lineMid(3,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(4,1),lineMid(4,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                grid on; axis equal;

                % circle boundary
                circIdx(1) = find(leftInlierIdx,1,"last");
                circIdx(2) = size(scatterOri,2) - find(flipud(rightInlierIdx),1,"last");
            else
                % ------------least-square line fitting to remove outliers & lsc arc fitting------------
                % least square line fitting 
                leftPoly = polyfit(leftPts(1,:),leftPts(2,:),1);
                leftDist = pt2Line(scatterOri,leftPoly);
                leftInlierIdx = leftDist < lineFitMaxDist;
    
                rightPoly = polyfit(rightPts(1,:),rightPts(2,:),1);
                rightDist = pt2Line(scatterOri,rightPoly);
                rightInlierIdx = rightDist < lineFitMaxDist;
    
                % plot the line fitting process
                figShow = figure('Name','Line Fitting of the Tool Tip Arc');
                plot(scatterOri(1,leftInlierIdx),leftPoly(1)*scatterOri(1,leftInlierIdx) + leftPoly(2), ...
                    '.','Color',[0.8500    0.3250    0.0980],'MarkerSize',8); hold on;
                plot(scatterOri(1,rightInlierIdx),rightPoly(1)*scatterOri(1,rightInlierIdx) + rightPoly(2), ...
                    '.','Color',[0.929,0.694,0.1250],'MarkerSize',8);
                plot(scatterOri(1,:),scatterOri(2,:),'.','MarkerSize',1, ...
                    'LineWidth',0.5,'Color',[0    0.4470    0.7410]);
                yLim = get(gca,'YLim');
                line([lineMid(1,1),lineMid(1,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(2,1),lineMid(2,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(3,1),lineMid(3,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                line([lineMid(4,1),lineMid(4,1)],[yLim(1),yLim(2)],'Color',[0.4940    0.1840    0.5560]);
                grid on; axis equal;

                % Least Square Fitting Based on the Inliers
                circIdx(1) = find(leftInlierIdx,1,"last");
                circIdx(2) = find(rightInlierIdx,1,"first");
            end
            [isContinue,lineFitMaxDist] = checkextractfig(lineFitMaxDist,figShow);
        end

        scatterOri = scatterOri(:,circIdx(1):circIdx(2));
        hPlot1 = plot(scatterOri(1,:),scatterOri(2,:),'.','MarkerSize',8,'Color',[0.4660    0.6740    0.1880]);
        xlim = get(gca,'XLim');
        ylim = get(gca,'YLim');
        set(gca,"XLim",xlim,'YLim',ylim);
        isContinue = checkextractcir();
        while ~isContinue
            fprintf('Please select the useless area of original data: \n');
            disableDefaultInteractivity(gca);
            k = waitforbuttonpress; % 等待鼠标按下
            point1 = get(gca,'CurrentPoint'); % 鼠标按下了
            finalRect = rbbox; %
            point2 = get(gca,'CurrentPoint'); % 鼠标松开了
            point1 = point1(1,1:2); % 提取出两个点
            point2 = point2(1,1:2);
            p1 = min(point1,point2); % 计算位置
            p2 = max(point1,point2);
            offset = abs(point1 - point2); % offset(1)表示宽，offset(2)表示高
            boxDel(1,:) = [p1(1), p1(1) + offset(1), p1(1) + offset(1), p1(1), p1(1)];
            boxDel(2,:) = [p1(2), p1(2), p1(2) + offset(2), p1(2) + offset(2), p1(2)];
            hplot = plot(boxDel(1,:),boxDel(2,:),'r'); %在原图上显示截取的区域
            clearInd = (scatterOri(1,:) >= p1(1)) & (scatterOri(1,:) <= p2(1)) ...
                & (scatterOri(2,:) >= p1(2)) & scatterOri(2,:) <= p2(2);
            scatterOri(:,clearInd) = [];
            [isContinue] = checkextractcir(hplot);
            delete(hPlot1);
            hPlot1 = plot(scatterOri(1,:),scatterOri(2,:),'.','MarkerSize',8,'Color',[0.4660    0.6740    0.1880]);
        end
        fprintf('\nSuccessfully selected!\n\n');
        varargout{1} = lineFitMaxDist;
end

% tool tip arc fitting
[circ2D,RMSE] = arcFit2D(scatterOri, ...
    'arcFitMethod',options.arcFitMethod,'displayType',options.arcFitdisplayType);

% test whetger the ransac function defined above is true
function testRansac()
    cx0 = 0*1000; % unit: mu m
    cy0 = 0*1000; % unit: mu m
    r0 = 0.1*1000; % unit: mu m
    noise = 0.03;
    theta = (linspace(0,2*pi/3,200))';
    r = r0*(1 - noise + 2*noise*rand(length(theta),1));
    oriPts(1,:) = r.*cos(theta) + cx0;
    oriPts(2,:) = r.*sin(theta) + cy0;
    rmse0 = sqrt(...
                sum(...
                    ((oriPts(1,:) - cx0).^2 + (oriPts(2,:) - cy0).^2 - r0^2).^2) ...
            /length(theta));
    fitCirc = fitLineFcn(oriPts');
    figure('Name','Function Testification');
    plot(oriPts(1,:),oriPts(2,:),'.','Color',[0,0.45,0.74]); hold on;
    % plot the fitting center of the circle
    scatter(fitCirc.center(1),fitCirc.center(2),36,[0.6350,0.0780,0.1840],'filled');
    quiver(fitCirc.center(1),fitCirc.center(2), ...
        0.6*fitCirc.radius*fitCirc.startV(1),0.6*fitCirc.radius*fitCirc.endV(2), ...
        'filled','Color',[0.6350,0.0780,0.1840]);
    % plot the fitting circle
    scaThe = linspace(0,2*pi);
    scat(1,:) = fitCirc.radius*cos(scaThe) + fitCirc.center(1);
    scat(2,:) = fitCirc.radius*sin(scaThe) + fitCirc.center(2);
    plot(scat(1,:),scat(2,:),'k--','LineWidth',1);
    % plot the fitting arc
    R = rotz(atan2(fitCirc.startV(2),fitCirc.startV(1)) - 0.5*fitCirc.op);
    scaThe = linspace(0,fitCirc.openAng);
    scat(1,:) = fitCirc.radius*cos(scaThe);
    scat(2,:) = fitCirc.radius*sin(scaThe);
    circFit = R(1:2,1:2)*scat + fitCirc.center;
    plot(circFit(1,:),circFit(2,:),'Color',[0.8500,0.3250,0.0980],'LineWidth',3);
    hold off;
    grid on;
    axis equal
    evalPer = evalLineFcn(fitCirc,oriPts');
    evalSum = sum(evalLineFcn(fitCirc,toolOri));
end


%% rigid transform
n = size(scatterOri,2);
rotAng = -90 - atan2d(circ2D.arcVec(2),circ2D.arcVec(1));
rotMat = rotz(rotAng);
rotMat = rotMat(1:2,1:2);
scatterDst = rotMat*(scatterOri - ndgrid(circ2D.center,1:n));
end