function [planeNorm,varargout] = fitPlane(pts,options)
%FITPLANE 此处显示有关此函数的摘要
%   此处显示详细说明


arguments
    pts (:,3)
    options.fitMethod {mustBeMember(options.fitMethod,{ ...
        'SVD','regression',''})} = 'SVD'
    options.normalization logical = true
end

%% Decentralization
pts = pts(~isnan(pts(:,3)),:);
num = size(pts,1);
planePt = mean(pts,1);

switch options.fitMethod
    case 'SVD'
        % SVD to fit the plane
        pts1 = pts - meshgrid(planePt,1:num);
        % svd: singular values are nonnegative and returned in decreasing order.
        [~,~,V] = svd(pts1,0,'vector'); 
        planeNorm = V(:,3); % normVec remains the norm vector of th e fitting plane
    case 'regression'
        regress()
end

if options.normalization
    dtmp = mean(pts*planeNorm);
    planeNorm = planeNorm*sign(dtmp);
    planeNorm = planeNorm./norm(planeNorm);
end

switch nargout
    case {0,1}
        planeNorm(4) = planePt*planeNorm;
    case 2
        varargout{1} = planePt;
end

end