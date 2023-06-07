function surfMesh = zygo_datx(surfFilePath,fileUnit)
%ZYGO_DATX 此处显示有关此函数的摘要
%   此处显示详细说明

surfDataStruct = load(surfFilePath);
surfMesh = zeros([size(surfDataStruct.mh),3]);
[surfMesh(:,:,2),surfMesh(:,:,1)] = meshgrid(surfDataStruct.vx,surfDataStruct.vy);
surfMesh(:,:,3) = surfDataStruct.mh;

%% unit convertion
unitList = {'m','mm','um','nm'};
xyUnit0 = find(strcmp(unitList,surfDataStruct.XYunit),1); % default unit is milimeter
zUnit0 = find(strcmp(unitList,surfDataStruct.Zunit),1); % default unit is meter
aimUnit = find(strcmp(unitList,fileUnit),1);
surfMesh(:,:,1:2) = 1000^(aimUnit - xyUnit0)*surfMesh(:,:,1:2);
surfMesh(:,:,3) = 1000^(aimUnit - zUnit0)*surfMesh(:,:,3);

%% rigid transform
% surfMesh(:,:,1) = transpose(surfMesh(:,:,1));
% surfMesh(:,:,2) = transpose(surfMesh(:,:,2));
% surfMesh(:,:,3) = transpose(surfMesh(:,:,3));
sizeSurf = size(surfMesh);
surfPt = reshape(surfMesh,[],3);
surfPt = surfPt*transpose(rotz(-90));
surfMesh = reshape(surfPt,sizeSurf(1),sizeSurf(2),3);
% surfMesh(:,:,3) = surfMesh(:,end:-1:1,3);

end