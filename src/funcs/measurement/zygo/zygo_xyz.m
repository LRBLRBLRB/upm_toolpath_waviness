function surfMesh = zygo_xyz(filePath,fileUnit)
%ZYGO_XYZ to load the xyz file from Zygo Nexview™ NX2 Coherent scanning
% interferometric profiler, and to extract the surf mesh data

if nargin == 1
    fileUnit = 'um';
end

%% get the settings of the data
fid = fopen(filePath,'r');
numHeader = 0;
while ~feof(fid)
    tmpLine = fgetl(fid);
    numHeader = numHeader + 1;
    % if the line begins with %d%d or -%d, then break
    if strcmp(tmpLine,'#')
        break;
    end
end

%% surface data extraction
surfDataStruct = importdata(filePath,' ',numHeader);

tmp = sscanf(surfDataStruct.textdata{4},'%d %d %d %d');
xNum = tmp(3);
yNum = tmp(4);
tmp = sscanf(surfDataStruct.textdata{8},'%f %f %f %f %f %f %f %f');
dxy = tmp(7);

surfData(:,1:2) = dxy*surfDataStruct.data(:,1:2);

surfMesh(:,:,1) = reshape(surfData(:,1),[yNum,xNum])';
surfMesh(:,:,2) = reshape(surfData(:,2),[yNum,xNum])';
surfMesh(:,:,3) = reshape(surfDataStruct.data(:,3),[xNum,yNum])';

%% unit convertion
unitList = {'m','mm','um','nm'};
xyUnit0 = 1; % default unit is milimeter
zUnit0 = 3; % default unit is meter
aimUnit = find(strcmp(unitList,fileUnit),1);
surfMesh(:,:,1:2) = 1000^(aimUnit - xyUnit0)*surfMesh(:,:,1:2);
surfMesh(:,:,3) = 1000^(aimUnit - zUnit0)*surfMesh(:,:,3);

end