function structSurf = loadXYZ(filePath,xNum,yNum,dxy,xyUnit,zUnit)
%% Latest revision date: 2021.11.15
%% Form
%  struct = loadXYZ(filePath, xNum, yNum, dxy, xyUnit, zUnit)
%% Description
%  Get Mat data from a Zygo interferometer
%% Inputs
%  filePath      string             
%  xNum,yNum     1                  get from original file of XYZ
%  dxy           1                  get from original file of XYZ
%  xyUnit,zUnit  'um'/'mm'/'mm'     the unit in XYZ file is 'm' for XY, 'um' for Z
%% Outputs
%  structSurf    X,Y,Z  (yNum,xNum)   meshgrid XY and height data Z 

%%
    if (nargin<5)
        xyUnit = 'mm';
        zUnit = 'um';
    end
    
    matrixSurf = load(filePath); 
    X = 0:dxy:dxy*(xNum-1);
    for i=1:yNum
        startInd = (i-1)*xNum + 1 ;
        endInd = i*xNum ;
        matrixSurf(startInd:endInd,1) = X';
        matrixSurf(startInd:endInd,2) = ones(xNum,1)*(i-1)*dxy;
    end
    X = reshape(matrixSurf(:,1),xNum,yNum)';
    Y = reshape(matrixSurf(:,2),xNum,yNum)';
    Z = reshape(matrixSurf(:,3),xNum,yNum)';
    Z = Z(end:-1:1,:);
    
    switch xyUnit
        case 'um'
            structSurf.X = X*1e6;
            structSurf.Y = Y*1e6;
            structSurf.dXY = dxy*1e6;
        case 'mm'
            structSurf.X = X*1000;
            structSurf.Y = Y*1000;
            structSurf.dXY = dxy*1000;
        case 'm'
            structSurf.X = X;
            structSurf.Y = Y;
            structSurf.dXY = dxy;
    end

    switch zUnit
        case 'um'
            structSurf.Z = Z;
        case 'm'
            structSurf.Z = Z*1e-6;
        case 'mm'
            structSurf.Z = Z*1e-3;
    end
    
    structSurf.xy_unit = xyUnit;
    structSurf.Z_unit = zUnit;
end