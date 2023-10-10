function data = unitconversion(inputMode,data, ...
    inputLengthUnit,outputLengthUnit,inputAngleUnit,outputAngleUnit)
%UNITCONVERSION convert the unit of the input data

arguments
    inputMode cell
    data
    inputLengthUnit {mustBeMember(inputLengthUnit,{'m','mm','\mum','nm'})}
    outputLengthUnit {mustBeMember(outputLengthUnit,{'m','mm','\mum','nm'})}
    inputAngleUnit {mustBeMember(inputAngleUnit,{'deg','rad'})} = 'rad'
    outputAngleUnit {mustBeMember(outputAngleUnit,{'deg','rad'})} = 'rad'
end

lengthUnitList = {'m','mm','\mum','nm'};
inputLengthUnitNo = find(strcmp(lengthUnitList,inputLengthUnit),1);
outputLengthUnitNo = find(strcmp(lengthUnitList,outputLengthUnit),1);

angleUnitList = {'deg','rad'};
inputAngleUnitNo = find(strcmp(angleUnitList,inputAngleUnit),1);
outputAngleUnitNo = find(strcmp(angleUnitList,outputAngleUnit),1);

num = length(inputMode);
if iscell(data)
    for ii = 1:num
        switch inputMode{ii}
            case 'length'
                data{ii} = 1000^(outputLengthUnitNo - inputLengthUnitNo)*data{ii};
            case 'angle'
                data{ii} = (pi/180)^(outputAngleUnitNo - inputAngleUnitNo)*data{ii};
        end
    end
else
    for ii = 1:num
        switch inputMode{ii}
            case 'length'
                data(ii,:) = 1000^(outputLengthUnitNo - inputLengthUnitNo)*data(ii,:);
            case 'angle'
                data(ii,:) = (pi/180)^(outputAngleUnitNo - inputAngleUnitNo)*data(ii,:);
        end
    end
end

end