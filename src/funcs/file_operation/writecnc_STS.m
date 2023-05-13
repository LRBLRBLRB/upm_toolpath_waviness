function writecnc_STS(ncFile,workpiece,tool,axisC,axisX,axisZ,loop)
%WRITECNC_2D 此处显示有关此函数的摘要
%   此处显示详细说明

isLoop = isfield(loop,{'num','offset','step'});
if ~all(isLoop)
    error('Invalid struct LOOP! \n');
end

ncFid = fopen(ncFile,'w');

%% head
fprintf(ncFid,'#100 = %d\t\t( TOTAL PASSES OF Z)\n',loop.num);
fprintf(ncFid,'#103 = %f\t\t( DEPTH OF CUT Z )\n',loop.step);
fprintf(ncFid,'#555 = 0.001\t\t( FEED RATE )\n\n');

fprintf(ncFid,'\nG52 G63 G71 G103 G40 G18 G90\n');
fprintf(ncFid,'\n( WORKING COORDINATES )\n%s\n%s\n\n',workpiece,tool);

fprintf(ncFid,'G18 G40 G94\n');
fprintf(ncFid,'M78\n');
fprintf(ncFid,'G94 Z20 F500\n\n');

fprintf(ncFid,'M26.1\n\n');

fprintf(ncFid,'#5 = 0\t\t( COUNT VARIABLE )\n');
fprintf(ncFid,'#8 = %f\t\t( CUT OFFSET Z )\n',loop.offset);
fprintf(ncFid,'WHILE[#5 LT #100] DO 2\n\n');

fprintf(ncFid,'#8 = #8 - #103\t\t( Z AXIS CUTTING )\n');
fprintf(ncFid,'G52 X0 Y0 Z[#8]\n');

fprintf(ncFid,'( LINK BLOCK )\n');
fprintf(ncFid,'G94 X%f\n',axisX(1));
fprintf(ncFid,'G93 F0.001 C0\n\n');

fprintf(ncFid,'G94 Z5 F500\n');
fprintf(ncFid,'Z%f F200\n',ceil(axisZ(1)));
fprintf(ncFid,'Z%f F20\n',axisZ(1) + 0.1);
fprintf(ncFid,'Z%f F1\n\n',axisZ(1));

fprintf(ncFid,'( CUTTING BLOCK )\nG93 F[#555]\n');
for ii = 1:length(axisC)
    fprintf(ncFid,'C%f X%f Z%f\n',axisC(ii),axisX(ii),axisZ(ii));
%     fprintf(ncFid,'%f %f %f\n',axisC(ii),axisX(ii),axisZ(ii));
%     fprintf(ncFid,'GOTO/%f,%f,%f,%f,%f,%f\n', ...
%         spiralPath1(1,ii),spiralPath1(2,ii),spiralPath1(3,ii), ...
%         spiralNorm(1,ii),spiralNorm(2,ii),spiralNorm(3,ii));
end

fprintf(ncFid,'(LINKING BLOCK)\n');
fprintf(ncFid,'G94 Z5 F200\n\n');
fprintf(ncFid,'#5 = #5 + 1\n');
fprintf(ncFid,'END 2\n\n');
fprintf(ncFid,'G94 Z20 F800\n');
fprintf(ncFid,'M29\n');
fprintf(ncFid,'M30\n');

fclose(ncFid);

end