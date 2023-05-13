function writecnc_2D(ncFile,workpiece,tool,axisX,axisZ,loop)
%WRITECNC_2D To write the cnc file of 2D-turning process

isLoop = isfield(loop,{'num','offset','step'});
if ~isLoop
    error('Invalid struct LOOP! \n');
end

ncFid = fopen(ncFile,'w');

%% head
fprintf(ncFid,'#105=10		( CUTTING SPEED )\n\n');

fprintf(ncFid,'#201=1500	( SPINDLE SPEED )\n\n');

fprintf(ncFid,'#501=%d		( ROUGH TIMES )\n',loop.num);
fprintf(ncFid,'#502=0		( CURRENT TIMES )\n');
fprintf(ncFid,'#503=%f	    ( G52 OFFSET )\n',loop.offset);
fprintf(ncFid,'#504=%f	    ( CUTTING DEPTH )\n\n',loop.step);

fprintf(ncFid,'G92.1\n');
fprintf(ncFid,'G71 G90 G01 G40 G94 G52 G144\n\n');

fprintf(ncFid,'%s		    ( WORKPIECE COORDINATE )\n',workpiece);
fprintf(ncFid,'%s		    ( TOOL COORDINATE )\n\n',tool);

fprintf(ncFid,'G94 Z20 F200 \n');
fprintf(ncFid,'M26.1\n');
fprintf(ncFid,'M03 S[#201]\n\n');

fprintf(ncFid,'WHILE[#502LT#501]DO 1\n\n');

fprintf(ncFid,'G52 Z[#503]\n');
fprintf(ncFid,'Z%d F200\n',ceil(axisZ(1)) + 1);
fprintf(ncFid,'X%f F200\n',axisX(1));
fprintf(ncFid,'Z%f F100\n',ceil(axisZ(1)) + 0.5);
fprintf(ncFid,'Z%f F50\n',axisZ(1) + 0.1);
fprintf(ncFid,'Z%f F[#105]\n\n',axisZ(1));

fprintf(ncFid,'( CUTTING BLOCK )\n');

%% body
for ii = 1:length(axisX)
    fprintf(ncFid,'X%f Z%f\n',axisX(ii),axisZ(ii));
end

%% tail
fprintf(ncFid,'( LINKING BLOCK )\n\n');

fprintf(ncFid,'G94 Z0.5 F100\n');
fprintf(ncFid,'Z2 F200\n\n');

fprintf(ncFid,'#503=#503+#504\n');
fprintf(ncFid,'#502=#502+1\n\n');

fprintf(ncFid,'END 1\n\n');

fprintf(ncFid,'G94 Z20 F200\n\n');

fprintf(ncFid,'M05\n');
fprintf(ncFid,'M29\n');
fprintf(ncFid,'M30\n');

fclose(ncFid);

end