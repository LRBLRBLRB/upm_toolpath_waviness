function toolDirect = toolDirection(surfPt,surfCenter)
% usage: toolDirect = toolDirection(surfPt,surfCenter)
% to solve the turning direction of the tool edge, provide the tool
% direction of the tool edge at the exact points that are given.
switch nargin
    case 2
        % 思路一：对于回转体的同心圆算法（concentric），可以直接求回转轴指向曲面点的向量，然后投影到回转轴垂直的表面，再取x,y后右转90度。
        I = rotz(pi/2);
        toolDirect = (surfPt - surfCenter)*I';
    case 1
        % 思路二：直接拿surfMesh中的下一个点和上一个点的坐标，然后算方向：
        surfPt = [surfPt(1,:);surfPt;surfPt(end,:)];
        toolDirect = 0.5*(surfPt(3:end,:) + surfPt(1:end-2,:) - 2*surfPt(2:end-1,:));
    otherwise
        error("Invalid input");
end

end