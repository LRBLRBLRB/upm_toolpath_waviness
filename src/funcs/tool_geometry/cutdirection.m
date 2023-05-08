function toolDirect = cutdirection(surfPt,surfCenter,angular,options)
% to solve the turning direction of the tool edge, provide the tool
% direction of the tool edge at the exact points that are given.
% usage: 
%   toolDirect = cutDirection(surfPt,surfCenter)
% Inputs:
%   surfPt (3,:) the surface point array
%   surfCenter (3,1) the rotation center or the point on the rotation axis
% Outputs:
%   toolDirect (3,:) the cutting direction of each point in the array
%
%   toolDirect = cutDirection(surfPt)

arguments
    surfPt (3,:) double
    surfCenter (3,1) double {mustBeVector} = [inf;inf;inf]
    angular (1,:) double {mustBeVector} = [0,0]
    options.method {mustBeMember(options.method, ...
        {'concentric','adjacent','vertical'})} = 'adjacent'
end

switch options.method
    % this case is activated when the spindle is perpendicular to the plane
    case 'vertical'
        % if diffAngular < 0, it means the tool point is counterclockwise,
        % and that the direction vector should be [-sin(theta),cos(theta)],
        % otherwise it should be clockwise.
        diffAngular = sign(angular(2) - angular(1));
        if ~diffAngular
            error(['Invalid input: concentric angle should be given and ', ...
                'should not be identical when using the method ''vertical''']);
        end
        toolDirect = [-sin(angular)*diffAngular;
            cos(angular)*diffAngular;
            zeros(1,length(angular))];
        toolDirect = toolDirect./vecnorm(toolDirect,2,1);
    case 'concentric'
        % 思路一：对于回转体的同心圆算法（concentric），可以直接求回转轴指向曲面点的向量，然后投影到回转轴垂直的表面，再取x,y后右转90度。
        if any(surfCenter == inf)
            error(['Invalid input: surfcae center should be given and ', ...
                'should not be inf when using the method ''concentric''']);
        end
        Rot = rotz(90);
        toolDirect = Rot*(surfPt - surfCenter);
        toolDirect(3,:) = 0;
        toolDirect = toolDirect ./ vecnorm(toolDirect,2,1);
    case 'adjacent'
        % 思路二：直接拿surfMesh中的下一个点和上一个点的坐标，然后算方向：
        surfPt = [surfPt(:,1),surfPt,surfPt(:,1)];
        toolDirect = 0.5*(surfPt(:,3:end) - surfPt(:,1:end-2));
        toolDirect = toolDirect ./ vecnorm(toolDirect,2,1);
    otherwise
        error("Invalid input");
end

end