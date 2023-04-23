function dist = pt2Line(pt,varargin)
% to calculate the distance form the point pt to the line line in a 2D
% space


switch nargin
    case 2
        line = varargin{1};
        switch length(line)
            case 3
                A = line(1);
                B = line(2);
                C = line(3);
            case 2
                A = line(1);
                B = -1;
                C = line(2);
        end
        dist = abs(A*pt(1,:) + B*pt(2,:) + C)/sqrt(A^2 + B^2);
    case 3
        pt1 = varargin{1};
        pt2 = varargin{2};
        len = size(pt,2);
        dist = vecnorm(cross(pt - pt1,ndgrid(pt2 - pt1,1:len),1),2,1)./norm(pt2 - pt1);
end

end