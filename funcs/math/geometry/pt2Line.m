function dist = pt2Line(pt,line)
% to calculate the distance form the point pt to the line line in a 2D
% space

arguments
    pt (2,:)
    line {mustBeFinite}
end

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

end