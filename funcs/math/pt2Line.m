function dist = pt2Line(pt,line)

A = line(1);
B = line(2);
C = line(3);

dist = abs(A*pt(1) + B*pt(2) + C)/sqrt(A^2 + B^2);

end