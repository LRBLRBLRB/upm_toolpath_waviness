function pt = curveIntersectLineOnPlane(cirPt1,cirPt2,surfPt)
% To calculate the intersection point between a circular curve and a plane.
% 
% pt = lineIntersectPlane(linePt1,linePt2,surfPt,surfNorm)
% linePt1 & linePt2 (3,1) the two points on the line
% surfPt (3,1) the point on the plane

r = 0.5*(norm(cirPt1(1:2)) + norm(cirPt2(1:2)));
r1 = norm(surfPt(1:2));
pt = (r/r1)*surfPt(1:2);

t = (atan2(pt(2),pt(1)) - atan2(cirPt1(2),cirPt1(1))) ...
    /(atan2(cirPt2(2),cirPt2(1)) - atan2(cirPt1(2),cirPt1(1)));
pt(3) = cirPt1(3) + t*(cirPt2(3) - cirPt1(3));

end