function pt = lineIntersectPlane(linePt1,linePt2,surfPt,surfNorm)
% To calculate the intersection point between a line and a plane.
% 
% pt = lineIntersectPlane(linePt1,linePt2,surfPt,surfNorm)
% linePt1 & linePt2 (3,1) the two points on the line
% surfPt (3,1) the point on the plane
% surfNorm (3,1) the normal vector of the plane
t = dot(surfPt - linePt1,surfNorm)/dot(linePt2 - linePt1,surfNorm);
pt = (1 - t)*linePt1 + t*linePt2;
end