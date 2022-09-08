function pt = lineIntersectPlane(linePt1,linePt2,surfPt,surfNorm)
t = dot(surfPt - linePt1,surfNorm)/dot(linePt2 - linePt1,surfNorm);
pt = (1 - t)*linePt1 + t*linePt2;
end