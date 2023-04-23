function vec2 = vecperp(vec1,n)
% calculate the other component of vector vec1, while one of it is the n
% vector

vec2 = n - dot(vec1,n)/norm(vec1)*vec1;

end