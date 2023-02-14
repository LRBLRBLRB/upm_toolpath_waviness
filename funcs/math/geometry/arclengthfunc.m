function arcFunc = arclengthfunc(curvFunc)
%ARCLENGTHFUNC Generate the function of arc length 


syms t;
curvFunc_t = diff(curvFunc,t);
int(curvFunc_t,t)
arcFunc = matlabFunction(arcFunc,'Vars',t);

end

