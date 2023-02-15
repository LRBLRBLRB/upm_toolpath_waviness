function fDer = numericdiff(f,N,dim)
%NUMERICDIFF calculate the numerical differentiates of sequence
%   二阶导以上还有问题的！

arguments
    f
    N {mustBeGreaterThanOrEqual(N,1),mustBeInteger} = 1
    dim {mustBeMember(dim,[1,2])} = 1
end

fDer = f;

for n = 1:N
    xDiff = diff(fDer,1,dim); % the difference of xDer
    if dim == 1
        forw = [xDiff(1,:);xDiff];
        backw = [xDiff;xDiff(end,:)];
    else 
        forw = [xDiff(:,1),xDiff];
        backw = [xDiff,xDiff(:,end)];
    end
    fDer = (forw + backw)/2;
end

end