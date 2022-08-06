function Nip = bBasisFun(i,u,p,U)
% 计算B样条基函数在参数u处的值
% 输入：基函数下标i（1~m-p），基函数次数p，参数u，节点矢量U
% 输出：B样条基函数在u处的值

m = length(U); % which means (m+1)
if (i==1 && u==U(1)) || (i==m-p && u==U(m))
    Nip = 1;
    return;
end
if (u<U(i)) || (u>=U(i+p+1)) % none of the Ni0 in the interval remains 0
    Nip = 0;
    return;
end

% initialize the 0-order basis functions
N = zeros(p+1,1);
Nmin = find(u>=U);
Nmax = find(u<U);
N(Nmin:Nmax) = 1;

% calculate the N table
for k = 1:p
    if N(1)==0
        saved = 0;
    else
        saved = (u-U(i)*N(1))/(U(i+k)-U(i));
    end
    for j = 0:p-k+1
        Uleft = U(i+j+1);
        Uright = U(i+j+k+1);
        if N(j+2)==0
            N(j+1) = saved;
            saved = 0;
        else
            temp = N(j+2)/(Uright-Uleft);
            N(j+1) = saved + (Uright-u)*temp;
            saved = (u-Uleft)*temp;
        end
    end
end
Nip = N(1);
end
