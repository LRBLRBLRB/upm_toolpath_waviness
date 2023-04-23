function varargout = arithmeticsequences(a)
%ARITHMETICSEQUENSES get the arithmetic sequenses from the array a

len = length(a);
num = 0;
aArith = zeros(2,num);
%% double pointers method
% i = 1;
% j = 2;
% while i < len - 1 && j <= len
%     der = a(i + 1) - a(i);
%     count = 0;
%     for j = i + 1:len
%         if a(j) - a(j - 1) == der
%             count = count + 1;
%         else
%             num = num + 1;
%             aArith(:,num) = [];
%             i = j - 1;
%             break;
%         end
%     end
% end

%% vectorized method
if len
    if size(a,1) ~= 1
        a = a';
    end
    isArith = find(diff(a) ~= 1);
    % isArith(diff(isArith) == 1) = []
    num = length(isArith) + 1;
    isArith = [0,isArith,length(a)];
    aArith = [isArith(1:end - 1) + 1;isArith(2:end)];
    % aArith(:,aArith(1,:) == aArith(2,:)) = [];
end

varargout{1} = num;
varargout{2} = aArith;
end