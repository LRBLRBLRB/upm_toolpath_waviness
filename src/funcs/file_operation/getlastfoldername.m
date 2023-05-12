function foldername = getlastfoldername(path)
%GETLASTFOLDERNAME to get the name of the last folder
    [~, foldername, extName] = fileparts(path); % 在路径字符串的结尾添加了斜杠符号/
    foldername = [foldername,extName];
    if isempty(foldername) % 若% 在路径字符串的结尾没有添加斜杠符号/
        folders = split(path, filesep);
        foldername = folders{end};
    end
end