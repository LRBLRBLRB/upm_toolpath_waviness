function data = read_STS(filePath)
%IMPORTNC 此处显示有关此函数的摘要
%   此处显示详细说明

%% initialize the data
tic
sprintf("Finding the size of data file: \r\n%s",file)
fileID = fopen(filePath);

row = 0;
while ~feof(fileID)
    % read 10000 str, calculate the number of \n (char(10)=\n)
    % '*char' indicates 1 str per read, *indicates output str
    row = row + sum(fread(fileID,10000,'*char') == newline);
end
fclose(fileID);
fprintf('File row number is %d.\r\n',row)
toc
MaxCounter = row;

%% read data from file



end

