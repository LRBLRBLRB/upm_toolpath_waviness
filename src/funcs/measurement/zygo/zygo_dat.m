function data = zygo_dat
%ZYGO_DAT 此处显示有关此函数的摘要
%   此处显示详细说明

fileName = '030130-DEPTH5+FEED0.001-20X+ZOOM0.5+AVE3-RAW.datx';

fileID = fopen(fileName,'rb');
[data,count]=fread(fileID,Inf,"char");

fclose(fileID);

end