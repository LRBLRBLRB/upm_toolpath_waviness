[filename pathname]=uigetfile({'*.txt','txt-file(*.txt)';'*.*','All the files(*.*)'},'Choose a file');

if isequal(filename,0)||isequal(pathname,0);
    h=msgbox ('Please choose a file!','Warning','warn');
    return;
else
% [FileName,PathName] = uigetfile('*.KGF','Select the Data file'); 
    path2=fullfile(pathname,filename)
    fidin=fopen(path2); %打开文件
  
    if  ~feof(fidin)  
     tline=fgetl(fidin);  
     str=[tline,10]
    end
    while ~feof(fidin)  
     tline=fgetl(fidin);  
     str1=[tline,10]
     str=[str,str1]
    end
    set(handles.editCode,'string',str1)

end

fclose(fidin);