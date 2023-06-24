%% initialize data size
clear
clc
close

path = "..\workspace\20230510\01000-sts-r1000c0.09100-D906+0417-res1.0-arc0.02+deg1";
filename = "01000-sts-r1000c0.09100-D906+0417-res1.0-arc0spiralPath20230510T224816.nc";
file = path + "\" + filename;

tic
sprintf("Finding the size of data file: \r\n%s",file)
fid=fopen(file);

row=0;
while ~feof(fid)
    % read 10000 str, calculate the number of \n (char(10)=\n)
    % '*char' indicates 1 str per read, *indicates output str
    row=row+sum(fread(fid,10000,'*char')==newline);
end
% fclose(fid);
fprintf('File row number is %d.\r\n',row)
toc
MaxCounter = row;



%% read data from file
tic
fprintf('Importing data...\r\n')
% fid=fopen(file);

Wait_Title = waitbar(0,'Data importing...');

A=zeros(MaxCounter,1);
A=string(A);

for i = 1 : MaxCounter
    
    if mod(i,fix(MaxCounter/100))==0
        Display_Data = num2str(roundn(i/MaxCounter*100,-2));
        % Calculate percentage
        Display_Str = ['Import Progress: ',Display_Data,'%'] ;
        % Show Calculate State
        waitbar(i/MaxCounter,Wait_Title,Display_Str);
        % Progress bar dynamic display
    end
    
    
    str = fgetl(fid);
    A(i) = string(str);
    
end
close(Wait_Title);   % Close Progress bar window
fclose(fid);
fprintf('Data imported.\r\n')
toc

%% Determine startup position
for i = 1 : 100
    
    if A(i)=="( CUTTING BLOCK )" || A(i)=="( LINK BLOCK )"
        startup = i;
        fprintf('Start up from line %d.\r\n',i)
        break
    end
    
    %if handled
    if i==100
        startup=1;
    end
    
end


%% Parallel process data
wrong_line_number=0;
data=zeros(MaxCounter,3);
% parpool(6)
tic
fprintf('Parallel processing data...\r\n')

parfor i = startup : MaxCounter
    
    if (A(i)=="( CUTTING BLOCK )")
        data(i,:) = [1111 0 0]; % data(i,:)=[1111 0 0]
        continue;
    end
    
    if (A(i)=="( LINK BLOCK )")
        data(i,:) = [2222 0 0]; % data(i,:)=[2222 0 0]
        continue;
    end
    
    %if handled cutting/linking flag
    if (A(i)=="1111")
        data(i,:) = [1111 0 0]; % data(i,:)=[2222 0 0]
        continue;
    end
    
    if (A(i)=="2222")
        data(i,:) = [2222 0 0]; % data(i,:)=[2222 0 0]
        continue;
    end
    
    if (startsWith(A(i),'C')+startsWith(A(i),'X')+startsWith(A(i),'Z')) == 0
        fprintf('Line %d format error.\r\n',i)
        wrong_line_number=wrong_line_number+1;
        continue
    end
    
    splitStr = split(A(i)," ");
    TempData = 3333*ones(3,1);
    
    for j = 1:length(splitStr)
        
        if startsWith(splitStr(j),"C") == 1
            TempData(1)=sscanf(splitStr(j),'C%f');
        end
        
        if startsWith(splitStr(j),"X") == 1
            TempData(2)=sscanf(splitStr(j),'X%f');
        end
        
        if startsWith(splitStr(j),"Z") == 1
            TempData(3)=sscanf(splitStr(j),'Z%f');
        end
        
        
    end
    data(i,:) = TempData;
    continue
    
end

fprintf('Data processed.\r\n')
toc
delete(gcp('nocreate'))  %delete pool generated

%% eliminate [3333 3333 3333 3333 3333]
tic
fprintf('Eliminating [3333 3333 3333]...\r\n')

flag=(data(:,1)==3333&data(:,2)==3333&data(:,3)==3333);
data2=data(~flag,:);

fprintf('data eliminated and saved as data2.\r\n')
toc

%% complement data2
tic
fprintf('Complementing data2\r\n')
data3 = data2;

for i = 2 : length(data3)
    
    if data3(i-1,1) == 1111 || data3(i-1,1) == 2222
        if i==2
            continue
        end
    end
    
    if data3(i,1) == 3333
        if data3(i-1,1) == 1111 || data3(i-1,1) == 2222  %eliminate the influence of cutting/linkage code
            data3(i,1) = data3(i-2,1);
        else
            data3(i,1) = data3(i-1,1);
        end
    end
    
    if data3(i,2) == 3333
        if data3(i-1,1) == 1111 || data3(i-1,1) == 2222
            data3(i,2) = data3(i-2,2);
        else
            data3(i,2) = data3(i-1,2);
        end
    end
    
    if data3(i,3) == 3333
        if data3(i-1,1) == 1111 || data3(i-1,1) == 2222
            data3(i,3) = data3(i-2,3);
        else
            data3(i,3) = data3(i-1,3);
        end
    end
      
end
fprintf('data2 complemented and saved as data3\r\n')
toc

%% Output the completed array
tic
data3=data3';
fid2=fopen("211122_SelfRearch_completed.nc","wt");


fprintf(fid2,"C%.6f X%.6f Z%.6f\n",data3);

fclose(fid2);
fclose all;
toc