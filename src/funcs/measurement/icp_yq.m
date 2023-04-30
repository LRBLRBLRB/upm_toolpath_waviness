function [R,T] = icp_yq(P,Q)
% 读取源点云数据P
P0 = P;
%读取目标点云数据Q

difference = zeros(size(Q,1),3);   %创建记录坐标差值矩阵
mapPoint = zeros(size(P,1),3);    %创建与源点集对应的映射矩阵
distance = zeros(size(Q,1),1);    %创建记录所有距离的矩阵
n = size(P,1);

tic
j = 0; 
d = 100;
while d > 1
    j = j + 1;
%     fprintf("迭代次数：%d\n",j);
    %%
    %寻找点云数据P中点Pi的对应点Qi
    %使用欧式距离最近的点作为对应点
    parfor i=1:n
        difference=Q-P(i,:);       %坐标差值
        difference=difference.^2;
        distance=sum(difference,2);
        [~,minIndex]=min(distance);
        mapPoint(i,:)=Q(minIndex,:);    %保存对应点
    end
%     ShowData(P,Q);
%     hold on
%     for i = 1:n
%         plot3([P(i,1) mapPoint(i,1)], [P(i,2) mapPoint(i,2)],[P(i,3) mapPoint(i,3)], 'y--')
%     end
%     title('对应点');
    %%
    %计算旋转矩阵R和平移矩阵t的最优解，使用svd方法
    centerP=mean(P);    %P点集的质心点
    centerMap=mean(mapPoint);      %对应点集的质心点
    tempP=P-centerP;        %进行去中心化
    tempMapPoint=mapPoint-centerMap;
    
    H=tempP'*tempMapPoint;  %得到H矩阵
    [U,~,V]=svd(H);
    R=V*U';
    %     T=(centerP-centerMap)';
    T=-R*centerP'+centerMap';   %利用质心点求解T参数
    %%
    %使用R和T来得到新的点集
    P=(R*P'+T)';       %使用转换参数得到新的点集P
    d=sum(sum((P-mapPoint).^2,2))/n;	%计算新的点集P到对应点的平均距离
    if j>400
        break
    end
end
toc
ShowData(P,Q);

%%求解变换矩阵
mapPoint=P;
P=P0;
centerP=mean(P);    %P点集的质心点
centerMap=mean(mapPoint);      %对应点集的质心点
tempP=P-centerP;        %进行去中心化
tempMapPoint=mapPoint-centerMap;
H=tempP'*tempMapPoint;  %得到H矩阵
[U,~,V]=svd(H);
R=V*U';
%     T=(centerP-centerMap)';
T=-R*centerP'+centerMap';   %利用质心点求解T参数


%%理论变换矩阵
% mm=[1 0 0 70;
%     0 cos(pi/6) -sin(pi/6) 70;
%     0 sin(pi/6) cos(pi/6) 20;
%     0 0 0 1];
% mm=inv(mm)