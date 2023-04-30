function [R,T] = icp_yq(P,Q)
% ��ȡԴ��������P
P0 = P;
%��ȡĿ���������Q

difference = zeros(size(Q,1),3);   %������¼�����ֵ����
mapPoint = zeros(size(P,1),3);    %������Դ�㼯��Ӧ��ӳ�����
distance = zeros(size(Q,1),1);    %������¼���о���ľ���
n = size(P,1);

tic
j = 0; 
d = 100;
while d > 1
    j = j + 1;
%     fprintf("����������%d\n",j);
    %%
    %Ѱ�ҵ�������P�е�Pi�Ķ�Ӧ��Qi
    %ʹ��ŷʽ��������ĵ���Ϊ��Ӧ��
    parfor i=1:n
        difference=Q-P(i,:);       %�����ֵ
        difference=difference.^2;
        distance=sum(difference,2);
        [~,minIndex]=min(distance);
        mapPoint(i,:)=Q(minIndex,:);    %�����Ӧ��
    end
%     ShowData(P,Q);
%     hold on
%     for i = 1:n
%         plot3([P(i,1) mapPoint(i,1)], [P(i,2) mapPoint(i,2)],[P(i,3) mapPoint(i,3)], 'y--')
%     end
%     title('��Ӧ��');
    %%
    %������ת����R��ƽ�ƾ���t�����Ž⣬ʹ��svd����
    centerP=mean(P);    %P�㼯�����ĵ�
    centerMap=mean(mapPoint);      %��Ӧ�㼯�����ĵ�
    tempP=P-centerP;        %����ȥ���Ļ�
    tempMapPoint=mapPoint-centerMap;
    
    H=tempP'*tempMapPoint;  %�õ�H����
    [U,~,V]=svd(H);
    R=V*U';
    %     T=(centerP-centerMap)';
    T=-R*centerP'+centerMap';   %�������ĵ����T����
    %%
    %ʹ��R��T���õ��µĵ㼯
    P=(R*P'+T)';       %ʹ��ת�������õ��µĵ㼯P
    d=sum(sum((P-mapPoint).^2,2))/n;	%�����µĵ㼯P����Ӧ���ƽ������
    if j>400
        break
    end
end
toc
ShowData(P,Q);

%%���任����
mapPoint=P;
P=P0;
centerP=mean(P);    %P�㼯�����ĵ�
centerMap=mean(mapPoint);      %��Ӧ�㼯�����ĵ�
tempP=P-centerP;        %����ȥ���Ļ�
tempMapPoint=mapPoint-centerMap;
H=tempP'*tempMapPoint;  %�õ�H����
[U,~,V]=svd(H);
R=V*U';
%     T=(centerP-centerMap)';
T=-R*centerP'+centerMap';   %�������ĵ����T����


%%���۱任����
% mm=[1 0 0 70;
%     0 cos(pi/6) -sin(pi/6) 70;
%     0 sin(pi/6) cos(pi/6) 20;
%     0 0 0 1];
% mm=inv(mm)