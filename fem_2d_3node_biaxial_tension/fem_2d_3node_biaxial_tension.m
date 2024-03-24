clc;clear;
%% 前处理：生成单元、节点、材料参数
% 矩形板长1.2m，宽0.4m，厚度0.05m
L=1.2;W=0.4;thick=0.05;
%杨氏模量1.15e11，泊松比0.25
Ys=1.15e11;poss=0.25;
%单元划分
nx=24;      %横向24段
ny=8;       %纵向8段
%存储各个节点坐标,划分平面矩形单元，每个单元有3个节点
coordinate_node=zeros((nx+1)*(ny+1),2);
element=zeros(nx*ny*2,3);
%生成各个节点的坐标，在coordinate_node为（x，y）
%节点每横排（nx+1）个，每一列(ny+1)个
for i=1:ny+1
    coordinate_node(1+(nx+1)*(i-1):(nx+1)*i,1)=L/nx*(0:nx);
    coordinate_node(1+(nx+1)*(i-1):(nx+1)*i,2)=W/ny*(i-1)*ones(1,nx+1);
end
%给每个单元分配节点，节点以逆时针为顺序
for i=1:nx
    for j=1:ny 
        element(2*ny*(i-1)+j,1)=i+(nx+1)*(j-1);
        element(2*ny*(i-1)+j,2)=i+1+(nx+1)*(j-1);
        element(2*ny*(i-1)+j,3)=i+(nx+1)*j;
    end
    for j=1:ny
        element(2*ny*(i-1)+j+ny,1)=i+(nx+1)*j;
        element(2*ny*(i-1)+j+ny,2)=i+1+(nx+1)*(j-1);
        element(2*ny*(i-1)+j+ny,3)=i+1+(nx+1)*j;
    end
end

%% 生成单元刚度矩阵，组装总刚
%统计节点和单元数量
num_node=size(coordinate_node,1);
num_element=size(element,1);
%总刚矩阵记为KK
KK=zeros(num_node*2,num_node*2);
%生成单元刚度矩阵并组装到总刚中
for i=1:num_element
    xt=[coordinate_node(element(i,1),1), coordinate_node(element(i,2),1), coordinate_node(element(i,3),1)];
    yt=[coordinate_node(element(i,1),2), coordinate_node(element(i,2),2), coordinate_node(element(i,3),2)];
    %生成单元刚度矩阵
    Ke=stiffness_element_2d_3node(Ys,poss,thick,xt,yt);
    %组装单元刚度矩阵
    KK=stiffness_element_2d_3node_assemble(KK,Ke,element(i,:));
end

%% 面力载荷移置
%设置左右面受大小为P的均布拉伸面力
%R为载荷列阵
%P=20MPa
P=[20e6 -20e6];
R=zeros(num_node*2,1);
%RP为每个节点所受的力
RP=P*thick*(W/ny)/2;
%加载节点载荷
for i=1:ny
    R(element(i,1)*2-1)=R(element(i,1)*2-1)+RP(2);
    R(element(i,3)*2-1)=R(element(i,3)*2-1)+RP(2);
end
for i=(num_element-7):num_element
    R(element(i,3)*2-1)=R(element(i,3)*2-1)+RP(1);
    R(element(i,2)*2-1)=R(element(i,2)*2-1)+RP(1);
end

%% 固定对称面
%与y轴平行的对称面u为0，与x轴平行的对称面v为0
tty=zeros(num_node*2,1);
ttx=zeros(1,num_node*2);
%平行于y轴的对称面，u为0
for i=1:(ny+1)
   ii=(13+(nx+1)*(i-1))*2-1;
   %变换总刚
   KK(:,ii)=tty;
   KK(ii,:)=ttx;
   KK(ii,ii)=1;
   R(ii)=0;
end
%平行于x轴的对称面，v为0
for i=1:(nx+1)
   ii=(4*(nx+1)+i)*2;
   %变换总刚
   KK(:,ii)=tty;
   KK(ii,:)=ttx;
   KK(ii,ii)=1;
   R(ii)=0;
end

%% 求解位移
T=inv(KK)*R;
idex=[1 2 3];
U=zeros(num_node,2);
U=[T(1:2:num_node*2),T(2:2:num_node*2)];

%% 求解应力
%存储单元应力
Stress_all_element=zeros(num_element,3);
%存储总单元应力
Stress_all_node=zeros(num_node,2);
%求解各单元应力
for i=1:num_element
    %求解单个单元应力
    xt=[coordinate_node(element(i,1),1), coordinate_node(element(i,2),1), coordinate_node(element(i,3),1)];
    yt=[coordinate_node(element(i,1),2), coordinate_node(element(i,2),2), coordinate_node(element(i,3),2)];
    qt=[T(element(i,1)*2-1);T(element(i,1)*2);T(element(i,2)*2-1);T(element(i,2)*2);T(element(i,3)*2-1);T(element(i,3)*2)];
    Stresse=stress_element(Ys,poss,xt,yt,qt);
    %单个单元应力存储
    Stress_all_element(i,:)=Stresse';
end


%% 后处理
for i=1:num_element
     fill(coordinate_node(element(i,:),1),coordinate_node(element(i,:),2),Stress_all_element(element(i,:),1));
     hold on;
end
colorbar;
axis equal;
title('x方向应力云图');
ss=zeros(nx*2,1);
for i=1:nx*2
    ss(i)=Stress_all_element(1+(i-1)*8,1);
end
sst=zeros(nx,1);
for i=1:nx
    sst(i)=(ss(i*2)+ss(i*2-1))/2/1000000; 
end
xs=linspace(0.025,1.2-0.025,nx);
sstt=zeros(nx,1)+20;