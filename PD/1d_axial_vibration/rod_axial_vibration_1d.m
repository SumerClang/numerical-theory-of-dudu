clc;clear;
%% 前处理
%给定杆件的长度、横截面积、杨氏模量、泊松比、剪切模量、密度，物质点间隔、邻域作用范围,体积、都为标准国际单位
L=1;
A=1e-6;
E=200e9;
Poisson=0.25;
G=E/(2*(1+Poisson));
rou=7850;
dx=0.001;
delta=3.015*dx;
dV=dx*A;
%考虑到左侧有位移约束，引入虚拟边界层，厚度为delta，故总物质点数量+3
nx=L/dx+3;
%生成物质点坐标
coordinate_m=zeros(nx,1);
for i=1:nx
   coordinate_m(nx-(i-1),1)=L-0.5*dx-dx*(i-1);
end

%% 确定每个质点近场域内的其他质点
%确定每个质点近场域内其他质点的数量
no=zeros(nx,2);
for i=1:nx
    %判断该质点左边有几个点
    if i>3
        no(i,1)=3;
    else
        no(i,1)=i-1;
    end
    %判断该质点右边有几个点
    if i<=1000
        no(i,2)=no(i,2)+3;
    else
        no(i,2)=no(i,2)+(1003-i);
    end
end

%% 确定近场动力学参数
a_pd=0;
a2_pd=0;
a3_pd=0;
b_pd=E/(2*A*delta*delta*delta);
d_pd=1/(2*delta*delta*A);
c_pd=2*E/(A*delta*delta);

%% 确定每个物质点的表面修正系数
%假定一个单向拉伸后确定其位移和应变能密度
%假定位移梯度kexi=1e-3，1d问题中弹性力学中应变能密度为W=0.5*E*epsil*epsilon,pd应变能密度待求
displacement=1e-3*coordinate_m;
coordinate_my=coordinate_m+displacement;
energy_strain_continuum=0.5*E*1e-6;
energy_strain_pd=zeros(nx,1);
%计算得到各个质点的表面修正系数
%facv为体积修正因子,r=dx/2
r=dx/2;
for i=1:nx
    %先算左侧近场域
    if no(i,1)>0
        for j=1:no(i,1)
            if abs(coordinate_m(i-j,1)-coordinate_m(i,1))<=(delta-r)
                facv=1;
            elseif abs(coordinate_m(i-j,1)-coordinate_m(i,1))<=(delta+r)
                facv=(delta+r-abs(coordinate_m(i-j,1)-coordinate_m(i,1)))/(2*r);
            else
                facv=0;
            end
            %计算x坐标差值，y坐标差值，PD应变能
            cx=abs(coordinate_m(i-j,1)-coordinate_m(i,1));
            cy=abs(coordinate_my(i-j,1)-coordinate_my(i,1));
            ce=b_pd*delta/cx*(cy-cx)*(cy-cx)*dV*facv;
            energy_strain_pd(i,1)=energy_strain_pd(i,1)+ce;
        end 
    end
    %再算右侧近场域
    if no(i,2)>0
       for j=1:no(i,2)
            if abs(coordinate_m(i+j,1)-coordinate_m(i,1))<=(delta-r)
                facv=1;
            elseif abs(coordinate_m(i+j,1)-coordinate_m(i,1))<=(delta+r)
                facv=(delta+r-abs(coordinate_m(i+j,1)-coordinate_m(i,1)))/(2*r);
            else
                facv=0;
            end
            %计算x坐标差值，y坐标差值，PD应变能
            cx=abs(coordinate_m(i+j,1)-coordinate_m(i,1));
            cy=abs(coordinate_my(i+j,1)-coordinate_my(i,1));
            ce=b_pd*delta/cx*(cy-cx)*(cy-cx)*dV*facv;
            energy_strain_pd(i,1)=energy_strain_pd(i,1)+ce;
       end
    end
end
% factor用于记录表面修正系数
factor=energy_strain_continuum./energy_strain_pd;

%% 初始化并设置边界条件
% 初始速度为0，设置位移梯度，左边界位移为0（设置虚拟边界位移为0）
velocity=zeros(nx,1);
displacement=1e-3*coordinate_m;
displacement(1:3,1)=0;

%% 开始时域积分
%确定时间步数和时间步长
nt=26000;
dt=0.8 * sqrt(2.0*rou*dx/(2.0 * delta * A * c_pd));
%ax存储x方向加速度,force_pd存储每个物质点总pd力
ax=zeros(nx,1);
force_pd=zeros(nx,1);
%d_draw用于绘制中心点的位移曲线
d_draw=zeros(nt,1);
%a_draw用于时域积分途中计算解析解
a_draw=zeros(nt,1);
for it=1:nt
    %每次计算pd力先置0
    force_pd=zeros(nx,1);
    %输出计算进度
    if rem(it,100)==0
        fprintf("%d/%d\n", it, nt); 
    end
    %进行时域积分，计算每个物质点的PD力
    for i=4:nx
    %每次计算前将t置0
        t=0;
    %s为伸长率。kexi为相对位置，yita为相对位移矢量
    %计算第i个物质点所受到的PD力
    %先算左侧近场域
        if no(i,1)>0
            for j=1:no(i,1)
                if abs(coordinate_m(i-j,1)-coordinate_m(i,1))<=(delta-r)
                    facv=1;
                elseif abs(coordinate_m(i-j,1)-coordinate_m(i,1))<=(delta+r)
                    facv=(delta+r-abs(coordinate_m(i-j,1)-coordinate_m(i,1)))/(2*r);
                else
                    facv=0;
                end
                %计算kexi、yita、伸长率、
                kexi=coordinate_m(i-j,1)-coordinate_m(i,1);
                yita=displacement(i-j,1)-displacement(i,1);
                s=(abs(kexi+yita)-abs(kexi))/(abs(kexi));
                fact=(factor(i)+factor(i-j))/2;
                t=t+2*(kexi+yita)/abs(kexi+yita)*2*b_pd*delta*s*dV*facv*fact;
            end 
        end
        %再算右侧近场域
        if no(i,2)>0
            for j=1:no(i,2)
                if abs(coordinate_m(i+j,1)-coordinate_m(i,1))<=(delta-r)
                    facv=1;
                elseif abs(coordinate_m(i+j,1)-coordinate_m(i,1))<=(delta+r)
                    facv=(delta+r-abs(coordinate_m(i+j,1)-coordinate_m(i,1)))/(2*r);
                else
                    facv=0;
                end
                %计算kexi、yita、伸长率、
                kexi=coordinate_m(i+j,1)-coordinate_m(i,1);
                yita=displacement(i+j,1)-displacement(i,1);
                s=(abs(kexi+yita)-abs(kexi))/(abs(kexi));
                fact=(factor(i)+factor(i+j))/2;
                t=t+2*(kexi+yita)/abs(kexi+yita)*2*b_pd*delta*s*dV*facv*fact;
            end 
        end
        force_pd(i,1)=force_pd(i,1)+t;
    end
    %得到加速度，并迭代
    ax=force_pd/rou;
    velocity=velocity+ax*dt;
    displacement=displacement+velocity*dt;
    d_draw(it,1)=displacement(503,1);
end

%% 画图
%给定时间
t_draw=linspace(1,26000,26000)*dt;
%画图
plot(t_draw*1000,d_draw,'LineWidth',1.2);
xlabel("t/\mus");ylabel("u/m");grid on;title("中心点位移随时间变化");