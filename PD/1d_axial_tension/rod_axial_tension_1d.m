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
%考虑到左侧有位移约束，引入虚拟边界层，厚度为delta，右侧力作用处引入虚拟边界层，厚度dx，故总物质点数量+4
nx=L/dx+4;
%生成物质点坐标
coordinate_m=zeros(nx,1);
for i=1:nx
   coordinate_m(nx-i+1,1)=L-0.5*dx-dx*(i-2);
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
    if i<=1001
        no(i,2)=no(i,2)+3;
    else
        no(i,2)=no(i,2)+(1004-i);
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

%% 计算密度矩阵
%自适应动力松弛中dt一般为1
dt=1;
%5是安全因子
%减小安全因子以加快收敛
density=zeros(nx,1);
for i=1:nx
    density(i,1)=0.25*dt*dt*(2*A*delta)*c_pd/dx*5;
end


%% 初始化并设置边界条件
% 初始速度为0，初始位移为0，左边界位移为0（设置虚拟边界位移为0）
velocity=zeros(nx,1);
velocity_middle_old=zeros(nx,1);
velocity_middle=zeros(nx,1);
displacement=zeros(nx,1);

%% 施加载荷
% 右侧有200N的力，则体力密度为200/(1e-9)
force_body=zeros(nx,1);
force_body(nx,1)=200e9;

%% 开始迭代得到稳态解
%确定迭代步数
nt=12000;
%ax存储x方向加速度,force_pd存储每个物质点总pd力
ax=zeros(nx,1);
force_pd=zeros(nx,1);
force_pd_old=zeros(nx,1);
%u_draw用于画图
u_draw=zeros(nt,1);
cn1=0;
cn2=0;
cn=0;
for it=1:nt
    cn1=0;
    cn2=0;
    cn=0;
    %每次计算pd力先置0
    force_pd=zeros(nx,1);
    %输出计算进度
    if rem(it,100)==0
        fprintf("%d/%d\n", it, nt); 
        %下述语句用于在计算中监控
        %displacement_ana=200/A/E*coordinate_m(4:1003);
        %plot(coordinate_m(4:1003,1),displacement(4:1003),'g',coordinate_m(4:1003,1),displacement_ana,'r--','LineWidth',1.2);
        %grid on;
        %xlabel("x/m");ylabel("u/m");title("杆任意点位移的PD解与解析解对比");
        %legend("PD解","解析解",'Location','northwest');
        %getframe;
    end
    %进行时域积分，计算每个物质点的PD力
    for i=1:nx
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
    %自适应动态松弛
    %计算参数cn
    %先计算分子cn1
    %先计算刚度矩阵
    k=zeros(nx,nx);
    for i=4:nx
        if (velocity_middle_old(i, 1) ~= 0.0)
            cn1=cn1+displacement(i,1)*(force_pd_old(i,1)/density(i,1)-force_pd(i,1)/density(i,1))*displacement(i,1)/(dt*velocity_middle_old(i,1));
        end
        cn2=cn2+displacement(i,1)*displacement(i,1);
    end
    %cn1=displacement'*k*displacement;
    %再计算分母cn2
    %最后计算cn
    if (cn2 ~= 0.0)
        if ((cn1 / cn2) > 0.0)
            cn = 2.0 * sqrt(cn1/cn2);
        else
            cn = 0.0;
        end
    else
        cn = 0.0;
    end

    if (cn > 2.0)
        cn = 1.9;
    end
    %迭代得到速度、位移
    %先记录现在的速度与、位移、力
    for i=4:nx
        if it==1
            velocity_middle(i,1)=dt/density(i,1)*(force_pd(i,1)+force_body(i,1))/2;
        else
            velocity_middle(i,1)=((2-cn*dt)*velocity_middle_old(i,1)+2*dt/density(i,1)*(force_pd(i,1)+force_body(i,1)))/(2+cn*dt);
        end
        velocity(i,1)=0.5*(velocity_middle(i,1)+velocity_middle_old(i,1));
        displacement(i,1)=displacement(i,1)+velocity_middle(i,1)*dt;
        velocity_middle_old(i,1)=velocity_middle(i,1);
        force_pd_old(i,1)=force_pd(i,1);     
    end
    %固定边界条件
    displacement(1:3,1)=0;
    velocity(1:3,1)=0;
    velocity_middle(1:3,1)=0;
    velocity_middle_old(1:3,1)=0;
    %记录中点位移
    u_draw(it,1)=displacement(503,1);
end

%% 画图
%计算解析解
subplot(1,2,1);
displacement_ana=200/A/E*coordinate_m(4:1003);
plot(coordinate_m(4:1003,1),displacement(4:1003)*1000,'g',coordinate_m(4:1003,1),displacement_ana*1000,'r--','LineWidth',1.2);
grid on;
xlabel("x/m");ylabel("u/mm");title("杆任意点位移的PD解与解析解对比");
legend("PD解","解析解",'Location','northwest');
subplot(1,2,2);
plot(linspace(1,nt,nt),u_draw,'k','LineWidth',1.2);
grid on;
xlabel("iteration");ylabel("u/m");title("中点位置收敛情况");
