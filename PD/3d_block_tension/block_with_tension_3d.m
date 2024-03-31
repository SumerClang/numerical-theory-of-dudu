clc;clear;
%% 前处理，给定平面的物理参数
%给定杆件的长度、宽度、厚度、杨氏模量、泊松比、热膨胀系数、体积模量、剪切模量、密度，物质点间隔、横截面积、邻域作用范围,体积、都为标准国际单位
L=1;
W=0.1;
H=0.1;
E=200e9;
Poisson=0.25;
alpha=23e-6;
K=E/(3*(1-2*Poisson));
G=E/(2*(1+Poisson));
rou=7850;
dx=L/100;
dy=W/10;
dz=H/10;
A=dy*dx;
delta=3.015*dx;
dV=dx*dy*dz;

%% 生成质点并判断近场域内的物质点编号
%左侧增加虚拟边界，nx=L/dx+3,ny=W/dy,nz=H/dz
nx=L/dx+3;
ny=W/dy;
nz=H/dz;
number_of_point=nx*ny*nz;
%coordinate_x存储原始坐标
coordinate_x=zeros(number_of_point,3);
t=1;
%生成各个点的坐标
for i=1:nx
    for j=1:ny
        for k=1:nz 
            coordinate_x(t,1)=dx*(-2.5+(i-1));
            coordinate_x(t,2)=dy*(0.5+(j-1));
            coordinate_x(t,3)=dz*(0.5+(k-1));
            t=t+1;
        end
    end
end
%list_of_near_point存储相邻点的下标，最多只有122个相邻点
list_of_near_point=zeros(number_of_point,122);
%判断哪个点在近场域内
for i=1:number_of_point
    xi=coordinate_x(i,1);
    yi=coordinate_x(i,2);
    zi=coordinate_x(i,3);
    k=1;
    for j=1:number_of_point
        if j~=i
            xj=coordinate_x(j,1);
            yj=coordinate_x(j,2);
            zj=coordinate_x(j,3);
            distance=sqrt((xj-xi)^2+(yj-yi)^2+(zj-zi)^2);
            if distance<delta
                list_of_near_point(i,k)=j;
                k=k+1;
            end    
        end
    end    
end

%% 确定近场动力学参数
a_pd=0;
a2_pd=6*alpha*a_pd;
a3_pd=9*a_pd*alpha*alpha;
b_pd=15*G/(2*pi*delta^5);
c_pd=18*K/(pi*delta^4);
d_pd=9/(4*pi*delta^4);

%% 确定每个点的表面修正系数
%r用于确定体积修正系数
r=dx/2;
%假定三个方向的小位移用于确定表面修正因子,c_my存储两种变形后坐标
displacement1=zeros(number_of_point,3);
displacement1(:,1)=1e-3*coordinate_x(:,1);
displacement2=zeros(number_of_point,3);
displacement2(:,2)=1e-3*coordinate_x(:,2);
displacement3=zeros(number_of_point,3);
displacement3(:,3)=1e-3*coordinate_x(:,3);
coordinate_y1=coordinate_x+displacement1;
coordinate_y2=coordinate_x+displacement2;
coordinate_y3=coordinate_x+displacement3;
%存储应变能密度（连续介质力学和PD）
energy_strain_continuum=0.6*E*1e-6;
%energy_strain_continuum=1.5*E/(1-2*Poisson)*1e-6;
energy_strain_pd=zeros(number_of_point,3);
for i=1:number_of_point
    %近场域内最多有122个点
    for j=1:122
        %记录近场域内点的索引
        n=list_of_near_point(i,j);
        %只有存取的点不是0的时候才进入计算
        if n~=0
            xxi=coordinate_x(i,1);xyi=coordinate_x(i,2);xzi=coordinate_x(i,3);
            xxj=coordinate_x(n,1);xyj=coordinate_x(n,2);xzj=coordinate_x(n,3);
            yxi1=coordinate_y1(i,1);yyi1=coordinate_y1(i,2);yzi1=coordinate_y1(i,3);
            yxj1=coordinate_y1(n,1);yyj1=coordinate_y1(n,2);yzj1=coordinate_y1(n,3);
            yxi2=coordinate_y2(i,1);yyi2=coordinate_y2(i,2);yzi2=coordinate_y2(i,3);
            yxj2=coordinate_y2(n,1);yyj2=coordinate_y2(n,2);yzj2=coordinate_y2(n,3);
            yxi3=coordinate_y3(i,1);yyi3=coordinate_y3(i,2);yzi3=coordinate_y3(i,3);
            yxj3=coordinate_y3(n,1);yyj3=coordinate_y3(n,2);yzj3=coordinate_y3(n,3);
            %计算两个点之间的距离，以及变形后两个点之间的距离
            distancex=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi)+(xzj-xzi)*(xzj-xzi));
            distancey1=sqrt((yxj1-yxi1)*(yxj1-yxi1)+(yyj1-yyi1)*(yyj1-yyi1)+(yzj1-yzi1)*(yzj1-yzi1));
            distancey2=sqrt((yxj2-yxi2)*(yxj2-yxi2)+(yyj2-yyi2)*(yyj2-yyi2)+(yzj2-yzi2)*(yzj2-yzi2));
            distancey3=sqrt((yxj3-yxi3)*(yxj3-yxi3)+(yyj3-yyi3)*(yyj3-yyi3)+(yzj3-yzi3)*(yzj3-yzi3));
            %facv存储体积修正系数
            if distancex<=(delta-r)
                facv=1;
            elseif distancex<=(delta+r)
                facv=(delta+r-distancex)/(2*r);
            else
                facv=0;
            end
            %ce代表每次的应变能增量
            ce=zeros(1,3);
            ce(1,1)=0.25*c_pd/distancex*(distancex-distancey1)*(distancex-distancey1)*dV*facv;
            ce(1,2)=0.25*c_pd/distancex*(distancex-distancey2)*(distancex-distancey2)*dV*facv;
            ce(1,3)=0.25*c_pd/distancex*(distancex-distancey3)*(distancex-distancey3)*dV*facv;
            energy_strain_pd(i,:)=energy_strain_pd(i,:)+ce;
        else
            break;
        end
    end
end
%存储表面修正系数
fac=energy_strain_continuum./energy_strain_pd;

%% 计算密度矩阵
%自适应动力松弛中dt一般为1
dt=1;
%safe是安全因子，一般为5
%减小安全因子以加快收敛
safe=1;
density=zeros(number_of_point,3);
for i=1:number_of_point
    density(i,1)=0.25*dt*dt*(4/3*pi*delta^3)*c_pd/dx*safe;
    density(i,2)=0.25*dt*dt*(4/3*pi*delta^3)*c_pd/dx*safe;
    density(i,3)=0.25*dt*dt*(4/3*pi*delta^3)*c_pd/dx*safe;
end

%% 初始化并设置边界条件
velocity=zeros(number_of_point,3);
displacement=zeros(number_of_point,3);

%% 在右边界施加载荷
%载荷大小200MPa,体力密度2e10
force_body_rou=2e10;
%给予边界体力
force_body=zeros(number_of_point,3);
force_body(10201:10300,1)=force_body_rou;

%% 开始迭代得到稳态解
%确定迭代步数
nt=1500;
%force_pd存储每个物质点总pd力
force_pd=zeros(number_of_point,3);
force_pd_old=zeros(number_of_point,3);
%vmiddle和vmiddle_old存储计算过程中的中间量
velocity_middle=zeros(number_of_point,3);
velocity_middle_old=zeros(number_of_point,3);
%u_draw用于画图
u_draw=zeros(nt,3);
%计算阻尼系数，cn1为分子，cn2为分母，cn为阻尼系数，每个时间步开始前重置
cn1=0;
cn2=0;
cn=0;
for it=1:nt
    %cn1=0;
    %cn2=0;
    %cn=0;
    %每次计算前PD力先置零
    %计算每个点的PD力
    %进入循环前先先表示出变形后的坐标
    coordinate_y=coordinate_x+displacement;
    for i=1:number_of_point
        force_pd(i,:)=zeros(1,3);
        %临时的近场力置0
        t=zeros(1,3);
        for j=1:122
            %记录近场域内点的索引
            n=list_of_near_point(i,j);
            %只有存取的点不是0的时候才进入计算
            if n~=0
                xxi=coordinate_x(i,1);xyi=coordinate_x(i,2);xzi=coordinate_x(i,3);
                xxj=coordinate_x(n,1);xyj=coordinate_x(n,2);xzj=coordinate_x(n,3);
                yxi=coordinate_y(i,1);yyi=coordinate_y(i,2);yzi=coordinate_y(i,3);
                yxj=coordinate_y(n,1);yyj=coordinate_y(n,2);yzj=coordinate_y(n,3);
                %计算两个点之间的距离，以及变形后两个点之间的距离
                distancex=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi)+(xzj-xzi)*(xzj-xzi));
                distancey=sqrt((yxj-yxi)*(yxj-yxi)+(yyj-yyi)*(yyj-yyi)+(yzj-yzi)*(yzj-yzi));
                %facv存储体积修正系数
                if distancex<=(delta-r)
                    facv=1;
                elseif distancex<=(delta+r)
                    facv=(delta+r-distancex)/(2*r);
                else
                    facv=0;
                end
                %计算x、y、z方向的表面修正系数
                factx=0.5*(fac(i,1)+fac(n,1));
                facty=0.5*(fac(i,2)+fac(n,2));
                factz=0.5*(fac(i,3)+fac(n,3));
                %计算相对单位矢量
                e=[xxi-xxj,xyi-xyj,xzi-xzj];
                e=e/norm(e);
                %计算两个质点之间相互作用的表面修正系数
                fact=1/sqrt(e(1,1)*e(1,1)/(factx*factx)+e(1,2)*e(1,2)/(facty*facty)+e(1,3)*e(1,3)/(factz*factz));
                %计算kexi、yita、伸长率、以及质点之间的近场力
                %kexi=coordinate_x(n,:)-coordinate_x(i,:);
                %yita=displacement(n,:)-displacement(i,:);
                yj=[yxj,yyj,yzj];yi=[yxi,yyi,yzi];
                s=(distancey-distancex)/distancex;
                %t=t+2*(kexi+yita)/(norm(kexi+yita))*2*b_pd*delta*s*dV*facv*fact;
                t=t+c_pd*s*(yj-yi)/(norm(yj-yi))*dV*facv*fact;
            else
                break;
            end    
        end
        force_pd(i,:)=force_pd(i,:)+t;
    end
    %自适应动态松弛
    %计算参数cn
    %先计算分子cn1,再计算分母cn2
    for i=1:number_of_point
        if (velocity_middle_old(i, 1) ~= 0.0)
            cn1=cn1+displacement(i,1)*(force_pd_old(i,1)/density(i,1)-force_pd(i,1)/density(i,1))*displacement(i,1)/(dt*velocity_middle_old(i,1));
        end
        if (velocity_middle_old(i, 2) ~= 0.0)
            cn1=cn1+displacement(i,2)*(force_pd_old(i,2)/density(i,2)-force_pd(i,2)/density(i,2))*displacement(i,2)/(dt*velocity_middle_old(i,2));
        end
        if (velocity_middle_old(i, 3) ~= 0.0)
            cn1=cn1+displacement(i,3)*(force_pd_old(i,3)/density(i,3)-force_pd(i,3)/density(i,3))*displacement(i,3)/(dt*velocity_middle_old(i,3));
        end
        cn2=cn2+displacement(i,1)*displacement(i,1);
        cn2=cn2+displacement(i,2)*displacement(i,2);
        cn2=cn2+displacement(i,3)*displacement(i,3);  
    end
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
    %迭代计算速度，位移，加速度
    for i=1:number_of_point
        if it==1
            velocity_middle(i,1)=dt/density(i,1)*(force_pd(i,1)+force_body(i,1))/2;
            velocity_middle(i,2)=dt/density(i,2)*(force_pd(i,2)+force_body(i,2))/2;
            velocity_middle(i,3)=dt/density(i,3)*(force_pd(i,3)+force_body(i,3))/2;
        else
            velocity_middle(i,1)=((2-cn*dt)*velocity_middle_old(i,1)+2*dt/density(i,1)*(force_pd(i,1)+force_body(i,1)))/(2+cn*dt);
            velocity_middle(i,2)=((2-cn*dt)*velocity_middle_old(i,2)+2*dt/density(i,2)*(force_pd(i,2)+force_body(i,2)))/(2+cn*dt);
            velocity_middle(i,3)=((2-cn*dt)*velocity_middle_old(i,3)+2*dt/density(i,3)*(force_pd(i,3)+force_body(i,3)))/(2+cn*dt);
        end
        velocity(i,:)=0.5*(velocity_middle(i,:)+velocity_middle_old(i,:));
        %velocity(i,1)=0.5*(velocity_middle(i,1)+velocity_middle_old(i,1));
        %velocity(i,2)=0.5*(velocity_middle(i,2)+velocity_middle_old(i,2));
        %velocity(i,3)=0.5*(velocity_middle(i,3)+velocity_middle_old(i,3));
        displacement(i,:)=displacement(i,:)+velocity_middle(i,:)*dt;
        %displacement(i,1)=displacement(i,1)+velocity_middle(i,1)*dt;
        %displacement(i,2)=displacement(i,2)+velocity_middle(i,2)*dt;
        %displacement(i,3)=displacement(i,3)+velocity_middle(i,3)*dt;
        velocity_middle_old(i,:)=velocity_middle(i,:);
        %velocity_middle_old(i,1)=velocity_middle(i,1);
        %velocity_middle_old(i,2)=velocity_middle(i,2);
        %velocity_middle_old(i,3)=velocity_middle(i,3);
        force_pd_old(i,:)=force_pd(i,:);
        %force_pd_old(i,1)=force_pd(i,1);
        %force_pd_old(i,2)=force_pd(i,2);
        %force_pd_old(i,3)=force_pd(i,3);
    end
    %左边界位移和速度为0
    velocity(1:300,:)=zeros(300,3);
    displacement(1:300,:)=zeros(300,3);
    u_draw(it,:)=displacement(10223,:);
    fprintf("%d/%d\n", it, nt); 
end

%% 画图
%绘制x向中线位移ux
subplot(2,2,1);
cx=zeros(100,1);
ux=zeros(100,1);
for i=1:(nx-3)
    cx(i,1)=(0.5+(i-1))*dx;
end
for i=1:100
    ux(i,1)=0.25*(displacement(345+(i-1)*100,1)+displacement(346+(i-1)*100,1)+displacement(355+(i-1)*100,1)+displacement(356+(i-1)*100,1)); 
end
uxa=1e-3*cx;
plot(cx,1000*ux,'g',cx,1000*uxa,'r--','LineWidth',1.2);
xlabel("x/m");ylabel("ux/mm");
legend("PD","解析解",'Location','northwest');
grid on;
title("x方向中线位移ux");

%绘制y向中线位移uy
subplot(2,2,2);
cy=zeros(10,1);
uy=zeros(10,1);
for i=1:10
    cy(i,1)=(0.5+(i-1))*dy;
end
uya=-Poisson*1e-3*(cy-W/2);
for i=1:10
    uy(i,1)=0.25*(displacement(5205+(i-1)*10,2)+displacement(5206+(i-1)*10,2)+displacement(5305+(i-1)*10,2)+displacement(5306+(i-1)*10,2)); 
end
plot(cy,1000*uy,'g',cy,1000*uya,'r--','LineWidth',1.2);
xlabel("y/m");ylabel("uy/mm");
legend("PD","解析解",'Location','northeast');
grid on;
title("y方向中线位移uy");

%绘制z向中线位移uz
subplot(2,2,3);
cz=zeros(10,1);
uz=zeros(10,1);
for i=1:10
    cz(i,1)=(0.5+(i-1))*dz;
end
uza=-Poisson*1e-3*(cy-H/2);
for i=1:10
    uz(i,1)=0.25*(displacement(5241+(i-1),3)+displacement(5251+(i-1),3)+displacement(5341+(i-1),3)+displacement(5351+(i-1),3)); 
end
plot(cz,1000*uz,'g',cz,1000*uza,'r--','LineWidth',1.2);
xlabel("z/m");ylabel("uz/mm");
legend("PD","解析解",'Location','northeast');
grid on;
title("z方向中线位移uz");

%绘制点(0.995,0.025,0.025)收敛情况
subplot(2,2,4);
ite=linspace(1,nt,nt);
plot(ite,u_draw(:,1)*1000,'r',ite,u_draw(:,2)*1000,'g',ite,u_draw(:,3)*1000,'b','LineWidth',1.2);
xlabel("iteration");ylabel("u/mm");
legend("ux","uy","uz",'Location','northeast');
title("点(0.995,0.025,0.025)的收敛情况");
grid on;