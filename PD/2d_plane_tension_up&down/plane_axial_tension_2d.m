clc;clear;
%% 前处理，给定平面的物理参数
%给定杆件的长度、宽度、厚度、杨氏模量、泊松比、热膨胀系数、体积模量、剪切模量、密度，物质点间隔、横截面积、邻域作用范围,体积、都为标准国际单位
L=1;
W=0.5;
H=0.01;
E=200e9;
Poisson=1/3;
alpha=23e-6;
K=E/(2*(1-Poisson));
G=E/(2*(1+Poisson));
rou=7850;
dx=L/100;
dy=W/50;
A=dy*dx;
delta=3.015*dx;
dV=dx*dy*H;

%% 生成质点并判断近场域内的物质点编号
% 上下各加一层虚拟边界，x方向有100个物质点，y方向有52个物质点
nx=L/dx;ny=W/dy+2;
nom=nx*ny;
%第一列存储x坐标，第二列存储y坐标
coordinate_mx=zeros(nom,2);
for i=1:ny
    coordinate_mx((1+nx*(i-1)):(nx*i),1)=linspace(dx/2,L-dx/2,nx);
    coordinate_mx((1+nx*(i-1)):(nx*i),2)=ones(nx,1)*dy*i-dy*1.5;    
end
%用coordinate_mpd存储近场域内的质点的索引，该问题中一个质点近场域内最多只有7*7-1-4*5=28个质点
coordinate_mpd=zeros(nx*ny,28);
%对每一个质点开始遍历近场域内的质点下标
for i=1:nom
    xxi=coordinate_mx(i,1);
    xyi=coordinate_mx(i,2);
    k=1;
    for j=1:nx*ny
        if j~=i
            xxj=coordinate_mx(j,1);
            xyj=coordinate_mx(j,2);
            distance=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi));
            if (distance)<delta
                coordinate_mpd(i,k)=j;
                k=k+1;
            end
        end
    end
end

%% 确定近场动力学参数
a_pd=0;
a2_pd=4*alpha*a_pd;
a3_pd=4*a_pd*alpha*alpha;
b_pd=6*G/(pi*H*delta^4);
c_pd=24*G/(pi*H*delta^3);
d_pd=2/(pi*H*delta^3);

%% 确定每个物质点的表面修正系数
%r用于确定体积修正系数
r=dx/2;
%假定两个方向的小位移用于确定表面修正因子,c_my存储两种变形后坐标
displacement1=zeros(nx*ny,2);
displacement1(:,1)=1e-3*coordinate_mx(:,1);
displacement2=zeros(nx*ny,2);
displacement2(:,2)=1e-3*coordinate_mx(:,2);
coordinate_my1=coordinate_mx+displacement1;
coordinate_my2=coordinate_mx+displacement2;
%存储应变能密度（连续介质力学和PD）
energy_strain_continuum=9/16*E*1e-6;
energy_strain_pd=zeros(nx*ny,2);
for i=1:nom
    for j=1:28
        %记录近场域内点的索引
        n=coordinate_mpd(i,j);
        %只有存取的点不是0的时候才进入计算
        if n~=0
            xxi=coordinate_mx(i,1);xyi=coordinate_mx(i,2);
            xxj=coordinate_mx(n,1);xyj=coordinate_mx(n,2);
            yxi1=coordinate_my1(i,1);yyi1=coordinate_my1(i,2);
            yxj1=coordinate_my1(n,1);yyj1=coordinate_my1(n,2);
            yxi2=coordinate_my2(i,1);yyi2=coordinate_my2(i,2);
            yxj2=coordinate_my2(n,1);yyj2=coordinate_my2(n,2);
            %计算两个点之间的距离，以及变形后两个点之间的距离
            distancex=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi));
            distancey1=sqrt((yxj1-yxi1)*(yxj1-yxi1)+(yyj1-yyi1)*(yyj1-yyi1));
            distancey2=sqrt((yxj2-yxi2)*(yxj2-yxi2)+(yyj2-yyi2)*(yyj2-yyi2));
            %facv存储体积修正系数
            if distancex<=(delta-r)
                facv=1;
            elseif distancex<=(delta+r)
                facv=(delta+r-distancex)/(2*r);
            else
                facv=0;
            end
            %ce代表每次的应变能增量
            ce=zeros(1,2);
            ce(1,1)=b_pd*delta/distancex*(distancex-distancey1)*(distancex-distancey1)*dV*facv;
            ce(1,2)=b_pd*delta/distancex*(distancex-distancey2)*(distancex-distancey2)*dV*facv;
            energy_strain_pd(i,:)=energy_strain_pd(i,:)+ce;
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
safe=5;
density=zeros(nom,2);
for i=1:nom
    density(i,1)=0.25*dt*dt*(pi*delta^2*H)*c_pd/dx*safe;
    density(i,2)=0.25*dt*dt*(pi*delta^2*H)*c_pd/dx*safe;
end

%% 初始化并设置边界条件
velocity=zeros(nom,2);
displacement=zeros(nom,2);

%% 在上下边界施加载荷
%载荷大小200MPa,体力密度2e10
force_body_rou=2e10;
%给予边界体力
force_body=zeros(nom,2);
force_body(1:nx,2)=-1*force_body_rou*ones(nx,1);
force_body((nx*ny-nx+1):nx*ny,2)=force_body_rou*ones(nx,1);

%% 开始迭代得到稳态解
%确定迭代步数
nt=1500;
%force_pd存储每个物质点总pd力
force_pd=zeros(nx*ny,2);
force_pd_old=zeros(nx*ny,2);
%vmiddle和vmiddle_old存储计算过程中的中间量
velocity_middle=zeros(nx*ny,2);
velocity_middle_old=zeros(nx*ny,2);
%u_draw用于画图
u_draw=zeros(nt,2);
%计算阻尼系数，cn1为分子，cn2为分母，cn为阻尼系数，每个时间步开始前重置
cn1=0;
cn2=0;
cn=0;
for it=1:nt
    %cn1=0;
    %cn2=0;
    %cn=0;
    if rem(it,100)==0
        %下述语句用于在计算中监控
        %plot(linspace(1,it-1,it-1),-2.05e-4*ones(it-1,1),'g',linspace(1,it-1,it-1),u_draw(1:(it-1),2),'r','LineWidth',1.2);
        %legend("uy_{ana}","uy",'Location','northeast');
        %x=linspace(0.005,0.995,100);ux=-Poisson*1e-3*(x-0.5);
        %plot(x,ux,'r--',x,displacement(2501:2600,1),'g','LineWidth',1.2);
        %y=linspace(0.005,0.495,50);uy_ana=1e-3*(y'-0.25);
        %uy=zeros(1,50);
        %for i=1:50
        %    uy(:,i)=(displacement(150+(i-1)*100,2)+displacement(151+(i-1)*100,2))/2; 
        %end
        %subplot(2,2,1);
        %plot(y,uy_ana,'r--',y,uy,'g','LineWidth',1.2)
        %xlabel("y/m");ylabel("uy/m");
        %legend("解析解","PD解",'Location','northwest');
        %title("纵向中线位移uy");
        %grid on;
        %subplot(2,2,2);
        %x=linspace(0.005,0.995,100);ux=-Poisson*1e-3*(x-0.5);
        %plot(x,ux,'r--',x,0.5*(displacement(2501:2600,1)+displacement(2601:2700,1)),'g','LineWidth',1.2);
        %xlabel("x/m");ylabel("ux/m");
        %legend("解析解","PD解",'Location','northeast');
        %title("横向中线位移ux");
        %grid on;
        %getframe;
    end
    %每次计算前PD力先置零
    %计算每个点的PD力
    %进入循环前先先表示出变形后的坐标
    coordinate_my=coordinate_mx+displacement;
    for i=1:nom
        force_pd(i,:)=zeros(1,2);
        %临时的近场力置0
        t=zeros(1,2);
        for j=1:28
            %记录近场域内点的索引
            n=coordinate_mpd(i,j);
            %只有存取的点不是0的时候才进入计算
            if n~=0
                xxi=coordinate_mx(i,1);xyi=coordinate_mx(i,2);
                xxj=coordinate_mx(n,1);xyj=coordinate_mx(n,2);
                yxi=coordinate_my(i,1);yyi=coordinate_my(i,2);
                yxj=coordinate_my(n,1);yyj=coordinate_my(n,2);
                %计算两个点之间的距离，以及变形后两个点之间的距离
                distancex=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi));
                distancey=sqrt((yxj-yxi)*(yxj-yxi)+(yyj-yyi)*(yyj-yyi));
                %facv存储体积修正系数
                if distancex<=(delta-r)
                    facv=1;
                elseif distancex<=(delta+r)
                    facv=(delta+r-distancex)/(2*r);
                else
                    facv=0;
                end
                %计算x和y方向的表面修正系数
                factx=0.5*(fac(i,1)+fac(n,1));
                facty=0.5*(fac(i,2)+fac(n,2));
                %计算相对单位矢量
                e=[xxi-xxj,xyi-xyj];
                e=e/norm(e);
                %计算两个质点之间相互作用的表面修正系数
                fact=1/sqrt(e(1,1)*e(1,1)/(factx*factx)+e(1,2)*e(1,2)/(facty*facty));
                %计算kexi、yita、伸长率、以及质点之间的近场力
                kexi=coordinate_mx(n,:)-coordinate_mx(i,:);
                yita=displacement(n,:)-displacement(i,:);
                s=(distancey-distancex)/distancex;
                t=t+2*(kexi+yita)/(norm(kexi+yita))*2*b_pd*delta*s*dV*facv*fact;
                %flag=zeros(1,2);
                %if yxj<yxi
                %    flag(1,1)=-1;
                %else
                %    flag(1,1)=1;
                %end
                %if yyj<yyi
                %    flag(1,2)=-1;
                %else
                %    flag(1,2)=1;
                %end
                %t1=c_pd*s*dV*fact*facv*flag(1,1);
                %t2=c_pd*s*dV*fact*facv*flag(1,2);
                %t(1,1)=t(1,1)+t1;
                %t(1,2)=t(1,2)+t2;
            end    
        end
        force_pd(i,:)=force_pd(i,:)+t;
    end
    %自适应动态松弛
    %计算参数cn
    %先计算分子cn1,再计算分母cn2
    for i=1:nom
        if (velocity_middle_old(i, 1) ~= 0.0)
            cn1=cn1+displacement(i,1)*(force_pd_old(i,1)/density(i,1)-force_pd(i,1)/density(i,1))*displacement(i,1)/(dt*velocity_middle_old(i,1));
        end
        if (velocity_middle_old(i, 2) ~= 0.0)
            cn1=cn1+displacement(i,2)*(force_pd_old(i,2)/density(i,2)-force_pd(i,2)/density(i,2))*displacement(i,2)/(dt*velocity_middle_old(i,2));
        end
        cn2=cn2+displacement(i,1)*displacement(i,1);
        cn2=cn2+displacement(i,2)*displacement(i,2);        
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
    for i=1:nx*ny
        if it==1
            velocity_middle(i,1)=dt/density(i,1)*(force_pd(i,1)+force_body(i,1))/2;
            velocity_middle(i,2)=dt/density(i,2)*(force_pd(i,2)+force_body(i,2))/2;
        else
            velocity_middle(i,1)=((2-cn*dt)*velocity_middle_old(i,1)+2*dt/density(i,1)*(force_pd(i,1)+force_body(i,1)))/(2+cn*dt);
            velocity_middle(i,2)=((2-cn*dt)*velocity_middle_old(i,2)+2*dt/density(i,2)*(force_pd(i,2)+force_body(i,2)))/(2+cn*dt);
        end
        velocity(i,1)=0.5*(velocity_middle(i,1)+velocity_middle_old(i,1));
        velocity(i,2)=0.5*(velocity_middle(i,2)+velocity_middle_old(i,2));
        displacement(i,1)=displacement(i,1)+velocity_middle(i,1)*dt;
        displacement(i,2)=displacement(i,2)+velocity_middle(i,2)*dt;
        velocity_middle_old(i,1)=velocity_middle(i,1);
        velocity_middle_old(i,2)=velocity_middle(i,2);
        force_pd_old(i,1)=force_pd(i,1);
        force_pd_old(i,2)=force_pd(i,2);
    end
    u_draw(it,:)=displacement(503,:);
    fprintf("%d/%d\n", it, nt); 
end

%% 画图
%任意处uy=1e-3(y-0.25),ux=-poisson*1e-3*(x-0.5)
y=linspace(0.005,0.495,50);uy_ana=1e-3*(y'-0.25);
uy=zeros(1,50);
for i=1:50
    uy(:,i)=(displacement(150+(i-1)*100,2)+displacement(151+(i-1)*100,2))/2; 
end
subplot(1,3,1);
plot(y,uy_ana*1000,'g',y,uy*1000,'r--','LineWidth',1.2)
xlabel("y/m");ylabel("uy/mm");
legend("解析解","PD解",'Location','northwest');
title("纵向中线位移uy");
grid on;
subplot(1,3,2);
x=linspace(0.005,0.995,100);ux=-Poisson*1e-3*(x-0.5);
plot(x,ux*1000,'g',x,1000*0.5*(displacement(2501:2600,1)+displacement(2601:2700,1)),'r--','LineWidth',1.2);
xlabel("x/m");ylabel("ux/mm");
legend("解析解","PD解",'Location','northeast');
title("横向中线位移ux");
grid on;
subplot(1,3,3)
ite=linspace(1,1500,1500);
plot(ite,u_draw(:,1)*1000,'r',ite,u_draw(:,2)*1000,'g','LineWidth',1.2);
xlabel("iteration");ylabel("u/mm");
legend("ux","uy",'Location','northeast');
title("选择的点的收敛情况");
grid on;