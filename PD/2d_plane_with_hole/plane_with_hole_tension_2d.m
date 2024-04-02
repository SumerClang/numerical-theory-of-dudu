clc;clear;
%% 前处理，给定平面的物理参数
%给定杆件的长度、宽度、厚度、中心圆直径、杨氏模量、泊松比、热膨胀系数、体积模量、剪切模量、密度，物质点间隔、横截面积、邻域作用范围,体积、都为标准国际单位
L=50/1000;
W=50/1000;
H=0.5/1000;
D=10/1000;
E=192e9;
Poisson=1/3;
alpha=23e-6;
K=E/(2*(1-Poisson));
G=E/(2*(1+Poisson));
rou=8000;
dx=L/100;
dy=W/100;
A=dy*dx;
delta=3.015*dx;
dV=dx*dy*H;

%% 生成质点并判断近场域内的物质点编号
% 上下为位移边界，故各加一层虚拟边界，x方向有100个物质点，y方向有100+3+3个物质点
nx=L/dx;ny=W/dy+6;
nom=nx*ny;
%先通过tx存储临时坐标，该坐标系中未挖去中心的圆
%第一列存储x坐标，第二列存储y坐标，第三列存储是否需要去掉该点,第三列为1时代表在圆外
number_of_point=nx*ny;
coordinate_tx=zeros(number_of_point,3);
t=1;
%生成tx中各个点的坐标
coordinate_tx(:,3)=ones(nx*ny,1);
for i=1:ny
    for j=1:nx
        coordinate_tx(t,1)=(0.5+(j-1))*dx;
        coordinate_tx(t,2)=(-2.5+(i-1))*dy;
        t=t+1;
    end
end
%判断该点到原点的距离，距离小于D/2则需要除去
for i=1:number_of_point
    tx=coordinate_tx(i,1);
    ty=coordinate_tx(i,2);
    tr=sqrt((tx-25/1000)*(tx-25/1000)+(ty-25/1000)*(ty-25/1000));
    if tr<(D/2)
        coordinate_tx(i,3)=0; 
    end
end
%总质点数等于第三列非0的个数
coordinate_x=zeros(sum(coordinate_tx(:,3)),2);
t=1;
%开始为coordinate_x赋值坐标
for i=1:number_of_point
    if coordinate_tx(i,3)~=0
        coordinate_x(t,1)=coordinate_tx(i,1);
        coordinate_x(t,2)=coordinate_tx(i,2); 
        t=t+1;
    end
end
%总质点数等于第三列非0的个数
number_of_point=sum(coordinate_tx(:,3));
%list_of_near_point存储相邻点的下标，最多只有28个相邻点，mu存储相互作用是否断开,1表示未断开
list_of_near_point=zeros(number_of_point,28);
mu=zeros(number_of_point,28);
%开始遍历寻找下标
for i=1:number_of_point
    xi=coordinate_x(i,1);
    yi=coordinate_x(i,2);
    k=1;
   for j=1:number_of_point
       if j~=i
            xj=coordinate_x(j,1);
            yj=coordinate_x(j,2);
            distance=sqrt((xj-xi)^2+(yj-yi)^2);
            if (distance)<delta
                list_of_near_point(i,k)=j;
                mu(i,k)=1;
                k=k+1;
            end
       end
   end
end
damage=zeros(number_of_point,1);

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
displacement1=zeros(number_of_point,2);
displacement1(:,1)=1e-3*coordinate_x(:,1);
displacement2=zeros(number_of_point,2);
displacement2(:,2)=1e-3*coordinate_x(:,2);
coordinate_y1=coordinate_x+displacement1;
coordinate_y2=coordinate_x+displacement2;
%存储应变能密度（连续介质力学和PD）
energy_strain_continuum=9/16*E*1e-6;
energy_strain_pd=zeros(number_of_point,2);
for i=1:number_of_point
    for j=1:28
        %记录近场域内点的索引
        n=list_of_near_point(i,j);
        %只有存取的点不是0的时候才进入计算
        if n~=0
            xxi=coordinate_x(i,1);xyi=coordinate_x(i,2);
            xxj=coordinate_x(n,1);xyj=coordinate_x(n,2);
            yxi1=coordinate_y1(i,1);yyi1=coordinate_y1(i,2);
            yxj1=coordinate_y1(n,1);yyj1=coordinate_y1(n,2);
            yxi2=coordinate_y2(i,1);yyi2=coordinate_y2(i,2);
            yxj2=coordinate_y2(n,1);yyj2=coordinate_y2(n,2);
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
density=zeros(number_of_point,2);
for i=1:nom
    density(i,1)=0.25*dt*dt*(pi*delta^2*H)*c_pd/dx*safe;
    density(i,2)=0.25*dt*dt*(pi*delta^2*H)*c_pd/dx*safe;
end

%% 初始化并设置边界条件
velocity=zeros(number_of_point,2);
displacement=zeros(number_of_point,2);

%% 在右边界施加载荷
%载荷大小200MPa,体力密度2e10
%force_body_rou=2e10;
%给予边界体力
force_body=zeros(number_of_point,2);

%% 开始迭代得到稳态解
%确定迭代步数
nt=1500;
%force_pd存储每个物质点总pd力
force_pd=zeros(number_of_point,2);
force_pd_old=zeros(number_of_point,2);
%vmiddle和vmiddle_old存储计算过程中的中间量
velocity_middle=zeros(number_of_point,2);
velocity_middle_old=zeros(number_of_point,2);
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
    %每次计算前PD力先置零
    %计算每个点的PD力
    %进入循环前先先表示出变形后的坐标
    coordinate_y=coordinate_x+displacement;
    for i=1:number_of_point
        force_pd(i,:)=zeros(1,2);
        %临时的近场力置0
        t=zeros(1,2);
        for j=1:28
            %记录近场域内点的索引
            n=list_of_near_point(i,j);
            %只有存取的点不是0的时候才进入计算
            if n~=0
                xxi=coordinate_x(i,1);xyi=coordinate_x(i,2);
                xxj=coordinate_x(n,1);xyj=coordinate_x(n,2);
                yxi=coordinate_y(i,1);yyi=coordinate_y(i,2);
                yxj=coordinate_y(n,1);yyj=coordinate_y(n,2);
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
                %计算x、y、z方向的表面修正系数
                factx=0.5*(fac(i,1)+fac(n,1));
                facty=0.5*(fac(i,2)+fac(n,2));
                %计算相对单位矢量
                e=[xxi-xxj,xyi-xyj];
                e=e/norm(e);
                %计算两个质点之间相互作用的表面修正系数
                fact=1/sqrt(e(1,1)*e(1,1)/(factx*factx)+e(1,2)*e(1,2)/(facty*facty));
                %计算kexi、yita、伸长率、以及质点之间的近场力
                %kexi=coordinate_x(n,:)-coordinate_x(i,:);
                %yita=displacement(n,:)-displacement(i,:);
                yj=[yxj,yyj];yi=[yxi,yyi];
                s=(distancey-distancex)/distancex;
                if s>=0.02
                    mu(i,j)=0;
                end
                t=t+c_pd*s*(yj-yi)/(norm(yj-yi))*dV*facv*fact*mu(i,j);
            else
                break;
            end    
        end
        force_pd(i,:)=force_pd(i,:)+t;
        %计算局部损伤
        sum1=0;
        sum2=0;
        for j=1:28
            sum1=sum1+mu(i,j);
            if list_of_near_point(i,j)~=0
               sum2=sum2+1; 
            end
        end
        damage(i,1)=1-sum1/sum2;
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
    for i=1:number_of_point
        %固定左侧和右侧边界的速度
        if it==1
            velocity_middle(i,1)=dt/density(i,1)*(force_pd(i,1)+force_body(i,1))/2;
            velocity_middle(i,2)=dt/density(i,2)*(force_pd(i,2)+force_body(i,2))/2;
        else
            velocity_middle(i,1)=((2-cn*dt)*velocity_middle_old(i,1)+2*dt/density(i,1)*(force_pd(i,1)+force_body(i,1)))/(2+cn*dt);
            velocity_middle(i,2)=((2-cn*dt)*velocity_middle_old(i,2)+2*dt/density(i,2)*(force_pd(i,2)+force_body(i,2)))/(2+cn*dt);
        end
    end
    velocity_middle(1:300,2)=-2.7541e-7;
    velocity_middle((number_of_point-300+1):number_of_point,2)=2.7541e-7;
    for i=1:number_of_point
        velocity(i,:)=0.5*(velocity_middle(i,:)+velocity_middle_old(i,:));
        displacement(i,:)=displacement(i,:)+velocity_middle(i,:)*dt;
        velocity_middle_old(i,:)=velocity_middle(i,:);
        force_pd_old(i,:)=force_pd(i,:);        
    end
    u_draw(it,:)=displacement(5000,:);
    fprintf("%d/%d\n", it, nt);
    if rem(it,100)==0
        scale = 1;
        scatter(coordinate_x(:,1)+scale*displacement(:,1), coordinate_x(:,2)+scale*displacement(:,2),10, damage, "filled")
        grid on;
        colormap jet;
        colorbar;
        getframe;
    end
end