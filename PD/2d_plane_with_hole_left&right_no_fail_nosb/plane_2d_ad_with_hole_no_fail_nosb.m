clear;clc;
close all;
format long;
%addpath 'D:\code\matlab\bin\math theory\PD\state_base\brazil_plate\common_functions';
%公式来源可以参考https://doi.org/10.1007/s10409-019-00873-y
%% 根据算例具体信息在basic_information中设置并导入
%导入时间步长，划分密度，算例维度，时间步长，时间步数，热膨胀系数，杨氏模量，体积模量，剪切模量，材料矩阵，近场域大小
[basic,geometry,material] = basic_information();
dt=basic.dt;
dx=geometry.dx;
dy=geometry.dy;
Thick=geometry.Thick;
dV=geometry.dV;
dimen=basic.dimen;
dt=basic.dt;
nt=basic.nt;
alpha=material.alpha;
E=material.E;
K=material.K;
G=material.G;
%sigma=C*epsilon，二维问题中C为3*3时，epsilon需要为3*1，排列为[xx，yy，xy]
%此时得到的sigma也为3*1，排列为[xx，yy，xy]
%需要注意[xx yy xy]和[xx xy;xy yy]之间的转换
C=material.C;
delta = basic.delta;

%% 生成质点坐标
[coordinate_x,number_of_point]=generate_point(geometry,basic);

%% 根据实际情况设置边界条件
[coordinate_x,number_of_point,force_body,node_left,node_right]=determine_boundary(coordinate_x,number_of_point,basic,dimen);

%% 判断质点近场域内的其他点，并初步设置相互作用
[list_of_near_point,mu,damage] = determine_near_point(coordinate_x,number_of_point,delta);

%% 确定近场动力学参数
a_pd=0.5*(K-2*G);
a2_pd=4*alpha*a_pd;
a3_pd=4*a_pd*alpha*alpha;
b_pd=6*G/(pi*Thick*delta^4);
c_pd=24*G/(pi*Thick*delta^3);
d_pd=2/(pi*Thick*delta^3);

%% 计算自适应动态松弛的密度矢量
%safe是密度矢量的安全因子
safe=1.5;
[density] = generate_massvector(number_of_point,dt,delta,Thick,dx,safe,dimen,c_pd);

%% 根据具体情况初始化位移，速度，加速度
[displacement,velocity,acceleration] = initial_u_v_a(number_of_point,dimen);

%% 进行计算计算前的最后一部分准备
%定义临界条件
f_critical=basic.f_critical;
s_critical=basic.s_critical;
%force_pd存储每个物质点总pd力
force_pd=zeros(number_of_point,2);
force_pd_old=zeros(number_of_point,2);
%vmiddle和vmiddle_old存储计算过程中的中间量
velocity_middle=zeros(number_of_point,2);
velocity_middle_old=zeros(number_of_point,2);
%先计算一次原始的形状张量并得到其逆矩阵
[tensor_shape_kb]=calculate_shape_tensor(mu,coordinate_x,list_of_near_point,displacement,number_of_point,dV,delta,dimen,0);
[inv_tensor_shape_kb]=calculate_inv_tensor(tensor_shape_kb,number_of_point,dimen);
%标记是否需要零能模式，0表示需要零能模式，不附加力态，1表示附加力态消除零能模式
flag_zero=1;


%% 开始时域积分
%记录图片编号
imt=1;
%计算阻尼系数，cn1为分子，cn2为分母，cn为阻尼系数
cn1=0;
cn2=0;
cn=0;
for it=1:nt
    %每个时间步先计算变形后的形状张量
    [tensor_shape_ka]=calculate_shape_tensor(mu,coordinate_x,list_of_near_point,displacement,number_of_point,dV,delta,dimen,1);
    %得到变形后的形状张量后计算变形梯度
    [deformation_grad]=calculate_deformation_grad(tensor_shape_ka,inv_tensor_shape_kb,dimen,number_of_point);
    %得到变形梯度后计算格林应变张量
    [tensor_strain] = calculate_strain_tensor(deformation_grad,number_of_point,dimen);
    %通过材料的本构计算柯西应力张量
    [tensor_stress_cauthy]=calculate_cauthy_stress_tensor(C,tensor_strain,number_of_point,dimen);
    %以伸长率为准则来破坏键，并计算损伤
    %[mu,damage]=judge_damage(coordinate_x,displacement,number_of_point,list_of_near_point,mu,damage,s_critical,dV);
    %通过柯西应力张量，变形梯度计算第一类PK应力张量
    [tensor_stress_PKI]=calculate_PKI_stress_tensor(tensor_stress_cauthy,deformation_grad,number_of_point,dimen);
    %通过上述计算得到的各种量，计算PD力
    [force_pd]=calculate_force_pd(force_body,mu,tensor_stress_PKI,inv_tensor_shape_kb,number_of_point,coordinate_x,list_of_near_point,dimen,delta,dV);
    %计算需要添加的力态（第一类沙漏力）
    %[force_fsand]=calculate_force_sand(mu,list_of_near_point,number_of_point,dimen,coordinate_x,displacement,E,dx,dV,delta);
    %计算需要添加的力态（参考DOI 10.1007/s10704-016-0126-6中的18）
    [force_fsand]=calculate_force_sand(c_pd,deformation_grad,mu,list_of_near_point,number_of_point,dimen,coordinate_x,displacement,dV,delta);
    %根据是否需要零能模式选择是否添加沙漏力，计算最终的PD力
    force_pd=force_pd+flag_zero*force_fsand;

    
    %进行自适应动态松弛法的计算
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
    %得到速度
    for i=1:number_of_point
        if it==1
            velocity_middle(i,1)=dt/density(i,1)*(force_pd(i,1))/2;
            velocity_middle(i,2)=dt/density(i,2)*(force_pd(i,2))/2;
        else
            velocity_middle(i,1)=((2-cn*dt)*velocity_middle_old(i,1)+2*dt/density(i,1)*(force_pd(i,1)))/(2+cn*dt);
            velocity_middle(i,2)=((2-cn*dt)*velocity_middle_old(i,2)+2*dt/density(i,2)*(force_pd(i,2)))/(2+cn*dt);
        end
    end
    %时域积分
    for i=1:number_of_point
        velocity(i,:)=0.5*(velocity_middle(i,:)+velocity_middle_old(i,:));
        displacement(i,:)=displacement(i,:)+velocity_middle(i,:)*dt;
        velocity_middle_old(i,:)=velocity_middle(i,:);
        force_pd_old(i,:)=force_pd(i,:);        
    end
    
    %计算中画图
    if rem(it,25)==0
        scale = 5;
        subplot(1,2,1);
        scatter(coordinate_x(:,1)+scale*displacement(:,1),coordinate_x(:,2)+scale*displacement(:,1),5, displacement(:,1), "filled");
        grid on;
        xlabel("x/m");ylabel("y/m");
        title("iteration=",num2str(it));
        colormap jet;
        colorbar;
        %caxis([0 1]);
        axis equal;
        axis([-0.06 0.06 -0.06 0.06]);
        subplot(1,2,2);
        scatter(coordinate_x(:,1)+scale*displacement(:,1),coordinate_x(:,2)+scale*displacement(:,1),5, tensor_stress_cauthy(:,1), "filled");
        grid on;
        xlabel("x/m");ylabel("y/m");
        title("iteration=",num2str(it));
        colormap jet;
        colorbar;
        %caxis([0 1]);
        axis equal;
        axis([-0.06 0.06 -0.06 0.06]);
        set(gcf,'Position',[350,100,900,600]);
        fig=gcf;
        frame=getframe(fig);
        im{imt} = frame2im(frame);
        imt=imt+1;
        picname="2d_plane_with_hole_nosb_fail_control";
        picname=strcat(picname,'-');
        picname=strcat(picname,num2str(it));
        picname=strcat(picname,".jpg");
        saveas(gcf,picname,'jpg');
    end
    fprintf("已完成%d/%d\n", it, nt);    
end

%% 导出GIF图
filename = 'D:\code\matlab\bin\math theory\PD\state_base\2d_plane_with_hole_left&right_no_fail_nosb\2d_plane_with_hole_nosb_fail_control.gif';
for idx = 1:imt
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end