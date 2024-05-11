function [basic,geometry,material] = basic_information()
% 此函数用于给定算例的几何，材料等基本信息
%% 先给定几何信息，包括圆盘半径，厚度，物质点密度
geometry.Length=0.1;
geometry.Width=0.1;
geometry.Thick=0.5/1000;
geometry.dx=0.5/1000;
geometry.dy=0.5/1000;
geometry.A=geometry.dy*geometry.dx;
geometry.dV=geometry.dx*geometry.dy*geometry.Thick;

%% 给定算例材料信息，包括热膨胀系数，杨氏模量，泊松比，热膨胀系数，体积模量，剪切模量，密度，应变转化为柯西应力张量所需的材料矩阵
material.alpha=23e-6;
material.E=200e9;
material.Poisson=0.2;
material.alpha=23e-6;
material.K=material.E/(2*(1-material.Poisson));
material.G=material.E/(2*(1+material.Poisson));
material.rou=1900;
material.C=[material.K+material.G material.K-material.G 0;material.K-material.G material.K+material.G 0;0 0 material.G];

%% 给定算例基本信息，包括是否使用动态松弛，时间步长，时间步数，近场域大小，算例的维度，存在何种边界条件
basic.adr=1;%使用动态松弛
basic.dt=1;%使用动态松弛时时间步长为1
basic.nt=1000;
basic.delta=1.05*geometry.dx;
basic.dimen=2;
basic.force_p=[1 0];%力边界条件
%basic.uva=[0 1];%位移边界
fprintf("已成功设置基本信息\n");
end

