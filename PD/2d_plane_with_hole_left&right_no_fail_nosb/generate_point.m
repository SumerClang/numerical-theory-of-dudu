function [coordinate_x,number_of_point] = generate_point(geometry,basic)
%根据geometry，basic，material中的信息生成质点坐标
%coordinate_x为质点坐标，number_of_point为质点数量
Length=geometry.Length;
Width=geometry.Width;
Thick=geometry.Thick;
R=geometry.R;
dx=geometry.dx;
dy=geometry.dy;
nx=Length/dx;
ny=Width/dy;
dimen=basic.dimen;
%tx临时存储坐标
%第一列存储x坐标，第二列存储y坐标，第三列存储是否需要去掉该点,第三列为1时代表在圆外
number_of_point=nx*ny;
coordinate_tx=zeros(number_of_point,dimen+1);
t=1;
%生成tx中各个点的坐标
coordinate_tx(:,3)=ones(nx*ny,1);
for i=1:ny
    for j=1:nx
        coordinate_tx(t,1) = -0.5*Length+0.5*dx+(j-1)*dx;
        coordinate_tx(t,2) = -0.5*Width+0.5*dx+(i-1)*dx;
        t=t+1;
    end
end
%判断该点到原点的距离，距离大于R则需要除去
for i=1:number_of_point
    tx=coordinate_tx(i,1);
    ty=coordinate_tx(i,2);
    tr=sqrt(tx*tx+ty*ty);
    if tr<R
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
fprintf("已成功生成质点坐标\n");
end

