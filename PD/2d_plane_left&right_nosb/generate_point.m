function [coordinate_x,number_of_point] = generate_point(geometry,basic)
%根据geometry，basic，material中的信息生成质点坐标
%coordinate_x为质点坐标，number_of_point为质点数量
Length=geometry.Length;
Width=geometry.Width;
Thick=geometry.Thick;
dx=geometry.dx;
dy=geometry.dy;
nx=Length/dx;
ny=Width/dy;
dimen=basic.dimen;
%tx临时存储坐标
%第一列存储x坐标，第二列存储y坐标，第三列存储是否需要去掉该点,第三列为1时代表在圆外
number_of_point=nx*ny;
coordinate_x=zeros(number_of_point,dimen);
int=0;
if dimen == 2
    for i = 1:ny
        for j = 1:nx
            int = int+1;
            coordinate_x(int,1) = -0.5*Length+0.5*dx+(j-1)*dx;
            coordinate_x(int,2) = -0.5*Width+0.5*dx+(i-1)*dx;
        end
    end
elseif dimen == 3
    for i = 1:nz
        for j = 1:ny
            for k =1:nz
                int = int+1;
                coordinate_x(int,1) = -0.5*Length+0.5*dx+(k-1)*dx;
                coordinate_x(int,2) = -0.5*Width+0.5*dx+(j-1)*dx;
                coordinate_x(int,3) = -0.5*Thick+0.5*dx+(i-1)*dx;
            end
        end
    end
end
fprintf("已成功生成质点坐标\n");
end

