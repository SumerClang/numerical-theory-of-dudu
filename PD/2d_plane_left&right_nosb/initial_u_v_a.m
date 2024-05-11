function [displacement,velocity,acceleration] = initial_u_v_a(number_of_point,dimen)
%INITIAL_U_V_A 根据具体情况对位移，速度，加速度进行初始化
displacement=zeros(number_of_point,dimen);
velocity=zeros(number_of_point,dimen);
acceleration=zeros(number_of_point,dimen);
fprintf("已成功初始化u,v,a\n");
end

