function [density] = generate_massvector(number_of_point,dt,delta,Thick,dx,safe,dimen,c_pd)
%该函数用于生成自适应动态松弛中的质量密度矢量
    density=zeros(number_of_point,dimen);
    for i=1:number_of_point
        density(i,1)=0.25*dt*dt*(pi*delta^2*Thick)*c_pd/dx*safe;
        density(i,2)=0.25*dt*dt*(pi*delta^2*Thick)*c_pd/dx*safe;
    end
    fprintf("已成功生成自适应动态松弛中的质量密度矢量\n");
end

