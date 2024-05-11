function [list_of_near_point,mu,damage] = determine_near_point(coordinate_x,number_of_point,delta)
%判断质点近场域内的其他点，并初步设置相互作用
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
fprintf("已成功判断近场域内质点并初步设置了连接情况\n");
end

