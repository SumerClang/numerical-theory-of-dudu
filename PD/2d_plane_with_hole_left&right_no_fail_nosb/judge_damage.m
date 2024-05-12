function [mu,damage] = judge_damage(coordinate_x,displacement,number_of_point,list_of_near_point,mu,damage,s_critical,dV)
%根据应力判断键是否断开，并计算损伤
%此时的柯西应力存储形式为[xx xy xy yy];
coordinate_y=coordinate_x+displacement;
for ii=1:number_of_point
    %计算局部损伤
    sum1=0;
    sum2=0;
    for jj=1:28
    %记录近场域内点的索引
    nn=list_of_near_point(ii,jj);
    %只有存取的点不是0的时候才进入计算
        if nn~=0
            xxi=coordinate_x(ii,1);xyi=coordinate_x(ii,2);
            xxj=coordinate_x(nn,1);xyj=coordinate_x(nn,2);
            yxi=coordinate_y(ii,1);yyi=coordinate_y(ii,2);
            yxj=coordinate_y(nn,1);yyj=coordinate_y(nn,2);
            %计算两个点之间的距离，以及变形后两个点之间的距离
            distancex=sqrt((xxj-xxi)*(xxj-xxi)+(xyj-xyi)*(xyj-xyi));
            distancey=sqrt((yxj-yxi)*(yxj-yxi)+(yyj-yyi)*(yyj-yyi));
            %计算伸长率
            st=(distancey-distancex)/distancex;
            %当两个点之间的伸长率大于临界伸长率时，断开键
            if st>s_critical
               mu(ii,jj)=0; 
            end
            sum1=sum1+mu(ii,jj)*dV;
            sum2=sum2+1*dV; 
        end
    end
    damage(ii)=1-sum1/sum2;
end
fprintf("已成功判断键的情况并计算损伤\n");
end