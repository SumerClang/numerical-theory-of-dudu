function [mu,damage] = judge_damage(number_of_point,list_of_near_point,mu,damage,tensor_stress_cauthy,f_critical,dV)
%根据应力判断键是否断开，并计算损伤
%此时的柯西应力存储形式为[xx xy xy yy];
for ii=1:number_of_point
    %计算局部损伤
    sum1=0;
    sum2=0;
    %首先计算点i的主应力
    %先存储i点的柯西应力，此时的存储形式为[xx yy xy]
    sigmai=[tensor_stress_cauthy(ii,1) tensor_stress_cauthy(ii,4) tensor_stress_cauthy(ii,2)];
    %计算i点的第一主应力
    sigmait=max([(sigmai(1)+sigmai(2))/2+sqrt(((sigmai(1)-sigmai(2))/2)^2+sigmai(3)^2) (sigmai(1)+sigmai(2))/2-sqrt(((sigmai(1)-sigmai(2))/2)^2+sigmai(3)^2)]);    
    for jj=1:28
    %记录近场域内点的索引
    nn=list_of_near_point(ii,jj);
    %只有存取的点不是0的时候才进入计算
        if nn~=0
            %先存储j点的柯西应力，此时的存储形式为[xx yy xy]
            sigmaj=[tensor_stress_cauthy(nn,1) tensor_stress_cauthy(nn,4) tensor_stress_cauthy(nn,2)];
            %计算j点的第一主应力
            sigmajt=max([(sigmaj(1)+sigmaj(2))/2+sqrt(((sigmaj(1)-sigmaj(2))/2)^2+sigmaj(3)^2) (sigmaj(1)+sigmaj(2))/2-sqrt(((sigmaj(1)-sigmaj(2))/2)^2+sigmaj(3)^2)]);
            sigmat=(sigmait+sigmajt)/2;
            %当两个点的最大主应力的平均值大于临界压力时，断开键
            if sigmat>f_critical
%               mu(ii,jj)=0; 
            end
            sum1=sum1+mu(ii,jj)*dV;
            sum2=sum2+1*dV; 
        end
    end
    damage(ii)=1-sum1/sum2;
end
fprintf("已成功判断键的情况并计算损伤\n");
end