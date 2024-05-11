function [tensor_stress_cauthy] = calculate_cauthy_stress_tensor(C,tensor_strain,number_of_point,dimen)
%通过格林应变张量和本构关系计算柯西应力张量
%最终柯西应力张量以[xx xy xy yy]的形式存储
tensor_stress_cauthy=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %计算各个节点时首先需要将排列为[xx xy;xy yy]的应变写成[xx yy xy]的形式
    tensor_strain_it=[tensor_strain(ii,1);tensor_strain(ii,4);tensor_strain(ii,2)];
    tensor_stress_cauthy_it=C*tensor_strain_it;
    %上一步计算得到的柯西应力以[xx yy xy]的形式存在
    %最终以[xx xy xy yy]的形式存储
    tensor_stress_cauthy(ii,:)=[tensor_stress_cauthy_it(1) tensor_stress_cauthy_it(3) tensor_stress_cauthy_it(3) tensor_stress_cauthy_it(2)];
end
fprintf("已成功计算柯西应力张量\n");
end

