function [tensor_stress_PKI] = calculate_PKI_stress_tensor(tensor_stress_cauthy,deformation_grad,number_of_point,dimen)
%通过柯西应变张量和变形梯度求出第一类PK应力张量
%前两者初始以1*4的[xx xy xy yy]的形式存储，需要先重构
%最终PKI也以这种形式存储
tensor_stress_PKI=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %先将其重构为2*2的矩阵形式
    tensor_stress_cauthy_it=reshape(tensor_stress_cauthy(ii,:),dimen,dimen);
    deformation_grad_it=reshape(deformation_grad(ii,:),dimen,dimen);
    Jt=det(deformation_grad_it);
    %计算得到2*2的第一类PK应力张量，并将其以1*4[xx xy xy yy]的形式存储
    tensor_stress_PKI_it=Jt*tensor_stress_cauthy_it*(inv(deformation_grad_it)');
    tensor_stress_PKI(ii,:)=[tensor_stress_PKI_it(1,1) tensor_stress_PKI_it(1,2) tensor_stress_PKI_it(2,1) tensor_stress_PKI_it(2,2)];
end
fprintf("已成功计算PKI应力张量\n"); 
end

