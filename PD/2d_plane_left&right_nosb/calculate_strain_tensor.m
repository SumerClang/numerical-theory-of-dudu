function [tensor_strain] = calculate_strain_tensor(deformation_grad,number_of_point,dimen)
%通过变形梯度计算应变张量
%应变张量以1*4的形式存储，存储顺序为xx，xy，xy，yy
tensor_strain=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %将1*4的变形梯度转化为2*2
    deformation_grad_it=reshape(deformation_grad(ii,:),dimen,dimen);
    %此时得到的是2*2的应变张量，存储顺序为[xx xy;xy yy]
    tensor_strain_it=0.5*(deformation_grad_it'*deformation_grad_it-eye(dimen));
    %将2*2的应变张量转化为1*4的形式存储
    tensor_strain(ii,:)=[tensor_strain_it(1,1) tensor_strain_it(1,2) tensor_strain_it(2,1) tensor_strain_it(2,2)];   
end
    fprintf("已成功计算应变张量\n");
end

