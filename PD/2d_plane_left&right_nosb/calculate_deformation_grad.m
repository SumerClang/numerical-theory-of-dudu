function [deformation_grad] = calculate_deformation_grad(tensor_shape_ka,inv_tensor_shape_kb,dimen,number_of_point)
%此函数用于计算各个物质点的变形梯度，并将各个点的2*2的变形梯度转化为1*4存起来
%计算需要变形前的形状张量的逆矩阵和变形后的形状张量
%ka为变形后的，kb为变形前的
deformation_grad=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %先将形状张量和原形状张量的逆从1*4转化为2*2的形式
    tensor_shapa_ka_it=reshape(tensor_shape_ka(ii,:),dimen,dimen);
    inv_tensor_shapa_kb_it=reshape(inv_tensor_shape_kb(ii,:),dimen,dimen);
    %计算得到每个点的变形梯度（2*2形式），然后转化为1*4存储起来
    deformation_grad_it=tensor_shapa_ka_it*inv_tensor_shapa_kb_it;
    deformation_grad(ii,:)=[deformation_grad_it(1,1) deformation_grad_it(1,2) deformation_grad_it(2,1) deformation_grad_it(2,2)];    
end
fprintf("已成功计算变形梯度\n");
end

