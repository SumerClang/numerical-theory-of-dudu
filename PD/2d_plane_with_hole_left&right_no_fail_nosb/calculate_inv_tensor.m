function [inv_tensor] = calculate_inv_tensor(tensor,number_of_point,dimen)
%计算并输出非正常张量排列方式的逆矩阵
inv_tensor=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %将tensor中的1*4矩阵转化为2*2再做运算，最后转化为1*4存起来
    tensor_it=reshape(tensor(ii,:),dimen,dimen);
    invtensor_ti=inv(tensor_it);
    inv_tensor(ii,:)=[invtensor_ti(1,1) invtensor_ti(1,2) invtensor_ti(2,1) invtensor_ti(2,2)];    
end
fprintf("已成功计算逆矩阵\n");
end

