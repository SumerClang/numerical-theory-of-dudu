function [force_pd] = calculate_force_pd(force_body,mu,tensor_stress_PKI,inv_tensor_shape_kb,number_of_point,coordinate_x,list_of_near_point,dimen,delta,dV)
%考虑损伤，体力等，计算总的PD力
%PD力先置0
force_pd=zeros(number_of_point,dimen);
for ii=1:number_of_point
    force_pd_it=zeros(1,dimen);
    for jj=1:28
        nn=list_of_near_point(ii,jj);
        if nn~=0
            %先计算影响函数
            %临时存储坐标
            xxi=coordinate_x(ii,1);xyi=coordinate_x(ii,2);
            xxj=coordinate_x(nn,1);xyj=coordinate_x(nn,2);
            kexi=[xxj-xxi xyj-xyi];
            %计算影响函数，有两种形式，一种是omega=exp(-norm(kexi)*norm(kexi)/delta/delta)，一种是omega=delta/norm(kexi)
            %omega=delta/norm(kexi);
            omega=exp(-norm(kexi)*norm(kexi)/delta/delta);
            %当键未断开时才进行计算
            if mu(ii,jj)~=0
                %首先将原来存储的1*4转化为2*2
                inv_tensor_kb_it=reshape(inv_tensor_shape_kb(ii,:),dimen,dimen);
                inv_tensor_kb_jt=reshape(inv_tensor_shape_kb(nn,:),dimen,dimen);
                tensor_stress_PKI_it=reshape(tensor_stress_PKI(ii,:),dimen,dimen);
                tensor_stress_PKI_jt=reshape(tensor_stress_PKI(nn,:),dimen,dimen);
                force_pd_it=force_pd_it+(omega*tensor_stress_PKI_it*inv_tensor_kb_it*(kexi')*dV+omega*tensor_stress_PKI_jt*inv_tensor_kb_jt*(kexi')*dV)';                
            end
        end       
    end
    force_pd(ii,:)=force_pd_it;
end
force_pd=force_pd+force_body;
fprintf("已成功计算不含沙漏力的PD力\n");
end