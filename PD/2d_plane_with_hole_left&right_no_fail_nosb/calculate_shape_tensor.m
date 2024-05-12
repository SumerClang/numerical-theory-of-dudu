function [tensor_shape_k] = calculate_shape_tensor(mu,coordinate_x,list_of_near_point,displacement,number_of_point,dV,delta,dimen,flagk)
%计算并存储各个质点的形状张量,flagk为0表示计算初始的形状张量，flag为1表示计算变形后的形状张量
%b=before,a=after
tensor_shape_k=zeros(number_of_point,dimen*dimen);
for ii=1:number_of_point
    %临时的形状张量置0
    tensor_shape_kt=zeros(dimen,dimen);
    for jj=1:28
        nn=list_of_near_point(ii,jj);
         if nn~=0
                %临时存储坐标
                xxi=coordinate_x(ii,1);xyi=coordinate_x(ii,2);
                xxj=coordinate_x(nn,1);xyj=coordinate_x(nn,2);
                uxi=displacement(ii,1);uyi=displacement(ii,2);
                uxj=displacement(nn,1);uyj=displacement(nn,2);
                kexi=[xxj-xxi xyj-xyi];yita=[uxj-uxi uyj-uyi];
                Y=kexi+yita;X=kexi;
                %计算影响函数，有两种形式，一种是omega=exp(-norm(kexi)*norm(kexi)/delta/delta)，一种是omega=delta/norm(kexi)
                omega=exp(-norm(kexi)*norm(kexi)/delta^2);
                %omega=delta/norm(kexi);
                %求形状张量
                if mu(ii,jj)~=0
                    if flagk==0
                        tensor_shape_kt=tensor_shape_kt+omega*[X(1)*X(1) X(1)*X(2);X(2)*X(1) X(2)*X(2)]*dV;
                    elseif flagk==1
                        tensor_shape_kt=tensor_shape_kt+omega*[Y(1)*X(1) Y(1)*X(2);Y(2)*X(1) Y(2)*X(2)]*dV;                    
                    end
                end
         end
        %将2*2矩阵转换为1*4进行存储
        tensor_shape_k(ii,:)=[tensor_shape_kt(1,1) tensor_shape_kt(1,2) tensor_shape_kt(2,1) tensor_shape_kt(2,2)];        
    end    
end

if flagk==0
    fprintf("已成功计算初始形状张量\n");  
elseif flagk==1
    fprintf("已成功计算变形后形状张量\n");
end
end

