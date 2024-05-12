function [force_fsand] = calculate_force_sand(c_pd,deformation_grad,mu,list_of_near_point,number_of_point,dimen,coordinate_x,displacement,dV,delta)
%考虑损伤，体力等，计算总的沙漏力
for ii=1:number_of_point
    force_sand_it=zeros(1,dimen);
    for jj=1:28
        nn=list_of_near_point(ii,jj);
        if nn~=0
            %先计算影响函数
            %临时存储坐标
            xxi=coordinate_x(ii,1);xyi=coordinate_x(ii,2);
            xxj=coordinate_x(nn,1);xyj=coordinate_x(nn,2);
            uxi=displacement(ii,1);uyi=displacement(ii,2);
            uxj=displacement(nn,1);uyj=displacement(nn,2);
            kexi=[xxj-xxi xyj-xyi];yita=[uxj-uxi uyj-uyi];
            %计算影响函数，有两种形式，一种是omega=exp(-norm(kexi)*norm(kexi)/delta/delta)，一种是omega=delta/norm(kexi)
            %omega=delta/norm(kexi);
            omega=exp(-norm(kexi)*norm(kexi)/delta/delta);
            Yi=kexi+yita;Yj=-Yi;
            Fi=[deformation_grad(ii,1) deformation_grad(ii,2);deformation_grad(ii,3) deformation_grad(ii,4)];
            Fj=[deformation_grad(jj,1) deformation_grad(jj,2);deformation_grad(jj,3) deformation_grad(jj,4)];
            Zi=Yi-(Fi*(kexi'))';Zj=Yj-(Fj*((-kexi)'))';
            C0=c_pd*[kexi(1)*kexi(1) kexi(1)*kexi(2);kexi(2)*kexi(1) kexi(2)*kexi(2)];
            %当键未断开时才进行计算
            if mu(ii,jj)~=0              
                force_sand_it=force_sand_it+((omega*C0*(Zi')*dV-omega*C0*(Zj')*dV)');                
            end
        end       
    end
    force_fsand(ii,:)=force_sand_it;
end
    fprintf("已成功计算沙漏力\n");
end