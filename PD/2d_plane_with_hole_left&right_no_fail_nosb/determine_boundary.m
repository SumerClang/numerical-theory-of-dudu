function [coordinate_x,number_of_point,force_body,node_left,node_right] = determine_boundary(coordinate_x,number_of_point,basic,dimen)
%根据basic information中的条件确定算例中的边界条件
%若存在位移边界，增加虚拟边界，若存在体力，设置force_body
%先判断位移边界
if isfield(basic,'uva')

    
end
%判断力边界
%存在体力条件时force_body对应处体力不为0
force_body=zeros(number_of_point,dimen);
node_left=zeros(200,1);
node_right=zeros(200,1);
if isfield(basic,'force_p')
    il=1;
    ir=1;
    for i=1:number_of_point
        %左边界
       if coordinate_x(i,1)==-0.04975
           force_body(i,1)=-4e10;
           node_left(il)=i;
           il=il+1;
       end
        %右边界
        if coordinate_x(i,1)==0.04975
            force_body(i,1)=4e10;
            node_right(ir)=i;
            ir=ir+1;
        end
    end    
end
fprintf("已成功确定边界条件\n");
end

