function [coordinate_x,number_of_point,force_body,number_of_virtual_boundary,node_left,node_right] = determine_boundary(coordinate_x,number_of_point,basic,dimen,dx,dy)
%根据basic information中的条件确定算例中的边界条件
%若存在位移边界，增加虚拟边界，若存在体力，设置force_body
%先判断位移边界
%本算例中上下最上和最下的四个节点存在位移，增加三层虚拟边界
t=number_of_point;
if isfield(basic,'uva')
    for i=1:number_of_point
        %当y坐标为0.025375时，在其上方施加三层虚拟边界
       if coordinate_x(i,2)==0.025375
           for j=1:3
                coordinate_x(t+j,1)=coordinate_x(i,1)+dx*0;
                coordinate_x(t+j,2)=coordinate_x(i,2)+j*dy;
           end
           t=t+3;
       end
        %当y坐标为-0.025375时，在其下方施加三层虚拟边界
       if coordinate_x(i,2)==-0.025375
           for j=1:3
                coordinate_x(t+j,1)=coordinate_x(i,1)+dx*0;
                coordinate_x(t+j,2)=coordinate_x(i,2)-j*dy;
           end
           t=t+3;
       end        
    end    
end
number_of_virtual_boundary=t-number_of_point;
number_of_point=t;
%判断需要施加向上和向下的位移的节点
node_down=zeros(number_of_virtual_boundary/6,1);
node_up=zeros(number_of_virtual_boundary/6,1);
tu=1;td=1;
for i=1:number_of_point
   if coordinate_x(i,2)>0.025875 
    node_down(td,1)=i;
    td=td+1;
   end    
   if coordinate_x(i,2)<-0.025875
    node_up(tu,1)=i;
    tu=tu+1;
   end
end


%判断力边界
%存在体力条件时force_body对应处体力不为0
force_body=zeros(number_of_point,dimen);
if isfield(basic,'force_p')
    %此算例中最开始的200个和最后的200个受到拉力
    node_left=linspace(1,39801,200);
    node_right=linspace(200,40000,200);
    %左边向左，右边向右
    force_body(node_left,1)=-4e11*ones(200,1);
    force_body(node_right,1)=4e11*ones(200,1);
end
fprintf("已成功确定边界条件\n");
end

